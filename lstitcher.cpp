#include "lstitcher.h"
//#include "modules/lightmapper_cpu/lightmapper_cpu.h"

using namespace LM;

void Stitcher::StitchObjectSeams(const MeshInstance &p_mi, LightImage<FColor> &r_image, float distance_threshold, float normal_threshold, bool p_visualize_seams) {

	// important - the distance threshold is assumed to be in world space,
	// but we are using the model in model space here. So we must apply a scaling
	// to roughly get them to match up. Note this assumes uniform scaling. With non-uniform
	// scaling your mileage may vary...
	m_fDistanceThreshold = distance_threshold;
	Vector3 scale3 = p_mi.get_global_transform().basis.get_scale();
	float scale = (scale3.x + scale3.y + scale3.z) / 3.0f;

	// prevent divide by zero
	if (scale < 0.00001f)
		return;

	m_fDistanceThreshold *= 1.0f / scale;

	// position epsilon must take account scale
	m_fPositionEpsilon = 0.001;
	m_fPositionEpsilon *= 1.0f / scale;
	m_fPositionEpsilon *= m_fPositionEpsilon;

	// normals should be same whatever the scale
	m_fNormalThreshold = normal_threshold;

	m_bVisualizeSeams = p_visualize_seams;

	// seam optimizer wants unindexed lists of tri vertex positions and uvs.
	// some godot jiggery pokery to get the mesh verts in local space
	Ref<Mesh> rmesh = p_mi.get_mesh();

	if (rmesh->get_surface_count() == 0)
		return;

	Array arrays = rmesh->surface_get_arrays(0);
	PoolVector<Vector3> p_vertices = arrays[VS::ARRAY_VERTEX];
	PoolVector<Vector2> p_uvs = arrays[VS::ARRAY_TEX_UV2];
	PoolVector<Vector3> p_normals = arrays[VS::ARRAY_NORMAL];
	PoolVector<int> p_indices = arrays[VS::ARRAY_INDEX];

	if (!p_vertices.size())
		return;
	if (!p_uvs.size())
		return;
	if (!p_normals.size())
		return;

	int num_inds = p_indices.size();

	LocalVector<Vector3> verts;
	verts.resize(num_inds);
	LocalVector<Vector2> uvs;
	uvs.resize(num_inds);
	LocalVector<Vector3> norms;
	norms.resize(num_inds);

	for (int n = 0; n < num_inds; n++) {
		int i = p_indices[n];
		verts[n] = p_vertices[i];
		uvs[n] = p_uvs[i];
		norms[n] = p_normals[i];
	}

	Vector3 *pImageData = (Vector3 *)r_image.Get(0);
	Vector2i lm_size = Vector2i(r_image.GetWidth(), r_image.GetHeight());

	//////////////////////////////////

	LocalVector<UVSeam> seams;

	_compute_seams(verts, uvs, norms, lm_size, seams);

	if (!seams.size())
		return;

	print_line("seam stitching " + p_mi.get_name() + ", " + itos(seams.size()) + " seams.");

	_fix_seams(seams, pImageData, lm_size);
}

void Stitcher::_compute_seams(const Vector<Vector3> &points, const Vector<Vector2> &uv2s, const Vector<Vector3> &normals, Vector2i lm_size, LocalVector<UVSeam> &r_seams) {
	float max_uv_distance = 1.0f / MAX(lm_size.x, lm_size.y);
	max_uv_distance *= max_uv_distance; // We use distance_to_squared(), so wee need to square the max distance as well

	//	float max_pos_distance = 0.00025f;
	//	float max_normal_distance = 0.05f;

	// the checks are against square distance
	float max_pos_distance = m_fDistanceThreshold * m_fDistanceThreshold;

	// normals are now measured by dot product
	float min_dot_threshold = Math::cos(Math::deg2rad(m_fNormalThreshold));

	//	float max_normal_distance = m_fNormalThreshold * m_fNormalThreshold;

	//	const Vector<Vector3> &points = p_mesh.data.points;
	//	const Vector<Vector2> &uv2s = p_mesh.data.uv2;
	//	const Vector<Vector3> &normals = p_mesh.data.normal;

	LocalVector<SeamEdge> edges;
	edges.resize(points.size()); // One edge per vertex

	for (int i = 0; i < points.size(); i += 3) {
		Vector3 triangle_vtxs[3] = { points[i + 0], points[i + 1], points[i + 2] };
		Vector2 triangle_uvs[3] = { uv2s[i + 0], uv2s[i + 1], uv2s[i + 2] };
		Vector3 triangle_normals[3] = { normals[i + 0], normals[i + 1], normals[i + 2] };

		for (int k = 0; k < 3; k++) {
			int idx[2];
			idx[0] = k;
			idx[1] = (k + 1) % 3;

			if (triangle_vtxs[idx[1]] < triangle_vtxs[idx[0]]) {
				SWAP(idx[0], idx[1]);
			}

			SeamEdge e;
			for (int l = 0; l < 2; ++l) {
				e.pos[l] = triangle_vtxs[idx[l]];
				e.uv[l] = triangle_uvs[idx[l]];
				e.normal[l] = triangle_normals[idx[l]];
			}
			edges[i + k] = e;
		}
	}

	edges.sort();

	for (unsigned int j = 0; j < edges.size(); j++) {
		const SeamEdge &edge0 = edges[j];

		// should this be editable? not yet implemented
		const float uv_epsilon = 0.001 * 0.001; // 0.001

		// position epsilon must take account scale
		//		const float pos_epsilon = 0.001 * 0.001; // 0.001
		const float pos_epsilon = m_fPositionEpsilon;

		// get rid of zero area edges
		if (edge0.uv[0].distance_squared_to(edge0.uv[1]) < uv_epsilon) {
			continue;
		}

		if (edge0.pos[0].distance_squared_to(edge0.pos[1]) < pos_epsilon) {
			continue;
		}

		for (unsigned int k = j + 1; k < edges.size() && edges[k].pos[0].x < (edge0.pos[0].x + max_pos_distance * 1.1f); k++) {
			const SeamEdge &edge1 = edges[k];

			// reject if the second edge is too small
			if (edge1.uv[0].distance_squared_to(edge1.uv[1]) < uv_epsilon) {
				continue;
			}

			if (edge1.pos[0].distance_squared_to(edge1.pos[1]) < pos_epsilon) {
				continue;
			}

			///////////////////////////////

			float uv_dist_0 = edge0.uv[0].distance_squared_to(edge1.uv[0]);
			float uv_dist_1 = edge0.uv[1].distance_squared_to(edge1.uv[1]);

			// if the edges are not seams (i.e. internal)
			if (uv_dist_0 < max_uv_distance && uv_dist_1 < max_uv_distance) {
				continue;
			}

			// if the edge positions are too far in space to be connected
			float dist_0 = edge0.pos[0].distance_squared_to(edge1.pos[0]);
			float dist_1 = edge0.pos[1].distance_squared_to(edge1.pos[1]);
			if (dist_0 > max_pos_distance || dist_1 > max_pos_distance) {
				continue;
			}

			float dot0 = edge0.normal[0].dot(edge1.normal[0]);
			float dot1 = edge0.normal[1].dot(edge1.normal[1]);

			if ((dot0 < min_dot_threshold) || (dot1 < min_dot_threshold)) {
				continue;
			}

			// old normal using distance
			//			if (edge0.normal[0].distance_squared_to(edge1.normal[0]) > max_normal_distance || edge0.normal[1].distance_squared_to(edge1.normal[1]) > max_normal_distance) {
			//				continue;
			//			}

			UVSeam s;
			s.edge0[0] = edge0.uv[0];
			s.edge0[1] = edge0.uv[1];
			s.edge1[0] = edge1.uv[0];
			s.edge1[1] = edge1.uv[1];
			r_seams.push_back(s);
		}
	}
}

void Stitcher::_fix_seams(const LocalVector<UVSeam> &p_seams, Vector3 *r_lightmap, Vector2i p_size) {
	LocalVector<Vector3> extra_buffer;
	extra_buffer.resize(p_size.x * p_size.y);

	copymem(extra_buffer.ptr(), r_lightmap, p_size.x * p_size.y * sizeof(Vector3));

	Vector3 *read_ptr = extra_buffer.ptr();
	Vector3 *write_ptr = r_lightmap;

	for (int i = 0; i < 5; i++) {
		for (unsigned int j = 0; j < p_seams.size(); j++) {
			_fix_seam(p_seams[j].edge0[0], p_seams[j].edge0[1], p_seams[j].edge1[0], p_seams[j].edge1[1], read_ptr, write_ptr, p_size);
			_fix_seam(p_seams[j].edge1[0], p_seams[j].edge1[1], p_seams[j].edge0[0], p_seams[j].edge0[1], read_ptr, write_ptr, p_size);
		}
		copymem(read_ptr, write_ptr, p_size.x * p_size.y * sizeof(Vector3));
	}
}

void Stitcher::_fix_seam(const Vector2 &p_pos0, const Vector2 &p_pos1, const Vector2 &p_uv0, const Vector2 &p_uv1, const Vector3 *p_read_buffer, Vector3 *r_write_buffer, const Vector2i &p_size) {
	Vector2 line[2];
	line[0] = p_pos0 * p_size;
	line[1] = p_pos1 * p_size;

	const Vector2i start_pixel = line[0].floor();
	const Vector2i end_pixel = line[1].floor();

	Vector2 seam_dir = (line[1] - line[0]).normalized();
	Vector2 t_delta = Vector2(1.0f / Math::abs(seam_dir.x), 1.0f / Math::abs(seam_dir.y));
	Vector2i step = Vector2(seam_dir.x > 0 ? 1 : (seam_dir.x < 0 ? -1 : 0), seam_dir.y > 0 ? 1 : (seam_dir.y < 0 ? -1 : 0));

	Vector2 t_next = Vector2(Math::fmod(line[0].x, 1.0f), Math::fmod(line[0].y, 1.0f));

	if (step.x == 1) {
		t_next.x = 1.0f - t_next.x;
	}

	if (step.y == 1) {
		t_next.y = 1.0f - t_next.y;
	}

	t_next.x /= Math::abs(seam_dir.x);
	t_next.y /= Math::abs(seam_dir.y);

	if (Math::is_nan(t_next.x)) {
		t_next.x = 1e20f;
	}

	if (Math::is_nan(t_next.y)) {
		t_next.y = 1e20f;
	}

	Vector2i pixel = start_pixel;
	Vector2 start_p = start_pixel;
	float line_length = line[0].distance_to(line[1]);

	if (line_length == 0.0f) {
		return;
	}

	while (start_p.distance_to(pixel) < line_length + 1.0f) {

		Vector2 current_point = Vector2(pixel) + Vector2(0.5f, 0.5f);
		current_point = Geometry::get_closest_point_to_segment_2d(current_point, line);
		float t = line[0].distance_to(current_point) / line_length;

		Vector2 current_uv = p_uv0 * (1.0 - t) + p_uv1 * t;
		Vector2i sampled_point = (current_uv * p_size).floor();

		Vector3 current_color = r_write_buffer[pixel.y * p_size.x + pixel.x];
		Vector3 sampled_color = p_read_buffer[sampled_point.y * p_size.x + sampled_point.x];

		if (m_bVisualizeSeams)
			r_write_buffer[pixel.y * p_size.x + pixel.x] = Vector3(1.0f, 0.0f, 0.0f);
		else
			r_write_buffer[pixel.y * p_size.x + pixel.x] = current_color * 0.6f + sampled_color * 0.4f;

		if (pixel == end_pixel) {
			break;
		}

		if (t_next.x < t_next.y) {
			pixel.x += step.x;
			t_next.x += t_delta.x;
		} else {
			pixel.y += step.y;
			t_next.y += t_delta.y;
		}
	}
}
