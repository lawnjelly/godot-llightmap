#include "llightmapper.h"
#include "core/os/threaded_array_processor.h"
#include "llightprobe.h"
#include "llightscene.h"
#include "lmerger.h"
#include "lscene_saver.h"
#include "lunmerger.h"
#include "scene/resources/packed_scene.h"

namespace LM {

bool LightMapper::uv_map_meshes(Spatial *pRoot) {
	if (m_Settings_UVFilename == "") {
		ShowWarning("UV Filename is not set, aborting.");
		return false;
	}

	CalculateQualityAdjustedSettings();

	bool replace_mesh_scene = false;

	if (!pRoot)
		return false;

	// can't do uv mapping if not tools build, as xatlas isn't compiled
	if (bake_begin_function) {
		bake_begin_function(4);
	}

	if (bake_step_function) {
		bake_step_function(0, String("Saving uvmap_backup.tscn"));
	}

	// first back up the existing meshes scene.
	SceneSaver saver;
	saver.SaveScene(pRoot, "res://uvmap_backup.tscn");

	if (bake_step_function) {
		bake_step_function(1, String("Merging to proxy"));
	}

	Merger m;
	MeshInstance *pMerged = m.Merge(pRoot, m_Settings_UVPadding);
	if (!pMerged) {
		if (bake_end_function) {
			bake_end_function();
		}
		return false;
	}

	// test save the merged mesh
	//saver.SaveScene(pMerged, "res://merged_test.tscn");

	if (bake_step_function) {
		bake_step_function(2, String("Unmerging"));
	}

	// unmerge
	UnMerger u;
	bool res = u.UnMerge(m, *pMerged);

	// for debug save the merged version
	saver.SaveScene(pMerged, "res://merged_proxy.tscn");

	pMerged->queue_delete();

	if (bake_step_function) {
		bake_step_function(2, String("Saving"));
	}

	// if we specified an output file, save
	if (m_Settings_UVFilename != "") {
		Node *pOrigOwner = pRoot->get_owner();

		if (!saver.SaveScene(pRoot, m_Settings_UVFilename, true)) {
			ShowWarning("Error saving UV mapped scene. Does the folder exist?\n\n" + m_Settings_UVFilename);
		}

		if (replace_mesh_scene) {
			// delete the orig
			Node *pParent = pRoot->get_parent();

			// rename the old scene to prevent naming conflict
			pRoot->set_name("ToBeDeleted");
			pRoot->queue_delete();

			// now load the new file
			//ResourceLoader::import(m_Settings_UVFilename);

			Ref<PackedScene> ps = ResourceLoader::load(m_Settings_UVFilename, "PackedScene");
			if (ps.is_null()) {
				if (bake_end_function) {
					bake_end_function();
				}
				return res;
			}

			Node *pFinalScene = ps->instance();
			if (pFinalScene) {
				pParent->add_child(pFinalScene);

				// set owners
				saver.SetOwnerRecursive(pFinalScene, pFinalScene);
				pFinalScene->set_owner(pOrigOwner);
			}

		} // if replace the scene
		else {
			// delete, to save confusion
			pRoot->queue_delete();
		}
	}

	if (bake_end_function) {
		bake_end_function();
	}

	return res;
}

bool LightMapper::lightmap_mesh(Spatial *pMeshesRoot, Spatial *pLR, Image *pIm_Lightmap, Image *pIm_AO, Image *pIm_Combined) {
	CalculateQualityAdjustedSettings();

	// get the output dimensions before starting, because we need this
	// to determine number of rays, and the progress range
	m_iWidth = pIm_Combined->get_width();
	m_iHeight = pIm_Combined->get_height();
	m_iNumRays = m_AdjustedSettings.m_Forward_NumRays;

	int nTexels = m_iWidth * m_iHeight;

	// set num rays depending on method
	if (m_Settings_Mode == LMMODE_FORWARD) {
		// the num rays / texel. This is per light!
		m_iNumRays *= nTexels;
	}

	// do twice to test SIMD
	uint32_t beforeA = OS::get_singleton()->get_ticks_msec();
	m_Scene.m_bUseSIMD = true;
	m_Scene.m_Tracer.m_bSIMD = true;

	bool res = LightmapMesh(pMeshesRoot, *pLR, *pIm_Lightmap, *pIm_AO, *pIm_Combined);

	uint32_t afterA = OS::get_singleton()->get_ticks_msec();
	print_line("Overall took : " + itos(afterA - beforeA));

	return res;
}

void LightMapper::Reset() {
	m_Lights.clear(true);
	m_Scene.Reset();
	RayBank_Reset();
	Base_Reset();
}

void LightMapper::Refresh_Process_State() {
	// defaults
	m_Logic_Process_Lightmap = true;
	m_Logic_Process_AO = true;
	m_Logic_Reserve_AO = true;
	m_Logic_Process_Probes = false;
	m_Logic_Output_Final = true;

	// process states
	switch (m_Settings_BakeMode) {
		case LMBAKEMODE_PROBES: {
			m_Logic_Process_Probes = true;
			m_Logic_Process_Lightmap = false;
			m_Logic_Process_AO = false;
			m_Logic_Output_Final = false;
		} break;
		case LMBAKEMODE_LIGHTMAP: {
			m_Logic_Process_AO = false;
			m_Logic_Reserve_AO = false;
		} break;
		case LMBAKEMODE_AO: {
			m_Logic_Process_Lightmap = false;
		} break;
		case LMBAKEMODE_MERGE: {
			m_Logic_Process_Lightmap = false;
			m_Logic_Process_AO = false;
		} break;
		case LMBAKEMODE_UVMAP: {
			m_Logic_Process_Lightmap = false;
			m_Logic_Process_AO = false;
			m_Logic_Reserve_AO = false;
			m_Logic_Output_Final = false;
		} break;
		default: {
		} break;
	}
}

bool LightMapper::LightmapMesh(Spatial *pMeshesRoot, const Spatial &light_root, Image &out_image_lightmap, Image &out_image_ao, Image &out_image_combined) {
	// print out settings
	print_line("Lightmap mesh");
	print_line("\tnum_directional_bounces " + itos(m_AdjustedSettings.m_NumDirectionalBounces));
	print_line("\tdirectional_bounce_power " + String(Variant(m_Settings_DirectionalBouncePower)));

	Refresh_Process_State();

	Reset();
	m_bCancel = false;

	if (m_iWidth <= 0)
		return false;
	if (m_iHeight <= 0)
		return false;

	// create stuff used by everything
	m_Image_L.Create(m_iWidth, m_iHeight);
	m_Image_L_mirror.Create(m_iWidth, m_iHeight);

	// whether we need storage
	if (m_Logic_Reserve_AO)
		m_Image_AO.Create(m_iWidth, m_iHeight);

	if (m_Settings_BakeMode != LMBAKEMODE_MERGE) {
		if (m_Logic_Process_AO)
			m_QMC.Create(m_AdjustedSettings.m_AO_Samples);

		uint32_t before, after;
		FindLights_Recursive(&light_root);
		print_line("Found " + itos(m_Lights.size()) + " lights.");

		m_Image_ID_p1.Create(m_iWidth, m_iHeight);
		//m_Image_ID2_p1.Create(m_iWidth, m_iHeight);

		if (m_Logic_Process_AO) {
			m_Image_TriIDs.Create(m_iWidth, m_iHeight);
			m_TriIDs.clear(true);
		}

		m_Image_Barycentric.Create(m_iWidth, m_iHeight);

		//m_Image_Cuts.Create(m_iWidth, m_iHeight);
		//m_CuttingTris.clear(true);

		print_line("Scene Create");
		before = OS::get_singleton()->get_ticks_msec();
		if (!m_Scene.Create(pMeshesRoot, m_iWidth, m_iHeight, m_Settings_VoxelDensity, m_AdjustedSettings.m_Max_Material_Size, m_AdjustedSettings.m_EmissionDensity))
			return false;

		PrepareLights();

		RayBank_Create();

		after = OS::get_singleton()->get_ticks_msec();
		print_line("SceneCreate took " + itos(after - before) + " ms");

		if (m_bCancel)
			return false;

		print_line("PrepareImageMaps");
		before = OS::get_singleton()->get_ticks_msec();
		PrepareImageMaps();
		after = OS::get_singleton()->get_ticks_msec();
		print_line("PrepareImageMaps took " + itos(after - before) + " ms");

		if (m_bCancel)
			return false;

		if (m_Logic_Process_AO) {
			print_line("ProcessAO");
			before = OS::get_singleton()->get_ticks_msec();
			ProcessAO();
			after = OS::get_singleton()->get_ticks_msec();
			print_line("ProcessAO took " + itos(after - before) + " ms");

			//Convolve_AO();
		}

		if (m_Logic_Process_Lightmap) {
			print_line("ProcessTexels");
			before = OS::get_singleton()->get_ticks_msec();
			if (m_Settings_Mode == LMMODE_BACKWARD)
				ProcessTexels();
			else {
				ProcessLights();
				ProcessEmissionTris();
			}
			DoAmbientBounces();
			after = OS::get_singleton()->get_ticks_msec();
			print_line("ProcessTexels took " + itos(after - before) + " ms");
		}

		if (m_Logic_Process_Probes) {
			// calculate probes
			print_line("ProcessProbes");
			if (LoadLightmap(out_image_lightmap)) {
				ProcessLightProbes();
			}
		}

		if (m_bCancel)
			return false;

		if (!m_Logic_Process_Probes) {
			WriteOutputImage_Lightmap(out_image_lightmap);
			WriteOutputImage_AO(out_image_ao);

			// test convolution
			//			Merge_AndWriteOutputImage_Combined(out_image_combined);
			//			out_image_combined.save_png("before_convolve.png");

			//			WriteOutputImage_AO(out_image_ao, true);
		}

	} // if not just merging
	else {
		// merging
		// need the meshes list for seam stitching
		m_Scene.FindMeshes(pMeshesRoot);

		// load the lightmap and ao from disk
		LoadLightmap(out_image_lightmap);
		LoadAO(out_image_ao);
	}

	//	print_line("WriteOutputImage");
	//	before = OS::get_singleton()->get_ticks_msec();
	if (m_Logic_Output_Final)
		Merge_AndWriteOutputImage_Combined(out_image_combined);
	//	after = OS::get_singleton()->get_ticks_msec();
	//	print_line("WriteOutputImage took " + itos(after -before) + " ms");

	// clear everything out of ram as no longer needed
	Reset();

	return true;
}

void LightMapper::ProcessTexels_AmbientBounce_Line_MT(uint32_t offset_y, int start_y) {
	int y = offset_y + start_y;

	for (int x = 0; x < m_iWidth; x++) {
		FColor power = ProcessTexel_AmbientBounce(x, y);

		// save the incoming light power in the mirror image (as the source is still being used)
		m_Image_L_mirror.GetItem(x, y) = power;
	}
}

void LightMapper::ProcessTexels_AmbientBounce(int section_size, int num_sections) {
	m_Image_L_mirror.Blank();

	// disable multithread
	//num_sections = 0;

	for (int s = 0; s < num_sections; s++) {
		int section_start = s * section_size;

		if (bake_step_function) {
			m_bCancel = bake_step_function(section_start, String("Process TexelsBounce: ") + " (" + itos(section_start) + ")");
			if (m_bCancel) {
				if (bake_end_function)
					bake_end_function();
				return;
			}
		}

		thread_process_array(section_size, this, &LightMapper::ProcessTexels_AmbientBounce_Line_MT, section_start);

		//		for (int n=0; n<section_size; n++)
		//		{
		//			ProcessTexel_Line_MT(n, section_start);
		//		}
	}

	int leftover_start = num_sections * section_size;

	for (int y = leftover_start; y < m_iHeight; y++) {
		if ((y % 10) == 0) {
			//			print_line("\tTexels bounce line " + itos(y));
			//			OS::get_singleton()->delay_usec(1);

			if (bake_step_function) {
				m_bCancel = bake_step_function(y, String("Process TexelsBounce: ") + " (" + itos(y) + ")");
				if (m_bCancel)
					return;
			}
		}

		ProcessTexels_AmbientBounce_Line_MT(y, 0);
		//		for (int x=0; x<m_iWidth; x++)
		//		{
		//			FColor power = ProcessTexel_Bounce(x, y);

		//			// save the incoming light power in the mirror image (as the source is still being used)
		//			m_Image_L_mirror.GetItem(x, y) = power;
		//		}
	}

	// merge the 2 luminosity maps
	for (int y = 0; y < m_iHeight; y++) {
		for (int x = 0; x < m_iWidth; x++) {
			FColor col = m_Image_L.GetItem(x, y);

			FColor col_add = m_Image_L_mirror.GetItem(x, y) * m_Settings_AmbientBouncePower;

			//			assert (col_add.r >= 0.0f);
			//			assert (col_add.g >= 0.0f);
			//			assert (col_add.b >= 0.0f);

			col += col_add;

			m_Image_L.GetItem(x, y) = col;
		}
	}
}

void LightMapper::Backward_TraceTriangles() {
	int nTris = m_Scene.m_Tris.size();

	if (bake_begin_function) {
		int progress_range = nTris;
		bake_begin_function(progress_range);
	}

	for (int n = 0; n < nTris; n++) {

		if ((n % 128) == 0) {
			if (bake_step_function) {
				m_bCancel = bake_step_function(n, String("Process Backward Tris: ") + " (" + itos(n) + ")");
				if (m_bCancel) {
					if (bake_end_function)
						bake_end_function();
					return;
				}
			}
		}

		Backward_TraceTriangle(n);
	}

	if (bake_end_function) {
		bake_end_function();
	}
}

void LightMapper::Backward_TraceTriangle(int tri_id) {

	const UVTri &tri = m_Scene.m_UVTris[tri_id];

	float area = tri.CalculateTwiceArea();
	if (area < 0.0f) area = -area;

	int nSamples = area * 400000.0f * m_AdjustedSettings.m_Backward_NumRays;

	for (int n = 0; n < nSamples; n++) {
		Backward_TraceSample(tri_id);
	}
}

void LightMapper::Backward_TraceSample(int tri_id) {
	const UVTri &uv_tri = m_Scene.m_UVTris[tri_id];
	const Tri &tri = m_Scene.m_Tris[tri_id];

	// choose random point on triangle
	Vector3 bary;
	RandomBarycentric(bary);

	// test, clamp the barycentric
	//	bary *= 0.998f;
	//	bary += Vector3(0.001f, 0.001f, 0.001f);

	// get position in world space
	Vector3 pos;
	tri.InterpolateBarycentric(pos, bary);

	// test, pull in a little
	//pos = pos.linear_interpolate(tri.GetCentre(), 0.001f);

	// get uv / texel position
	Vector2 uv;
	uv_tri.FindUVBarycentric(uv, bary);
	uv.x *= m_iWidth;
	uv.y *= m_iHeight;

	// round?
	int tx = uv.x;
	int ty = uv.y;
	//	int tx = Math::round(uv.x);
	//	int ty = Math::round(uv.y);

	// could be off the image
	FColor *pTexel = m_Image_L.Get(tx, ty);
	if (!pTexel)
		return;

	// apply bias
	// add epsilon to pos to prevent self intersection and neighbour intersection
	const Vector3 &plane_normal = m_Scene.m_TriPlanes[tri_id].normal;
	pos += plane_normal * m_Settings_SurfaceBias;

	// interpolated normal
	Vector3 normal;
	m_Scene.m_TriNormals[tri_id].InterpolateBarycentric(normal, bary.x, bary.y, bary.z);

	FColor temp;
	for (int l = 0; l < m_Lights.size(); l++) {
		BF_ProcessTexel_Light(Color(), l, pos, plane_normal, normal, temp, 1);
		*pTexel += temp;
	}

	// add emission
	Color emission_tex_color;
	Color emission_color;
	if (m_Scene.FindEmissionColor(tri_id, bary, emission_tex_color, emission_color)) {
		FColor femm;
		femm.Set(emission_tex_color);

		//		float power = m_Settings_Backward_RayPower * m_AdjustedSettings.m_Backward_NumRays * 128.0f;
		float power = m_Settings_Backward_RayPower * 128.0f;

		*pTexel += femm * power;
	}
}

void LightMapper::ProcessTexels() {
	//Backward_TraceTriangles();
	//return;

	//int nCores = OS::get_singleton()->get_processor_count();

	int section_size = m_iHeight / 64; //nCores;
	int num_sections = m_iHeight / section_size;
	int leftover_start = 0;

	if (bake_begin_function) {
		int progress_range = m_iHeight;
		bake_begin_function(progress_range);
	}

	m_iNumTests = 0;

	// prevent multithread
	//num_sections = 0;

	for (int s = 0; s < num_sections; s++) {
		int section_start = s * section_size;

		if (bake_step_function) {
			m_bCancel = bake_step_function(section_start, String("Process Texels: ") + " (" + itos(section_start) + ")");
			if (m_bCancel) {
				if (bake_end_function)
					bake_end_function();
				return;
			}
		}

		thread_process_array(section_size, this, &LightMapper::ProcessTexel_Line_MT, section_start);

		//		for (int n=0; n<section_size; n++)
		//		{
		//			ProcessTexel_Line_MT(n, section_start);
		//		}
	}

	leftover_start = num_sections * section_size;

	for (int y = leftover_start; y < m_iHeight; y++) {
		if ((y % 10) == 0) {
			//print_line("\tTexels line " + itos(y));
			//OS::get_singleton()->delay_usec(1);

			if (bake_step_function) {
				m_bCancel = bake_step_function(y, String("Process Texels: ") + " (" + itos(y) + ")");
				if (m_bCancel) {
					if (bake_end_function)
						bake_end_function();
					return;
				}
			}
		}

		ProcessTexel_Line_MT(y, 0);
	}

	//	m_iNumTests /= (m_iHeight * m_iWidth);
	print_line("num tests : " + itos(m_iNumTests));

	if (bake_end_function) {
		bake_end_function();
	}
}

void LightMapper::DoAmbientBounces() {
	if (!m_AdjustedSettings.m_NumAmbientBounces)
		return;

	int section_size = m_iHeight / 64; //nCores;
	int num_sections = m_iHeight / section_size;

	if (bake_begin_function) {
		int progress_range = m_iHeight;
		bake_begin_function(progress_range);
	}

	for (int b = 0; b < m_AdjustedSettings.m_NumAmbientBounces; b++) {
		ProcessTexels_AmbientBounce(section_size, num_sections);
	}

	if (bake_end_function) {
		bake_end_function();
	}
}

void LightMapper::ProcessTexel_Line_MT(uint32_t offset_y, int start_y) {
	int y = offset_y + start_y;

	for (int x = 0; x < m_iWidth; x++) {
		BF_ProcessTexel(x, y);
	}
}

bool LightMapper::Light_RandomSample(const LLight &light, const Vector3 &ptSurf, Ray &ray, Vector3 &ptLight, float &ray_length, float &multiplier) const {
	// for a spotlight, we can cull completely in a lot of cases.
	if (light.type == LLight::LT_SPOT) {
		Ray r;
		r.o = light.spot_emanation_point;
		r.d = ptSurf - r.o;
		r.d.normalize();
		float dot = r.d.dot(light.dir);
		//float angle = safe_acosf(dot);
		//if (angle >= light.spot_angle_radians)

		dot -= light.spot_dot_max;

		if (dot <= 0.0f)
			return false;
	}

	// default
	multiplier = 1.0f;

	switch (light.type) {
		case LLight::LT_DIRECTIONAL: {
			// we set this to maximum, better not to check at all but makes code simpler
			ray_length = FLT_MAX;

			ray.o = ptSurf;
			//r.d = -light.dir;

			// perturb direction depending on light scale
			//Vector3 ptTarget = light.dir * -2.0f;

			Vector3 offset;
			RandomUnitDir(offset);
			offset *= light.scale;

			offset += (light.dir * -2.0f);
			ray.d = offset.normalized();

			// disallow zero length (should be rare)
			if (ray.d.length_squared() < 0.00001f)
				return false;

			// don't allow from opposite direction
			if (ray.d.dot(light.dir) > 0.0f)
				ray.d = -ray.d;

			// this is used for intersection test to see
			// if e.g. sun is obscured. So we make the light position a long way from the origin
			ptLight = ptSurf + (ray.d * 10000000.0f);
		} break;
		case LLight::LT_SPOT: {
			// source
			float dot;

			while (true) {
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;

				ray.o = light.pos;
				ray.o += offset;

				// offset from origin to destination texel
				ray.d = ptSurf - ray.o;

				// normalize and find length
				ray_length = NormalizeAndFindLength(ray.d);

				dot = ray.d.dot(light.dir);
				//float angle = safe_acosf(dot);
				//if (angle >= light.spot_angle_radians)

				dot -= light.spot_dot_max;

				// if within cone, it is ok
				if (dot > 0.0f)
					break;
			}

			// store the light sample point
			ptLight = ray.o;

			// reverse ray for precision reasons
			ray.d = -ray.d;
			ray.o = ptSurf;

			dot *= 1.0f / (1.0f - light.spot_dot_max);
			multiplier = dot * dot;
			multiplier *= multiplier;
		} break;
		default: {
			Vector3 offset;
			RandomUnitDir(offset);
			offset *= light.scale;
			ptLight = light.pos + offset;

			// offset from origin to destination texel
			ray.o = ptSurf;
			ray.d = ptLight - ray.o;

			// normalize and find length
			ray_length = NormalizeAndFindLength(ray.d);

			// safety
			//assert (r.d.length() > 0.0f);

		} break;
	}

	// by this point...
	// ray should be set, ray_length, and ptLight
	return true;
}

// trace from the poly TO the light, not the other way round, to avoid precision errors
void LightMapper::BF_ProcessTexel_Light(const Color &orig_albedo, int light_id, const Vector3 &ptSource, const Vector3 &orig_face_normal, const Vector3 &orig_vertex_normal, FColor &color, int nSamples) //, uint32_t tri_ignore)
{
	const LLight &light = m_Lights[light_id];

	Ray r;

	float power = light.energy;
	power *= m_Settings_Backward_RayPower;

	//int nSamples = m_AdjustedSettings.m_Backward_NumRays;

	// total light hitting texel
	color.Set(0.0f);

	// for a spotlight, we can cull completely in a lot of cases.
	if (light.type == LLight::LT_SPOT) {
		r.o = light.spot_emanation_point;
		r.d = ptSource - r.o;
		r.d.normalize();
		float dot = r.d.dot(light.dir);
		//float angle = safe_acosf(dot);
		//if (angle >= light.spot_angle_radians)

		dot -= light.spot_dot_max;

		if (dot <= 0.0f)
			return;
	}

	// lower power of directionals - no distance falloff
	if (light.type == LLight::LT_DIRECTIONAL) {
		power *= 0.08f;
	}

	// each ray
	for (int n = 0; n < nSamples; n++) {
		Vector3 ptDest = light.pos;

		// allow falloff for cones
		float multiplier = 1.0f;

		// if the hit point is further that the distance to the light, then it doesn't count
		// must be set by the light type
		float ray_length;

		switch (light.type) {
			case LLight::LT_DIRECTIONAL: {
				// we set this to maximum, better not to check at all but makes code simpler
				ray_length = FLT_MAX;

				r.o = ptSource;
				//r.d = -light.dir;

				// perturb direction depending on light scale
				//Vector3 ptTarget = light.dir * -2.0f;

				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;

				offset += (light.dir * -2.0f);
				r.d = offset.normalized();

				// disallow zero length (should be rare)
				if (r.d.length_squared() < 0.00001f)
					continue;

				// don't allow from opposite direction
				if (r.d.dot(light.dir) > 0.0f)
					r.d = -r.d;
			} break;
			case LLight::LT_SPOT: {
				// source
				float dot;

				while (true) {
					Vector3 offset;
					RandomUnitDir(offset);
					offset *= light.scale;

					r.o = light.pos;
					r.o += offset;

					// offset from origin to destination texel
					r.d = ptSource - r.o;

					// normalize and find length
					ray_length = NormalizeAndFindLength(r.d);

					dot = r.d.dot(light.dir);
					//float angle = safe_acosf(dot);
					//if (angle >= light.spot_angle_radians)

					dot -= light.spot_dot_max;

					// if within cone, it is ok
					if (dot > 0.0f)
						break;
				}

				// reverse ray for precision reasons
				r.d = -r.d;
				r.o = ptSource;

				dot *= 1.0f / (1.0f - light.spot_dot_max);
				multiplier = dot * dot;
				multiplier *= multiplier;
			} break;
			default: {
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;
				ptDest += offset;

				// offset from origin to destination texel
				r.o = ptSource;
				r.d = ptDest - r.o;

				// normalize and find length
				ray_length = NormalizeAndFindLength(r.d);

				// safety
				//assert (r.d.length() > 0.0f);

			} break;
		}

		// only bother tracing if the light is in front of the surface normal
		float dot_light_surf = orig_vertex_normal.dot(r.d);
		if (dot_light_surf <= 0.0f)
			continue;

		//Vector3 ray_origin = r.o;
		FColor sample_color = light.color;
		int panic_count = 32;

		// for bounces
		Vector3 hit_face_normal = orig_face_normal;

		bool keep_tracing = true;
		while (keep_tracing) {
			keep_tracing = false;

			// collision detect
			float u, v, w, t;

			m_Scene.m_Tracer.m_bUseSDF = true;
			int tri = m_Scene.FindIntersect_Ray(r, u, v, w, t);
			//		m_Scene.m_Tracer.m_bUseSDF = false;
			//		int tri2 = m_Scene.IntersectRay(r, u, v, w, t, m_iNumTests);
			//		if (tri != tri2)
			//		{
			//			// repeat SDF version
			//			m_Scene.m_Tracer.m_bUseSDF = true;
			//			int tri = m_Scene.IntersectRay(r, u, v, w, t, m_iNumTests);
			//		}

			// further away than the destination?
			if (tri != -1) {
				if (t > ray_length)
					tri = -1;
				//				else
				//				{
				//					// we hit something, move the ray origin to the thing hit
				//					r.o += r.d * t;
				//				}
			}

			// nothing hit
			if (tri == -1)
			//if ((tri == -1) || (tri == (int) tri_ignore))
			{
				// for backward tracing, first pass, this is a special case, because we DO
				// take account of distance to the light, and normal, in order to simulate the effects
				// of the likelihood of 'catching' a ray. In forward tracing this happens by magic.
				float local_power;

				// no drop off for directional lights
				//float dist = (ptDest - ray_origin).length();
				float dist = ray_length;
				local_power = LightDistanceDropoff(dist, light, power);

				// take into account normal
				float dot = r.d.dot(orig_vertex_normal);
				dot = fabs(dot);

				local_power *= dot;

				// cone falloff
				local_power *= multiplier;

				// total color
				sample_color *= local_power;

				color += sample_color;

				// ONLY BOUNCE IF WE HIT THE TARGET SURFACE!!
				// start the bounces

				// probability of bounce based on alpha? NYI
				if (m_AdjustedSettings.m_NumDirectionalBounces) {
					// pass the direction as incoming (i.e. from the light to the surface, not the other way around)

					// the source ray is always the target surface (no matter whether we hit transparent or blockers)
					r.o = ptSource;
					r.d = -r.d;

					// we need to apply the color from the original target surface
					sample_color.r *= orig_albedo.r;
					sample_color.g *= orig_albedo.g;
					sample_color.b *= orig_albedo.b;

					// bounce ray direction
					if (BounceRay(r, orig_face_normal, false)) {
						// should this normal be plane normal or vertex normal?
						BF_ProcessTexel_LightBounce(m_AdjustedSettings.m_NumDirectionalBounces, r, sample_color);
					}
				}
			} else {
				// hit something, could be transparent

				// back face?
				bool bBackFace = HitBackFace(r, tri, Vector3(u, v, w), hit_face_normal);

				// first get the texture details
				Color albedo;
				bool bTransparent;
				m_Scene.FindPrimaryTextureColors(tri, Vector3(u, v, w), albedo, bTransparent);

				if (bTransparent) {
					bool opaque = bTransparent && (albedo.a > 0.999f);

					// push ray past the surface
					if (!opaque) {
						// hard code
						//						if (albedo.a > 0.0f)
						//							albedo.a = 0.3f;

						// position of potential hit
						Vector3 pos;
						const Tri &triangle = m_Scene.m_Tris[tri];
						triangle.InterpolateBarycentric(pos, u, v, w);

						float push = -m_Settings_SurfaceBias;
						if (bBackFace) push = -push;

						r.o = pos + (hit_face_normal * push);

						// apply the color to the ray
						CalculateTransmittance(albedo, sample_color);

						// this shouldn't happen, but if it did we'd get an infinite loop if we didn't break out
						panic_count--;
						if (!panic_count)
							break; // break out of while .. does this work?

						keep_tracing = true;
					}
				}
			}

		} // while keep tracing

	} // for n through samples

	// the color is returned in color
}

bool LightMapper::BounceRay(Ray &r, const Vector3 &face_normal, bool apply_epsilon) {
	float face_dot = face_normal.dot(r.d);

	// back face, don't bounce
	if (face_dot >= 0.0f)
		return false;

	// BOUNCING - mirror
	Vector3 mirror_dir = r.d - (2.0f * (face_dot * face_normal));

	// random hemisphere
	Vector3 hemi_dir;

	const float range = 1.0f;
	while (true) {
		hemi_dir.x = Math::random(-range, range);
		hemi_dir.y = Math::random(-range, range);
		hemi_dir.z = Math::random(-range, range);

		float sl = hemi_dir.length_squared();
		if (sl > 0.0001f) {
			break;
		}
	}
	// compare direction to normal, if opposite, flip it
	if (hemi_dir.dot(face_normal) < 0.0f)
		hemi_dir = -hemi_dir;

	r.d = hemi_dir.linear_interpolate(mirror_dir, m_Settings_Smoothness);

	// standard epsilon? NYI
	if (apply_epsilon)
		r.o += (face_normal * m_Settings_SurfaceBias); //0.01f);

	return true;
}

void LightMapper::BF_ProcessTexel_LightBounce(int bounces_left, Ray r, FColor ray_color) {
	// apply bounce power
	ray_color *= m_Settings_DirectionalBouncePower;

	float u, v, w, t;
	int tri = m_Scene.FindIntersect_Ray(r, u, v, w, t);

	// nothing hit
	if (tri == -1) {
		// terminate bounces (or detect sky)
		return;
	}

	//	if (bounces_left == 1)
	//	{
	//		print_line("test");
	//	}

	// hit the back of a face? if so terminate ray
	Vector3 vertex_normal;
	const Tri &triangle_normal = m_Scene.m_TriNormals[tri];
	triangle_normal.InterpolateBarycentric(vertex_normal, u, v, w);
	vertex_normal.normalize(); // is this necessary as we are just checking a dot product polarity?

	// first get the texture details
	Color albedo;
	bool bTransparent;
	m_Scene.FindPrimaryTextureColors(tri, Vector3(u, v, w), albedo, bTransparent);
	bool pass_through = bTransparent && (albedo.a < 0.001f);

	bool bBackFace = false;
	const Vector3 &face_normal = m_Scene.m_TriPlanes[tri].normal;

	//float face_dot = face_normal.dot(r.d);
	//	if (face_dot >= 0.0f)
	//		bBackFace = true;

	float vertex_dot = vertex_normal.dot(r.d);
	if (vertex_dot >= 0.0f)
		bBackFace = true;

	// if not transparent and backface, then terminate ray
	if (bBackFace && !bTransparent)
		return;

	// convert barycentric to uv coords in the lightmap
	Vector2 uv;
	m_Scene.FindUVsBarycentric(tri, uv, u, v, w);

	// texel address
	int tx = uv.x * m_iWidth;
	int ty = uv.y * m_iHeight;

	// could be off the image
	if (!m_Image_L.IsWithin(tx, ty))
		return;

	// position of potential hit
	Vector3 pos;
	const Tri &triangle = m_Scene.m_Tris[tri];
	triangle.InterpolateBarycentric(pos, u, v, w);

	// deal with tranparency
	if (bTransparent) {
		// if not passing through, because clear, chance of pass through
		if (!pass_through && !bBackFace) {
			pass_through = Math::randf() > albedo.a;
		}

		// if the ray is passing through, we want to calculate the color modified by the surface
		if (pass_through)
			CalculateTransmittance(albedo, ray_color);

		// if pass through
		if (bBackFace || pass_through) {
			// push the ray origin through the hit surface
			float push = -0.001f; // 0.001
			if (bBackFace) push = -push;

			r.o = pos + (face_normal * push);

			// call recursively
			BF_ProcessTexel_LightBounce(bounces_left, r, ray_color);
			return;
		}

	} // if transparent

	///////

	// if we got here, it is front face and either solid or no pass through,
	// so there is a hit
	float lambert = MAX(0.0f, -vertex_dot);
	// apply lambert diffuse
	ray_color *= lambert;

	// apply
	MT_SafeAddToTexel(tx, ty, ray_color);

	// any more bounces to go?
	bounces_left--;
	if (bounces_left <= 0)
		return;

	// move the ray origin up to the hit point
	// (we will use the barycentric derived pos here, rather than r.o + (r.d * t)
	// for more accuracy relative to the surface
	r.o = pos;

	// bounce and lower power
	FColor falbedo;
	falbedo.Set(albedo);

	// bounce color
	ray_color = ray_color * falbedo;

	// find the bounced direction from the face normal
	BounceRay(r, face_normal);

	// call recursively
	BF_ProcessTexel_LightBounce(bounces_left, r, ray_color);
}

void LightMapper::ProcessLightProbes() {
	LightProbes probes;
	int stages = probes.Create(*this);
	if (stages != -1) {
		if (bake_begin_function) {
			bake_begin_function(stages);
		}

		for (int n = 0; n < stages; n++) {
			if (bake_step_function) {
				m_bCancel = bake_step_function(n, String("Process LightProbes: ") + " (" + itos(n) + ")");
				if (m_bCancel)
					break;
			}

			probes.Process(n);
		}

		if (!m_bCancel)
			probes.Save();

		if (bake_end_function) {
			bake_end_function();
		}
	}
}

FColor LightMapper::ProcessTexel_AmbientBounce(int x, int y) {
	FColor total;
	total.Set(0.0f);

	// find triangle
	uint32_t tri_source = *m_Image_ID_p1.Get(x, y);
	if (!tri_source)
		return total;
	tri_source--; // plus one based

	// barycentric
	const Vector3 &bary = *m_Image_Barycentric.Get(x, y);

	Vector3 pos;
	m_Scene.m_Tris[tri_source].InterpolateBarycentric(pos, bary.x, bary.y, bary.z);

	Vector3 norm;
	const Tri &triangle_normal = m_Scene.m_TriNormals[tri_source];
	triangle_normal.InterpolateBarycentric(norm, bary.x, bary.y, bary.z);
	norm.normalize();

	// precalculate normal push and ray origin on top of the surface
	const Vector3 &plane_norm = m_Scene.m_TriPlanes[tri_source].normal;
	Vector3 normal_push = plane_norm * m_Settings_SurfaceBias;
	Vector3 ray_origin = pos + normal_push;

	int nSamples = m_AdjustedSettings.m_NumAmbientBounceRays;
	int samples_counted = nSamples;

	for (int n = 0; n < nSamples; n++) {
		if (!ProcessTexel_AmbientBounce_Sample(plane_norm, ray_origin, total)) {
			samples_counted--;
		}
	}

	// some samples may have been missed due to transparency
	if (samples_counted)
		return total / samples_counted;

	return total;
}

FColor LightMapper::Probe_CalculateIndirectLight(const Vector3 &pos) {
	FColor total;
	total.Set(0.0f);

	// we will reuse the same routine from texel bounce
	int nSamples = m_Settings_ProbeSamples;
	int samples_counted = nSamples;

	Vector3 ray_origin = pos;

	for (int n = 0; n < nSamples; n++) {
		// random normal
		Vector3 norm;
		RandomUnitDir(norm);

		if (!ProcessTexel_AmbientBounce_Sample(norm, ray_origin, total)) {
			samples_counted--;
		}
	}

	// some samples may have been missed due to transparency
	if (samples_counted)
		return total / samples_counted;

	return total;
}

bool LightMapper::ProcessTexel_AmbientBounce_Sample(const Vector3 &plane_norm, const Vector3 &ray_origin, FColor &total_col) {
	// first dot
	Ray r;
	r.o = ray_origin;

	// SLIDING
	//			Vector3 temp = r.d.cross(norm);
	//			new_ray.d = norm.cross(temp);

	// BOUNCING - mirror
	//new_ray.d = r.d - (2.0f * (dot * norm));

	// random hemisphere
	RandomUnitDir(r.d);

	// compare direction to normal, if opposite, flip it
	if (r.d.dot(plane_norm) < 0.0f)
		r.d = -r.d;

	// loop here just in case transparent
	while (true) {

		// collision detect
		Vector3 bary;
		float t;
		int tri_hit = m_Scene.FindIntersect_Ray(r, bary, t);

		// nothing hit
		//	if ((tri_hit != -1) && (tri_hit != (int) tri_source))
		if (tri_hit != -1) {
			// look up the UV of the tri hit
			Vector2 uvs;
			m_Scene.FindUVsBarycentric(tri_hit, uvs, bary);

			// find texel
			int dx = (uvs.x * m_iWidth); // round?
			int dy = (uvs.y * m_iHeight);

			// texel not on the UV map
			if (!m_Image_L.IsWithin(dx, dy))
				return true;

			// back face?
			Vector3 face_normal;
			bool bBackFace = HitBackFace(r, tri_hit, bary, face_normal);

			// the contribution is the luminosity at that spot and the albedo
			Color albedo;
			bool bTransparent;
			m_Scene.FindPrimaryTextureColors(tri_hit, bary, albedo, bTransparent);

			FColor falbedo;
			falbedo.Set(albedo);

			bool opaque = !(bTransparent && (albedo.a < 0.5f));

			// if the surface transparent we may want to downgrade the influence a little
			if (bTransparent)
				falbedo *= albedo.a;

			// see through has no effect on colour
			if (!opaque || (bBackFace && bTransparent)) {
				// if it is front facing, still apply the 'haze' from the colour
				// (this all counts as 1 sample, so could get super light in theory, if passing through several hazes...)
				if (!bBackFace)
					total_col += (m_Image_L.GetItem(dx, dy) * falbedo);

				// move the ray origin and do again
				// position of potential hit
				Vector3 pos;
				const Tri &triangle = m_Scene.m_Tris[tri_hit];
				triangle.InterpolateBarycentric(pos, bary);

				float push = -m_Settings_SurfaceBias;
				if (bBackFace) push = -push;

				r.o = pos + (face_normal * push);

				// and repeat the loop
			} else {
				total_col += (m_Image_L.GetItem(dx, dy) * falbedo);

				//assert (total_col.r >= 0.0f);
				break;
			}
		} else {
			// nothing hit, exit loop
			break;
		}

	} // while true (loop while transparent)

	return true;
}

//	if ((tx == 3) && (ty == 17))
//	{
//		print_line("test");
//	}

// find triangle
//	uint32_t tri = *m_Image_ID_p1.Get(tx, ty);
//	if (!tri)
//		return;
//	tri--; // plus one based

// may be more than 1 triangle on this texel
//	uint32_t tri2 = *m_Image_ID2_p1.Get(tx, ty);

// barycentric
//	const Vector3 &bary = *m_Image_Barycentric.Get(tx, ty);

//	Vector3 pos;
//	m_Scene.m_Tris[tri].InterpolateBarycentric(pos, bary.x, bary.y, bary.z);

//	Vector3 normal;
//	m_Scene.m_TriNormals[tri].InterpolateBarycentric(normal, bary.x, bary.y, bary.z);

void LightMapper::BF_ProcessTexel(int tx, int ty) {
	//		if ((tx == 13) && (ty == 284))
	//			print_line("testing");

	// find triangle
	uint32_t tri_id = *m_Image_ID_p1.Get(tx, ty);
	if (!tri_id)
		return;
	tri_id--; // plus one based

	// barycentric
	const Vector3 &bary = *m_Image_Barycentric.Get(tx, ty);
	Vector3 bary_clamped; // = bary;

	// new .. cap the barycentric to prevent edge artifacts
	const float clamp_margin = 0.0001f;
	const float clamp_margin_high = 1.0f - clamp_margin;

	bary_clamped.x = CLAMP(bary.x, clamp_margin, clamp_margin_high);
	bary_clamped.y = CLAMP(bary.y, clamp_margin, clamp_margin_high);
	bary_clamped.z = CLAMP(bary.z, clamp_margin, clamp_margin_high);

	// we will trace
	// FROM THE SURFACE TO THE LIGHT!!
	// this is very important, because the ray is origin and direction,
	// and there will be precision errors at the destination.
	// At the light end this doesn't matter, but if we trace the other way
	// we get artifacts due to precision loss due to normalized direction.
	Vector3 pos;
	m_Scene.m_Tris[tri_id].InterpolateBarycentric(pos, bary_clamped);

	// add epsilon to pos to prevent self intersection and neighbour intersection
	const Vector3 &plane_normal = m_Scene.m_TriPlanes[tri_id].normal;
	pos += plane_normal * m_Settings_SurfaceBias;

	Vector3 normal;
	m_Scene.m_TriNormals[tri_id].InterpolateBarycentric(normal, bary.x, bary.y, bary.z);

	//Vector2i tex_uv = Vector2i(x, y);

	// could be off the image
	if (!m_Image_L.IsWithin(tx, ty))
		return;

	// find the colors of this texel
	Color albedo;
	Color emission;
	bool transparent;
	bool emitter;
	m_Scene.FindAllTextureColors(tri_id, bary, albedo, emission, transparent, emitter);

	FColor texel_add;
	texel_add.Set(0.0f);

	//	FColor * pTexel = m_Image_L.Get(tx, ty);
	//	if (!pTexel)
	//		return;

	int nSamples = m_AdjustedSettings.m_Backward_NumRays;
	FColor temp;
	for (int l = 0; l < m_Lights.size(); l++) {
		BF_ProcessTexel_Light(albedo, l, pos, plane_normal, normal, temp, nSamples);
		texel_add += temp;
	}

	// add emission
	//	Color emission_tex_color;
	//	Color emission_color;
	//	if (m_Scene.FindEmissionColor(tri_id, bary, emission_tex_color, emission_color))
	if (emitter) {
		// Glow determines how much the surface itself is lighted (and thus the ratio between glow and emission)
		// emission density determines the number of rays and lighting effect
		FColor femm;
		femm.Set(emission);
		float power = m_Settings_Backward_RayPower * m_AdjustedSettings.m_Backward_NumRays * 32.0f;
		power *= m_Settings_Glow;
		texel_add += femm * power;

		// only if directional bounces are being used (could use ambient bounces for emission)
		if (m_AdjustedSettings.m_NumDirectionalBounces) {
			// needs to be adjusted according to size of texture .. as there will be more emissive texels
			int nSamples = m_AdjustedSettings.m_Backward_NumRays * 2 * m_AdjustedSettings.m_EmissionDensity;

			// apply the albedo to the emission color to get the color emanating
			femm.r *= albedo.r;
			femm.g *= albedo.g;
			femm.b *= albedo.b;

			Ray r;
			r.o = pos;
			femm *= m_Settings_Backward_RayPower * 128.0f * (1.0f / m_AdjustedSettings.m_EmissionDensity);

			for (int n = 0; n < nSamples; n++) {
				// send out emission bounce rays
				RandomUnitDir(r.d);

				float dot = plane_normal.dot(r.d);
				if (dot < 0.0f)
					r.d = -r.d;

				BF_ProcessTexel_LightBounce(m_AdjustedSettings.m_NumDirectionalBounces, r, femm); // always at least 1 emission ray
			}
		}
	}

	// safe write
	MT_SafeAddToTexel(tx, ty, texel_add);
}

void LightMapper::ProcessEmissionTris() {
	int num_sections = m_iNumRays / m_iRaysPerSection;

	if (!num_sections)
		num_sections = 1;

	float fraction = 1.0f / num_sections;

	if (bake_begin_function) {
		bake_begin_function(num_sections);
	}

	for (int s = 0; s < num_sections; s++) {
		if ((s % 1) == 0) {
			if (bake_step_function) {
				m_bCancel = bake_step_function(s, String("Process Emission Section: ") + " (" + itos(s) + ")");
				if (m_bCancel) {
					if (bake_end_function)
						bake_end_function();
					return;
				}
			}
		}

		ProcessEmissionTris_Section(fraction);

		while (!RayBank_AreVoxelsClear())
		//for (int b=0; b<m_AdjustedSettings.m_Forward_NumBounces+1; b++)
		{
			RayBank_Process();
			RayBank_Flush();
		} // for bounce
	} // for section

	// left over
	//	{
	//		int num_leftover = m_iNumRays - (num_sections * m_iRaysPerSection);
	//		ProcessLight(n, num_leftover);

	//		for (int b=0; b<m_AdjustedSettings.m_Forward_NumBounces+1; b++)
	//		{
	//			RayBank_Process();
	//			RayBank_Flush();
	//		} // for bounce
	//	}

	if (bake_end_function) {
		bake_end_function();
	}
}

void LightMapper::ProcessEmissionTris_Section(float fraction_of_total) {
	for (int n = 0; n < m_Scene.m_EmissionTris.size(); n++) {
		ProcessEmissionTri(n, fraction_of_total);
	}
}

void LightMapper::ProcessEmissionTri(int etri_id, float fraction_of_total) {
	const EmissionTri &etri = m_Scene.m_EmissionTris[etri_id];
	int tri_id = etri.tri_id;

	// get the material

	// positions
	const Tri &tri_pos = m_Scene.m_Tris[tri_id];

	// normals .. just use plane normal for now (no interpolation)
	Vector3 norm = m_Scene.m_TriPlanes[tri_id].normal;

	Ray ray;
	ray.d = norm;

	// use the area to get number of samples
	float rays_per_unit_area = m_iNumRays * m_AdjustedSettings.m_EmissionDensity * 0.12f * 0.5f;
	int nSamples = etri.area * rays_per_unit_area * fraction_of_total;

	// nSamples may be zero incorrectly for small triangles, maybe we need to adjust for this
	// NYI

	for (int s = 0; s < nSamples; s++) {
		// find a random barycentric coord
		Vector3 bary;
		RandomBarycentric(bary);

		// find point on this actual triangle
		tri_pos.InterpolateBarycentric(ray.o, bary);

		// shoot a ray from this pos and normal, using the emission color
		RandomUnitDir(ray.d);
		ray.d += norm;

		if (ray.d.dot(norm) < 0.0f)
			ray.d = -ray.d;

		ray.d.normalize();

		// get the albedo etc
		Color emission_tex_color;
		Color emission_color;
		m_Scene.FindEmissionColor(tri_id, bary, emission_tex_color, emission_color);

		FColor fcol;
		fcol.Set(emission_tex_color);
		RayBank_RequestNewRay(ray, m_AdjustedSettings.m_NumDirectionalBounces + 1, fcol);

		// special. For emission we want to also affect the emitting surface.
		// convert barycentric to uv coords in the lightmap
		Vector2 uv;
		m_Scene.FindUVsBarycentric(tri_id, uv, bary);

		// texel address
		int tx = uv.x * m_iWidth;
		int ty = uv.y * m_iHeight;

		// could be off the image
		if (!m_Image_L.IsWithin(tx, ty))
			continue;

		FColor *pTexelCol = m_Image_L.Get(tx, ty);

		//		fcol.Set(emission_color * 0.5f);
		*pTexelCol += fcol * m_Settings_Glow;
	}
}

void LightMapper::ProcessLights() {
	//	const int rays_per_section = 1024 * 16;

	int num_sections = m_iNumRays / m_iRaysPerSection;

	if (bake_begin_function) {
		bake_begin_function(num_sections);
	}

	for (int n = 0; n < m_Lights.size(); n++) {
		for (int s = 0; s < num_sections; s++) {
			// double check the voxels are clear
#ifdef DEBUG_ENABLED
			//RayBank_CheckVoxelsClear();
#endif

			if ((s % 1) == 0) {
				if (bake_step_function) {
					m_bCancel = bake_step_function(s, String("Process Light Section: ") + " (" + itos(s) + ")");
					if (m_bCancel) {
						if (bake_end_function)
							bake_end_function();
						return;
					}
				}
			}

			ProcessLight(n, m_iRaysPerSection);

			while (!RayBank_AreVoxelsClear())
			//for (int b=0; b<m_AdjustedSettings.m_Forward_NumBounces+1; b++)
			{
				RayBank_Process();
				RayBank_Flush();
			} // for bounce

			// if everything went correctly,
			// all the raybank should be empty
			//RayBank_AreVoxelsClear();

		} // for section

		// left over
		{
			int num_leftover = m_iNumRays - (num_sections * m_iRaysPerSection);
			ProcessLight(n, num_leftover);

			while (!RayBank_AreVoxelsClear())
			//for (int b=0; b<m_AdjustedSettings.m_Forward_NumBounces+1; b++)
			{
				RayBank_Process();
				RayBank_Flush();
			} // for bounce
		}

		// this is not really required, but is a protection against memory leaks
		RayBank_Reset(true);
	} // for light

	if (bake_end_function) {
		bake_end_function();
	}
}

void LightMapper::ProcessLight(int light_id, int num_rays) {
	const LLight &light = m_Lights[light_id];

	Ray r;
	//	r.o = Vector3(0, 5, 0);

	//	float range = light.scale.x;
	//	const float range = 2.0f;

	// the power should depend on the volume, with 1x1x1 being normal power
	//	float power = light.scale.x * light.scale.y * light.scale.z;
	//	float power = light.energy;
	//	power *= m_Settings_Forward_RayPower;

	//	light_color.r = power;
	//	light_color.g = power;
	//	light_color.b = power;

	// each ray

	//	num_rays = 1; // debug

	// new... allow the use of indirect energy to scale the number of samples
	num_rays *= light.indirect_energy;

	// compensate for the number of rays in terms of the power per ray
	float power = 0.01f; //m_Settings_Forward_RayPower;

	if (light.indirect_energy > 0.0001f)
		power *= 1.0f / light.indirect_energy;

	// for directional, we need a load more rays for it to work well - it is expensive
	if (light.type == LLight::LT_DIRECTIONAL) {
		num_rays *= 2;
		//		float area = light.dl_tangent_range * light.dl_bitangent_range;
		//		num_rays = num_rays * area;
		// we will increase the power as well, because daylight more powerful than light bulbs typically.
		power *= 4.0f;

		if (light.dir.y >= 0.0f)
			return; // not supported
	}

	FColor light_color = light.color * power;

	for (int n = 0; n < num_rays; n++) {
		//		if ((n % 10000) == 0)
		//		{
		//			if (bake_step_function)
		//			{
		//				m_bCancel = bake_step_function(n, String("Process Rays: ") + " (" + itos(n) + ")");
		//				if (m_bCancel)
		//					return;
		//			}
		//		}

		r.o = light.pos;

		switch (light.type) {
			case LLight::LT_DIRECTIONAL: {
				bool hit_bound = false;
				while (!hit_bound) {
					hit_bound = true;

					// use the precalculated source plane stored in the llight
					r.o = light.dl_plane_pt;
					r.o += light.dl_tangent * Math::random(0.0f, light.dl_tangent_range);
					r.o += light.dl_bitangent * Math::random(0.0f, light.dl_bitangent_range);

					bool direction_ok = false;
					while (!direction_ok) {
						RandomUnitDir(r.d);
						r.d *= light.scale;
						// must point down - reverse hemisphere if pointing up
						if (light.dir.dot(r.d) < 0.0f) {
							r.d = -r.d;
						}

						r.d += (light.dir * 2.0f);
						r.d.normalize();

						if (r.d.y < 0.0f) {
							direction_ok = true;
						}
					}

					// now we are starting from destination, we must trace back to origin
					const AABB &bb = GetTracer().GetWorldBound_expanded();

					// num units to top of map
					// plus a bit to make intersecting the world bound easier
					float units = (bb.size.y + 1.0f) / r.d.y;

					r.o += r.d * units;

					//				r.o.y = 3.2f;

					// hard code
					//				r.o = Vector3(0.2, 10, 0.2);
					//				r.d = Vector3(0, -1, 0);

					// special .. for dir light .. will it hit the AABB? if not, do a wraparound
					//					Vector3 clip;
					//					if (!GetTracer().IntersectRayAABB(r, GetTracer().m_SceneWorldBound_mid, clip))
					//					{
					//						hit_bound = false;
					//					}

					//					Vector3 ptHit = r.o + (r.d * 2.0f);
					//					Vector3 bb_min = bb.position;
					//					Vector3 bb_max = bb.position + bb.size;

					//					if (ptHit.x > bb_max.x)
					//						r.o.x -= bb.size.x;
					//					if (ptHit.x < bb_min.x)
					//						r.o.x += bb.size.x;
					//					if (ptHit.z > bb_max.z)
					//						r.o.z -= bb.size.z;
					//					if (ptHit.z < bb_min.z)
					//						r.o.z += bb.size.z;

					if (!GetTracer().GetWorldBound_expanded().intersects_ray(r.o, r.d)) {
						hit_bound = false;
					}
				} // while hit bound
			} break;
			case LLight::LT_SPOT: {
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;
				r.o += offset;

				r.d = light.dir;

				// random axis
				Vector3 axis;
				RandomAxis(axis);

				float falloff_start = 0.5f; // this could be adjustable;

				float ang_max = light.spot_angle_radians;
				float ang_falloff = ang_max * falloff_start;
				float ang_falloff_range = ang_max - ang_falloff;

				// random angle giving equal circle distribution
				float a = Math::random(0.0f, 1.0f);

				// below a certain proportion are the central, outside is falloff
				a *= 2.0f;
				if (a > 1.0f) {
					// falloff
					a -= 1.0f;
					//angle /= 3.0f;

					//a *= a;

					a = ang_falloff + (a * ang_falloff_range);
				} else {
					// central
					a = sqrtf(a);
					a *= ang_falloff;
				}

				//				if (angle > ang_falloff)
				//				{
				//					// in the falloff zone, use different math.
				//					float r = Math::random(0.0f, 1.0f);
				//					r = r*r;
				//					r *= ang_falloff_range;
				//					angle = ang_falloff + r;
				//				}

				Quat rot;
				rot.set_axis_angle(axis, a);

				// TURNED OFF FOR TESTING
				r.d = rot.xform(r.d);

				// this is the radius of the cone at distance 1
				//				float radius_at_dist_one = Math::tan(Math::deg2rad(light.spot_angle));

				//				float spot_ball_size = radius_at_dist_one;

				//				Vector3 ball;
				//				RandomSphereDir(ball, spot_ball_size);

				//				r.d += ball;
				//				r.d.normalize();
			} break;
			default: {
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;
				r.o += offset;

				RandomUnitDir(r.d);
			} break;
		}
		//r.d.normalize();

		RayBank_RequestNewRay(r, m_AdjustedSettings.m_NumDirectionalBounces + 1, light_color, 0);

		//		ProcessRay(r, 0, power);
	}

	//print_line("PASSED " + itos (passed) + ", MISSED " + itos(missed));
}

} // namespace LM
