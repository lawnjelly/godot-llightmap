#include "llightmapper.h"
#include "core/os/threaded_array_processor.h"
#include "lmerger.h"
#include "lunmerger.h"
#include "lscene_saver.h"
#include "scene/resources/packed_scene.h"
#include "llightscene.h"


namespace LM {


bool LightMapper::uv_map_meshes(Spatial * pRoot)
{
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
	MeshInstance * pMerged = m.Merge(pRoot, m_Settings_UVPadding);
	if (!pMerged)
	{
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
	if (m_Settings_UVFilename != "")
	{
		Node * pOrigOwner = pRoot->get_owner();

		saver.SaveScene(pRoot, m_Settings_UVFilename, true);

		if (replace_mesh_scene)
		{
			// delete the orig
			Node * pParent = pRoot->get_parent();

			// rename the old scene to prevent naming conflict
			pRoot->set_name("ToBeDeleted");
			pRoot->queue_delete();

			// now load the new file
			//ResourceLoader::import(m_Settings_UVFilename);

			Ref<PackedScene> ps = ResourceLoader::load(m_Settings_UVFilename, "PackedScene");
			if (ps.is_null())
			{
				if (bake_end_function) {
					bake_end_function();
				}
				return res;
			}

			Node * pFinalScene = ps->instance();
			if (pFinalScene)
			{
				pParent->add_child(pFinalScene);

				// set owners
				saver.SetOwnerRecursive(pFinalScene, pFinalScene);
				pFinalScene->set_owner(pOrigOwner);
			}

		} // if replace the scene
		else
		{
			// delete, to save confusion
			pRoot->queue_delete();
		}
	}

	if (bake_end_function) {
		bake_end_function();
	}

	return res;
}


bool LightMapper::lightmap_mesh(Spatial * pMeshesRoot, Spatial * pLR, Image * pIm_Lightmap, Image * pIm_AO, Image * pIm_Combined)
{
	CalculateQualityAdjustedSettings();

	// get the output dimensions before starting, because we need this
	// to determine number of rays, and the progress range
	m_iWidth = pIm_Combined->get_width();
	m_iHeight = pIm_Combined->get_height();
	m_iNumRays = m_AdjustedSettings.m_Forward_NumRays;

		int nTexels = m_iWidth * m_iHeight;

	// set num rays depending on method
		if (m_Settings_Mode == LMMODE_FORWARD)
		{
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

void LightMapper::Reset()
{
	m_Lights.clear(true);
	m_Scene.Reset();
	RayBank_Reset();
}

void LightMapper::Refresh_Process_State()
{
	// process states
	switch (m_Settings_BakeMode)
	{
	case LMBAKEMODE_LIGHTMAP:
		{
			m_Settings_Process_Lightmap = true;
			m_Settings_Process_AO = false;
		}
		break;
	case LMBAKEMODE_AO:
		{
			m_Settings_Process_Lightmap = false;
			m_Settings_Process_AO = true;
		}
		break;
	case LMBAKEMODE_MERGE:
		{
			m_Settings_Process_Lightmap = false;
			m_Settings_Process_AO = false;
		}
		break;
	default:
		{
			m_Settings_Process_Lightmap = true;
			m_Settings_Process_AO = true;
		}
		break;
	}
}


bool LightMapper::LightmapMesh(Spatial * pMeshesRoot, const Spatial &light_root, Image &out_image_lightmap, Image &out_image_ao, Image &out_image_combined)
{
	// print out settings
	print_line("Lightmap mesh");
	print_line("\tnum_bounces " + itos(m_AdjustedSettings.m_Forward_NumBounces));
	print_line("\tbounce_power " + String(Variant(m_Settings_Forward_BouncePower)));

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
	m_Image_AO.Create(m_iWidth, m_iHeight);

	if (m_Settings_BakeMode != LMBAKEMODE_MERGE)
	{

		m_QMC.Create(m_AdjustedSettings.m_AO_Samples);

		uint32_t before, after;
		FindLights_Recursive(&light_root);
		print_line("Found " + itos (m_Lights.size()) + " lights.");


		m_Image_ID_p1.Create(m_iWidth, m_iHeight);
		m_Image_ID2_p1.Create(m_iWidth, m_iHeight);

		m_Image_TriIDs.Create(m_iWidth, m_iHeight);
		m_TriIDs.clear(true);

		m_Image_Barycentric.Create(m_iWidth, m_iHeight);

		m_Image_Cuts.Create(m_iWidth, m_iHeight);
		//m_CuttingTris.clear(true);

		print_line("Scene Create");
		before = OS::get_singleton()->get_ticks_msec();
		if (!m_Scene.Create(pMeshesRoot, m_iWidth, m_iHeight, m_Settings_VoxelDensity, m_AdjustedSettings.m_Max_Material_Size, m_AdjustedSettings.m_Forward_Emission_Density))
			return false;

		PrepareLights();

		RayBank_Create();

		after = OS::get_singleton()->get_ticks_msec();
		print_line("SceneCreate took " + itos(after -before) + " ms");

		if (m_bCancel)
			return false;

		print_line("PrepareImageMaps");
		before = OS::get_singleton()->get_ticks_msec();
		PrepareImageMaps();
		after = OS::get_singleton()->get_ticks_msec();
		print_line("PrepareImageMaps took " + itos(after -before) + " ms");

		if (m_bCancel)
			return false;

		if (m_Settings_Process_AO)
		{
			print_line("ProcessAO");
			before = OS::get_singleton()->get_ticks_msec();
			ProcessAO();
			after = OS::get_singleton()->get_ticks_msec();
			print_line("ProcessAO took " + itos(after -before) + " ms");
		}


		if (m_Settings_Process_Lightmap)
		{
			print_line("ProcessTexels");
			before = OS::get_singleton()->get_ticks_msec();
			if (m_Settings_Mode == LMMODE_BACKWARD)
				ProcessTexels();
			else
			{
				ProcessLights();
				ProcessEmissionTris();
			}
			after = OS::get_singleton()->get_ticks_msec();
			print_line("ProcessTexels took " + itos(after -before) + " ms");
		}

		if (m_bCancel)
			return false;

		WriteOutputImage_Lightmap(out_image_lightmap);
		WriteOutputImage_AO(out_image_ao);

	} // if not just merging
	else
	{
		// merging, load the lightmap and ao from disk
		LoadLightmap(out_image_lightmap);
		LoadAO(out_image_ao);
	}

	//	print_line("WriteOutputImage");
	//	before = OS::get_singleton()->get_ticks_msec();
	Merge_AndWriteOutputImage_Combined(out_image_combined);
	//	after = OS::get_singleton()->get_ticks_msec();
	//	print_line("WriteOutputImage took " + itos(after -before) + " ms");

	// clear everything out of ram as no longer needed
	Reset();

	return true;
}

void LightMapper::ProcessTexels_Bounce()
{
	m_Image_L_mirror.Blank();


	for (int y=0; y<m_iHeight; y++)
	{
		if ((y % 10) == 0)
		{
			//			print_line("\tTexels bounce line " + itos(y));
			//			OS::get_singleton()->delay_usec(1);

			if (bake_step_function) {
				m_bCancel = bake_step_function(y, String("Process TexelsBounce: ") + " (" + itos(y) + ")");
				if (m_bCancel)
					return;
			}
		}

		for (int x=0; x<m_iWidth; x++)
		{
			FColor power = ProcessTexel_Bounce(x, y);

			// save the incoming light power in the mirror image (as the source is still being used)
			m_Image_L_mirror.GetItem(x, y) = power;
		}
	}

	// merge the 2 luminosity maps
	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			FColor col = m_Image_L.GetItem(x, y);

			col += (m_Image_L_mirror.GetItem(x, y) * m_Settings_Backward_BouncePower);

			m_Image_L.GetItem(x, y) = col;
		}
	}

}


void LightMapper::ProcessTexels()
{
	// set num rays depending on method
	//	if (m_Settings_Mode == LMMODE_FORWARD)
	//	{
	//		// the num rays / texel. This is per light!
	//		m_iNumRays *= nTexels;
	//		progress_range = m_iNumRays / m_iRaysPerSection;
	//	}

	if (bake_begin_function) {
		int progress_range = m_iHeight;
		bake_begin_function(progress_range);
	}


	m_iNumTests = 0;

	for (int y=0; y<m_iHeight; y++)
	{
		if ((y % 10) == 0)
		{
			//print_line("\tTexels line " + itos(y));
			//OS::get_singleton()->delay_usec(1);

			if (bake_step_function) {
				m_bCancel = bake_step_function(y, String("Process Texels: ") + " (" + itos(y) + ")");
				if (m_bCancel)
				{
					if (bake_end_function)
						bake_end_function();
					return;
				}
			}
		}

		for (int x=0; x<m_iWidth; x++)
		{
			ProcessTexel(x, y);
		}
	}



	//	m_iNumTests /= (m_iHeight * m_iWidth);
	print_line("num tests : " + itos(m_iNumTests));

	for (int b=0; b<m_AdjustedSettings.m_Backward_NumBounces; b++)
	{
		ProcessTexels_Bounce();
	}

	if (bake_end_function) {
		bake_end_function();
	}
}


void LightMapper::ProcessTexel_Light(int light_id, const Vector3 &ptDest, const Vector3 &ptNormal, FColor &color, uint32_t tri_ignore)
{
	const LLight &light = m_Lights[light_id];

	Ray r;
	//	r.o = Vector3(0, 5, 0);

	//	float range = light.scale.x;
	//	const float range = 2.0f;

	// the power should depend on the volume, with 1x1x1 being normal power
	//	float power = light.scale.x * light.scale.y * light.scale.z;
	float power = light.energy;
	power *= m_Settings_Backward_RayPower;

	int nSamples = m_AdjustedSettings.m_Backward_NumRays;

	// total light hitting texel
//	float fTotal = 0.0f;
	color.Set(0.0f);
	float total = 0.0f;


	// for a spotlight, we can cull completely in a lot of cases.
	if (light.type == LLight::LT_SPOT)
	{
		r.o = light.spot_emanation_point;
		r.d = ptDest - r.o;
		r.d.normalize();
		float dot = r.d.dot(light.dir);
		//float angle = safe_acosf(dot);
		//if (angle >= light.spot_angle_radians)

		dot -= light.spot_dot_max;

		if (dot <= 0.0f)
			return;
	}


	// each ray
	for (int n=0; n<nSamples; n++)
	{
		r.o = light.pos;

		// allow falloff for cones
		float multiplier = 1.0f;

		switch (light.type)
		{
		case LLight::LT_SPOT:
			{
				// source
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;
				r.o += offset;

				// offset from origin to destination texel
				r.d = ptDest - r.o;
				r.d.normalize();

				float dot = r.d.dot(light.dir);
				//float angle = safe_acosf(dot);
				//if (angle >= light.spot_angle_radians)

				dot -= light.spot_dot_max;

				if (dot <= 0.0f)
					continue;

				dot *= 1.0f / (1.0f - light.spot_dot_max);
				multiplier = dot * dot;
				multiplier *= multiplier;



				// direction
				//				r.d = light.dir;
				//				float spot_ball = 0.2f;
				//				float x = Math::random(-spot_ball, spot_ball);
				//				float y = Math::random(-spot_ball, spot_ball);
				//				float z = Math::random(-spot_ball, spot_ball);
				//				r.d += Vector3(x, y, z);
				//				r.d.normalize();
			}
			break;
		default:
			{
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;
				r.o += offset;

				// offset from origin to destination texel
				r.d = ptDest - r.o;
				r.d.normalize();

			}
			break;
		}


		// collision detect
		float u, v, w, t;

		m_Scene.m_Tracer.m_bUseSDF = true;
		int tri = m_Scene.FindIntersect_Ray(r, u, v, w, t, nullptr, m_iNumTests);
		//		m_Scene.m_Tracer.m_bUseSDF = false;
		//		int tri2 = m_Scene.IntersectRay(r, u, v, w, t, m_iNumTests);
		//		if (tri != tri2)
		//		{
		//			// repeat SDF version
		//			m_Scene.m_Tracer.m_bUseSDF = true;
		//			int tri = m_Scene.IntersectRay(r, u, v, w, t, m_iNumTests);
		//		}

		// nothing hit
		if ((tri == -1) || (tri == (int) tri_ignore))
		{
			// for backward tracing, first pass, this is a special case, because we DO
			// take account of distance to the light, and normal, in order to simulate the effects
			// of the likelihood of 'catching' a ray. In forward tracing this happens by magic.
			float dist = (r.o - ptDest).length();
			float local_power = power * InverseSquareDropoff(dist);

			// take into account normal
			float dot = r.d.dot(ptNormal);
			dot = fabs(dot);

			local_power *= dot;

			// cone falloff
			local_power *= multiplier;

			// albedo
			total += local_power;
//			FColor col = light.color * local_power;

//			color += col;
			//fTotal += local_power;
		}
	}


	color = light.color * total;
	// save in the texel
	//return fTotal;
}


FColor LightMapper::ProcessTexel_Bounce(int x, int y)
{
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

	int nSamples = m_AdjustedSettings.m_Backward_NumBounceRays;
	for (int n=0; n<nSamples; n++)
	{
		// bounce

		// first dot
		Ray r;

		// SLIDING
		//			Vector3 temp = r.d.cross(norm);
		//			new_ray.d = norm.cross(temp);

		// BOUNCING - mirror
		//new_ray.d = r.d - (2.0f * (dot * norm));

		// random hemisphere
		RandomUnitDir(r.d);

		// compare direction to normal, if opposite, flip it
		if (r.d.dot(norm) < 0.0f)
			r.d = -r.d;

		// add a little epsilon to prevent self intersection
		r.o = pos + (norm * 0.01f);
		//ProcessRay(new_ray, depth+1, power * 0.4f);

		// collision detect
		//r.d.normalize();
		float u, v, w, t;
		int tri_hit = m_Scene.FindIntersect_Ray(r, u, v, w, t, nullptr, m_iNumTests);

		// nothing hit
		if ((tri_hit != -1) && (tri_hit != (int) tri_source))
		{
			// look up the UV of the tri hit
			Vector2 uvs;
			m_Scene.FindUVsBarycentric(tri_hit, uvs, u, v, w);

			// find texel
			int dx = (uvs.x * m_iWidth); // round?
			int dy = (uvs.y * m_iHeight);

			if (m_Image_L.IsWithin(dx, dy))
			{
				// the contribution is the luminosity at that spot and the albedo
				Color albedo;
				m_Scene.FindPrimaryTextureColors(tri_hit, Vector3(u, v, w), albedo);
				FColor falbedo;
				falbedo.Set(albedo);

				total += (m_Image_L.GetItem(dx, dy) * falbedo);
			}

		}
	}

	return total / nSamples;
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


void LightMapper::ProcessTexel(int tx, int ty)
{
	// find triangle
	uint32_t tri = *m_Image_ID_p1.Get(tx, ty);
	if (!tri)
		return;
	tri--; // plus one based

	// barycentric
	const Vector3 &bary = *m_Image_Barycentric.Get(tx, ty);

	Vector3 pos;
	m_Scene.m_Tris[tri].InterpolateBarycentric(pos, bary.x, bary.y, bary.z);

	Vector3 normal;
	m_Scene.m_TriNormals[tri].InterpolateBarycentric(normal, bary.x, bary.y, bary.z);


	//Vector2i tex_uv = Vector2i(x, y);

	// could be off the image
	FColor * pTexel = m_Image_L.Get(tx, ty);
	if (!pTexel)
		return;

	FColor temp;
	for (int l=0; l<m_Lights.size(); l++)
	{
		ProcessTexel_Light(l, pos, normal, temp, tri);
		*pTexel += temp;
	}

}

void LightMapper::ProcessRay(LM::Ray r, int depth, float power, int dest_tri_id, const Vector2i * pUV)
{
	// unlikely
	if (r.d.x == 0.0f && r.d.y == 0.0f && r.d.z == 0.0f)
		return;

	// test
	//	r.d = Vector3(0, -1, 0);
	//	r.d = Vector3(-2.87, -5.0 + 0.226, 4.076);

	r.d.normalize();
	float u, v, w, t;
	int tri = m_Scene.FindIntersect_Ray(r, u, v, w, t, nullptr, m_iNumTests);

	// nothing hit
	if (tri == -1)
		return;

	// convert barycentric to uv coords in the lightmap
	Vector2 uv;
	m_Scene.FindUVsBarycentric(tri, uv, u, v, w);
	//	m_UVTris[tri].FindUVBarycentric(uvs, u, v, w);

	// texel address
	int tx = uv.x * m_iWidth;
	int ty = uv.y * m_iHeight;

	// override?
	if (pUV && tri == dest_tri_id)
	{
		tx = pUV->x;
		ty = pUV->y;
	}

	// could be off the image
	float * pf = &m_Image_L.Get(tx, ty)->r;
	if (!pf)
		return;

	// scale according to distance
	t /= 10.0f;
	t = 1.0f - t;
	if (t < 0.0f)
		t = 0.0f;
	t *= 2.0f;

	t = power;
	//	if (t > *pf)

	//	if (depth > 0)
	*pf += t;

	// bounce and lower power

	if (depth < m_Settings_Forward_NumBounces)
	{
		Vector3 pos;
		const Tri &triangle = m_Scene.m_Tris[tri];
		triangle.InterpolateBarycentric(pos, u, v, w);

		Vector3 norm;
		const Tri &triangle_normal = m_Scene.m_TriNormals[tri];
		triangle_normal.InterpolateBarycentric(norm, u, v, w);
		norm.normalize();

		// first dot
		float dot = norm.dot(r.d);
		if (dot <= 0.0f)
		{

			Ray new_ray;

			// SLIDING
			//			Vector3 temp = r.d.cross(norm);
			//			new_ray.d = norm.cross(temp);

			// BOUNCING - mirror
			Vector3 mirror_dir = r.d - (2.0f * (dot * norm));

			// random hemisphere
			const float range = 1.0f;
			Vector3 hemi_dir;
			while (true)
			{
				hemi_dir.x = Math::random(-range, range);
				hemi_dir.y = Math::random(-range, range);
				hemi_dir.z = Math::random(-range, range);

				float sl = hemi_dir.length_squared();
				if (sl > 0.0001f)
				{
					break;
				}
			}
			// compare direction to normal, if opposite, flip it
			if (hemi_dir.dot(norm) < 0.0f)
				hemi_dir = -hemi_dir;

			new_ray.d = hemi_dir.linear_interpolate(mirror_dir, m_Settings_Forward_BounceDirectionality);

			new_ray.o = pos + (norm * 0.01f);
			ProcessRay(new_ray, depth+1, power * m_Settings_Forward_BouncePower);
		} // in opposite directions
	}

}

void LightMapper::ProcessEmissionTris()
{
	int num_sections = m_iNumRays / m_iRaysPerSection;

	if (!num_sections)
		num_sections = 1;

	float fraction = 1.0f / num_sections;

	if (bake_begin_function) {
		bake_begin_function(num_sections);
	}

	for (int s=0; s<num_sections; s++)
	{
		if ((s % 1) == 0)
		{
			if (bake_step_function)
			{
				m_bCancel = bake_step_function(s, String("Process Emission Section: ") + " (" + itos(s) + ")");
				if (m_bCancel)
				{
					if (bake_end_function)
						bake_end_function();
					return;
				}
			}
		}

		ProcessEmissionTris_Section(fraction);

		for (int b=0; b<m_AdjustedSettings.m_Forward_NumBounces+1; b++)
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

void LightMapper::ProcessEmissionTris_Section(float fraction_of_total)
{
	for (int n=0; n<m_Scene.m_EmissionTris.size(); n++)
	{
		ProcessEmissionTri(n, fraction_of_total);
	}
}


void LightMapper::ProcessEmissionTri(int etri_id, float fraction_of_total)
{
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
	float rays_per_unit_area = m_iNumRays * m_AdjustedSettings.m_Forward_Emission_Density  * 0.12f * 0.05f;
	int nSamples = etri.area * rays_per_unit_area * fraction_of_total;

	// nSamples may be zero incorrectly for small triangles, maybe we need to adjust for this
	// NYI


	for (int s=0; s<nSamples; s++)
	{
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

		RayBank_RequestNewRay(ray, m_Settings_Forward_NumBounces + 1, fcol);

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

		FColor * pTexelCol = m_Image_L.Get(tx, ty);

//		fcol.Set(emission_color * 0.5f);
		*pTexelCol += fcol;

	}
}


void LightMapper::ProcessLights()
{
	//	const int rays_per_section = 1024 * 16;

	int num_sections = m_iNumRays / m_iRaysPerSection;

	if (bake_begin_function) {
		bake_begin_function(num_sections);
	}


	for (int n=0; n<m_Lights.size(); n++)
	{
		for (int s=0; s<num_sections; s++)
		{
			// double check the voxels are clear
#ifdef DEBUG_ENABLED
			//RayBank_CheckVoxelsClear();
#endif

			if ((s % 1) == 0)
			{
				if (bake_step_function)
				{
					m_bCancel = bake_step_function(s, String("Process Light Section: ") + " (" + itos(s) + ")");
					if (m_bCancel)
					{
						if (bake_end_function)
							bake_end_function();
						return;
					}
				}
			}


			ProcessLight(n, m_iRaysPerSection);

			for (int b=0; b<m_AdjustedSettings.m_Forward_NumBounces+1; b++)
			{
				RayBank_Process();
				RayBank_Flush();
			} // for bounce
		} // for section

		// left over
		{
			int num_leftover = m_iNumRays - (num_sections * m_iRaysPerSection);
			ProcessLight(n, num_leftover);

			for (int b=0; b<m_AdjustedSettings.m_Forward_NumBounces+1; b++)
			{
				RayBank_Process();
				RayBank_Flush();
			} // for bounce
		}

	} // for light

	if (bake_end_function) {
		bake_end_function();
	}

}

void LightMapper::ProcessLight(int light_id, int num_rays)
{
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
	float power = m_Settings_Forward_RayPower;

	if (light.indirect_energy > 0.001f)
		power *= 1.0f / light.indirect_energy;

	// for directional, we need a load more rays for it to work well - it is expensive
	if (light.type == LLight::LT_DIRECTIONAL)
	{
		num_rays *= 2;
//		float area = light.dl_tangent_range * light.dl_bitangent_range;
//		num_rays = num_rays * area;
		// we will increase the power as well, because daylight more powerful than light bulbs typically.
		power *= 4.0f;
	}


	FColor light_color = light.color * power;

	for (int n=0; n<num_rays; n++)
	{
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

		switch (light.type)
		{
		case LLight::LT_DIRECTIONAL:
			{
				// use the precalculated source plane stored in the llight
				r.o = light.dl_plane_pt;
				r.o += light.dl_tangent * Math::random(0.0f, light.dl_tangent_range);
				r.o += light.dl_bitangent * Math::random(0.0f, light.dl_bitangent_range);


//				Vector3 offset;
//				RandomUnitDir(offset);
//				offset *= light.scale;
//				r.o += offset;

//				r.d = light.dir;
				RandomUnitDir(r.d);
				r.d *= light.scale;

				// must point down - reverse hemisphere if pointing up
				if (light.dir.dot(r.d) < 0.0f)
				{
					r.d = -r.d;
				}

				//r.d = light.dir.linear_interpolate(r.d, fract);
				r.d += (light.dir * 2.0f);

//				r.d += (light.dir * 3.0f);
				r.d.normalize();
			}
			break;
		case LLight::LT_SPOT:
			{
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;
				r.o += offset;

				r.d = light.dir;


				// random axis
				Vector3 axis;
				RandomAxis(axis);

				float falloff_start = 0.5f;// this could be adjustable;

				float ang_max = light.spot_angle_radians;
				float ang_falloff = ang_max * falloff_start;
				float ang_falloff_range = ang_max - ang_falloff;

				// random angle giving equal circle distribution
				float a = Math::random(0.0f, 1.0f);

				// below a certain proportion are the central, outside is falloff
				a *= 2.0f;
				if (a > 1.0f)
				{
					// falloff
					a -= 1.0f;
					//angle /= 3.0f;

					//a *= a;

					a = ang_falloff + (a * ang_falloff_range);
				}
				else
				{
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
			}
			break;
		default:
			{
				Vector3 offset;
				RandomUnitDir(offset);
				offset *= light.scale;
				r.o += offset;

				RandomUnitDir(r.d);
			}
			break;
		}
		//r.d.normalize();

		RayBank_RequestNewRay(r, m_Settings_Forward_NumBounces + 1, light_color, 0);

		//		ProcessRay(r, 0, power);
	}
}

} // namespace
