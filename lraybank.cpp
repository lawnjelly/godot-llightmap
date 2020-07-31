#include "lraybank.h"
#include "llightscene.h"
#include "core/os/threaded_array_processor.h"

namespace LM {


void RayBank::RayBank_Data::Swap()
{
	if (m_MapRead)
	{
		m_MapRead = 0;
		m_MapWrite = 1;
	}
	else
	{
		m_MapRead = 1;
		m_MapWrite = 0;
	}
}


void RayBank::RayBank_Reset()
{
	m_Data_RB.m_Voxels[0].clear(true);
	m_Data_RB.m_Voxels[1].clear(true);

}

void RayBank::RayBank_Create()
{
	m_Data_RB.m_MapWrite = 0;
	m_Data_RB.m_MapRead = 1;
	int num_voxels = GetTracer().m_iNumVoxels;
	m_Data_RB.m_Voxels[0].resize(num_voxels);
	m_Data_RB.m_Voxels[1].resize(num_voxels);
}


// either we know the start voxel or we find it during this routine (or it doesn't cut the world)
FRay * RayBank::RayBank_RequestNewRay(Ray ray, int num_rays_left, const FColor &col, const Vec3i * pStartVoxel)
{
	// if we don't know the start voxel
	Vec3i ptStartVoxel;
	if (!pStartVoxel)
	{
		pStartVoxel = &ptStartVoxel;

		// if tracing from outside, try to trace to the edge of the world bound
		if (!GetTracer().m_SceneWorldBound.has_point(ray.o))
		{
			Vector3 clip;

			// if the ray starts outside, and doesn't hit the world, the ray is invalid
			if (!GetTracer().IntersectRayAABB(ray, GetTracer().m_SceneWorldBound_contracted, clip))
				return 0;

			// does hit the world bound
			ray.o = clip;
		}

		const AABB &world_bound = GetTracer().m_SceneWorldBound;
		const Vector3 &voxel_size = GetTracer().m_VoxelSize;

		// ray origin should now be in the bound
		Vector3 o_voxelspace = ray.o - world_bound.position;
		o_voxelspace.x /= voxel_size.x;
		o_voxelspace.y /= voxel_size.y;
		o_voxelspace.z /= voxel_size.z;

		ptStartVoxel.x = o_voxelspace.x;
		ptStartVoxel.y = o_voxelspace.y;
		ptStartVoxel.z = o_voxelspace.z;
	}

	// check start voxel is within bound
	if (!GetTracer().VoxelWithinBounds(*pStartVoxel))
	{
		// should not happen?
		return 0;
	}

	RB_Voxel &vox = RayBank_GetVoxelWrite(*pStartVoxel);

	FRay * fray = vox.m_Rays.request();

	// should not happen
	if (!fray)
		return 0;

	fray->ray = ray;
	fray->hit.SetNoHit();
	fray->num_rays_left = num_rays_left;
	fray->color = col;

	return fray;
}


// can be used from several threads
void RayBank::RayBank_Process()
{
	// swap the write and read
	m_Data_RB.Swap();

	int nVoxels = GetTracer().m_iNumVoxels;

	//thread_process_array(width, this, &VoxelLightBaker::_lightmap_bake_point, &lightmap_ptr[i * width]);

	int nCores = OS::get_singleton()->get_processor_count();

	for (int v=0; v<nVoxels; v++)
	{
		RB_Voxel & vox = m_Data_RB.GetVoxels_Read()[v];

//		if (bake_step_function)
//		{
//			m_bCancel = bake_step_function(v, String("RayBank_Process: ") + " (" + itos(v) + ")");
//			if (m_bCancel)
//				return;
//		}


		int num_rays = vox.m_Rays.size();
		if (!num_rays)
			continue;

		m_pCurrentThreadVoxel = &vox;

		//const int section_size = 1024 * 96;
		int section_size = num_rays / nCores;

		int leftover_start = 0;

		// not worth doing multithread below a certain size
		// because of threading overhead
#define RAYBANK_USE_THREADING
#ifdef RAYBANK_USE_THREADING
		if (section_size >= 64)
		{
			int num_sections = num_rays / section_size;
			for (int s=0; s<num_sections; s++)
			{
				int section_start = s * section_size;


				thread_process_array(section_size, this, &RayBank::RayBank_ProcessRay_MT, section_start);

//				for (int n=0; n<section_size; n++)
//				{
//					RayBank_ProcessRay_MT(n, section_start);
//				}
			}

			leftover_start = num_sections * section_size;
		}
#endif

		// leftovers
		for (int n=leftover_start; n<num_rays; n++)
		{
			RayBank_ProcessRay_MT(n, 0);
		}


/*
		for (int n=0; n<vox.m_Rays.size(); n++)
		{
//			if ((n % 10000) == 0)
//			{
//				if (bake_step_function)
//				{
//					m_bCancel = bake_step_function(n, String("RayBank_Process: ") + " (" + itos(n) + ")");
//					if (m_bCancel)
//						return;
//				}
//			}

			RayBank_ProcessRay(n, vox);
		}
		*/
	}
}

void RayBank::RayBank_CheckVoxelsClear()
{
#ifdef DEBUG_ENABLED
	int nVoxels = GetTracer().m_iNumVoxels;
	LVector<RB_Voxel> &voxelsr = m_Data_RB.GetVoxels_Read();
	LVector<RB_Voxel> &voxelsw = m_Data_RB.GetVoxels_Write();

	for (int v=0; v<nVoxels; v++)
	{
		RB_Voxel & voxr = voxelsr[v];
		assert (voxr.m_Rays.size() == 0);
		RB_Voxel & voxw = voxelsw[v];
		assert (voxw.m_Rays.size() == 0);
	}
#endif
}

// flush ray results to the lightmap
void RayBank::RayBank_Flush()
{
	int nVoxels = GetTracer().m_iNumVoxels;
	LVector<RB_Voxel> &voxels = m_Data_RB.GetVoxels_Read();

	for (int v=0; v<nVoxels; v++)
	{
		RB_Voxel & vox = voxels[v];

		// save results to lightmap
		for (int n=0; n<vox.m_Rays.size(); n++)
		{
			RayBank_FlushRay(vox, n);
		}

		// delete rays
		vox.m_Rays.clear();
	}

	// swap the write and read
//	m_Data_RB.Swap();
}


void RayBank::RayBank_FlushRay(RB_Voxel &vox, int ray_id)
{
	const FRay &fray= vox.m_Rays[ray_id];

	// bounces first
	if (fray.num_rays_left)
	{
//		RayBank_RequestNewRay(fray.ray, fray.num_rays_left, fray.color * m_Settings_Forward_BouncePower, 0);
		RayBank_RequestNewRay(fray.ray, fray.num_rays_left, fray.bounce_color, 0);
	}

	// now write the hit to the lightmap
	const FHit &hit = fray.hit;
	if (hit.IsNoHit())
		return;

	FColor * pf = m_Image_L.Get(hit.tx, hit.ty);
#ifdef DEBUG_ENABLED
	assert (pf);
#endif
//	if (!pf)
//		return; // should never happen

	*pf += fray.color;
}

//void RayBank::RayBank_ProcessRay(uint32_t ray_id, RB_Voxel &vox)

void RayBank::RayBank_ProcessRay_MT(uint32_t ray_id, int start_ray)
{
	ray_id += start_ray;
	RB_Voxel &vox = *m_pCurrentThreadVoxel;

	FRay &fray= vox.m_Rays[ray_id];
	Ray r = fray.ray;

	// each evaluation
	fray.num_rays_left -= 1;

	// unlikely
	if (r.d.x == 0.0f && r.d.y == 0.0f && r.d.z == 0.0f)
	{
		fray.num_rays_left = 0;
		return;
	}

	r.d.normalize();
	float u, v, w, t;
	int tri = m_Scene.FindIntersect_Ray(r, u, v, w, t, nullptr, m_iNumTests);

	// nothing hit
	if (tri == -1)
	{
		fray.num_rays_left = 0;
		return;
	}

	// hit the back of a face? if so terminate ray
	Vector3 face_normal;
	const Tri &triangle_normal = m_Scene.m_TriNormals[tri];
	triangle_normal.InterpolateBarycentric(face_normal, u, v, w);
	face_normal.normalize();

	// first dot
	float dot = face_normal.dot(r.d);
	if (dot >= 0.0f)
	{
		fray.num_rays_left = 0;
		return;
	}

	// convert barycentric to uv coords in the lightmap
	Vector2 uv;
	m_Scene.FindUVsBarycentric(tri, uv, u, v, w);

	// texel address
	int tx = uv.x * m_iWidth;
	int ty = uv.y * m_iHeight;

	// override?
//	if (pUV && tri == dest_tri_id)
//	{
//		tx = pUV->x;
//		ty = pUV->y;
//	}

	// could be off the image
	if (!m_Image_L.IsWithin(tx, ty))
	{
		fray.num_rays_left = 0;
		return;
	}


	// register the hit
	FHit &hit = fray.hit;
	hit.tx = tx;
	hit.ty = ty;
//	hit.power = fray.power;
//	fray.num_hits += 1;

	// max hits?
//	if (fray.num_hits == FRay::FRAY_MAX_HITS)
//	{
//		RayBank_EndRay(fray);
//		return false;
//	}

	/*
	float * pf = m_Image_L.Get(tx, ty);
	if (!pf)
		return;

	// scale according to distance
	t /= 10.0f;
	t = 1.0f - t;
	if (t < 0.0f)
		t = 0.0f;
	t *= 2.0f;

	t = fray.power;

//	if (depth > 0)
	*pf += t;
	*/

	// bounce and lower power


	if (fray.num_rays_left)
	{
		Vector3 pos;
		const Tri &triangle = m_Scene.m_Tris[tri];
		triangle.InterpolateBarycentric(pos, u, v, w);

		// get the albedo etc
		Color albedo;
		m_Scene.FindPrimaryTextureColors(tri, Vector3(u, v, w), albedo);
		FColor falbedo;
		falbedo.Set(albedo);

		// test
		//fray.color = falbedo;

		// pre find the bounce color here
		fray.bounce_color = fray.color * falbedo * m_Settings_Forward_BouncePower;
//		fray.bounce_color = fray.color * m_Settings_Forward_BouncePower;

//		Vector3 norm;
//		const Tri &triangle_normal = m_Scene.m_TriNormals[tri];
//		triangle_normal.InterpolateBarycentric(norm, u, v, w);
//		face_normal.normalize();

		// first dot
//		float dot = face_normal.dot(r.d);
//		if (dot <= 0.0f)
		{

			Ray new_ray;

			// SLIDING
//			Vector3 temp = r.d.cross(norm);
//			new_ray.d = norm.cross(temp);

			// BOUNCING - mirror
			Vector3 mirror_dir = r.d - (2.0f * (dot * face_normal));

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
			if (hemi_dir.dot(face_normal) < 0.0f)
				hemi_dir = -hemi_dir;

			new_ray.d = hemi_dir.linear_interpolate(mirror_dir, m_Settings_Forward_BounceDirectionality);

			new_ray.o = pos + (face_normal * 0.01f);

			// copy the info to the existing fray
			fray.ray = new_ray;
			//fray.power *= m_Settings_Forward_BouncePower;

			return;
//			return true;
//			RayBank_RequestNewRay(new_ray, fray.num_rays_left, fray.power * m_Settings_Forward_BouncePower, 0);
		} // in opposite directions
//		else
//		{ // if normal in same direction as ray
//			fray.num_rays_left = 0;
//		}
	} // if there are bounces left

//	return false;
}

//void RayBank::RayBank_EndRay(FRay &fray)
//{
//	// mark fray as done
//	fray.num_rays_left = 0;
//}


} // namespace
