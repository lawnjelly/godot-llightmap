#include "llighttracer.h"
#include "llightscene.h"
#include "core/os/os.h"



using namespace LM;

void LightTracer::Reset()
{
	m_Voxels.clear(true);
	m_VoxelBounds.clear(true);
	m_BFTrisHit.Blank();
	m_iNumTris = 0;


}


void LightTracer::Create(const LightScene &scene, int voxel_density)
{
	m_bUseSDF = true;

	m_pScene = &scene;
	m_iNumTris = m_pScene->m_Tris.size();

	CalculateWorldBound();
	CalculateVoxelDims(voxel_density);

	m_iNumVoxels = m_Dims.x * m_Dims.y * m_Dims.z;
	m_DimsXTimesY = m_Dims.x * m_Dims.y;

	m_Voxels.resize(m_iNumVoxels);
	m_VoxelBounds.resize(m_iNumVoxels);
	m_BFTrisHit.Create(m_iNumTris);

	m_VoxelSize.x = m_SceneWorldBound_expanded.size.x / m_Dims.x;
	m_VoxelSize.y = m_SceneWorldBound_expanded.size.y / m_Dims.y;
	m_VoxelSize.z = m_SceneWorldBound_expanded.size.z / m_Dims.z;

	// fill the bounds
	AABB aabb;
	aabb.size = m_VoxelSize;
	int count = 0;
	for (int z=0; z<m_Dims.z; z++)
	{
		aabb.position.z = m_SceneWorldBound_expanded.position.z + (z * m_VoxelSize.z);
		for (int y=0; y<m_Dims.y; y++)
		{
			aabb.position.y = m_SceneWorldBound_expanded.position.y + (y * m_VoxelSize.y);
			for (int x=0; x<m_Dims.x; x++)
			{
				aabb.position.x = m_SceneWorldBound_expanded.position.x + (x * m_VoxelSize.x);
				m_VoxelBounds[count++] = aabb;
			} // for x
		} // for y
	} // for z

	FillVoxels();
}

void LightTracer::FindNearestVoxel(const Vector3 &ptWorld, Vec3i &ptVoxel) const
{
	Vector3 pt = ptWorld;
	pt -= m_SceneWorldBound_expanded.position;
	pt.x /= m_VoxelSize.x;
	pt.y /= m_VoxelSize.y;
	pt.z /= m_VoxelSize.z;

	ptVoxel.Set(pt.x, pt.y, pt.z);
}

void LightTracer::GetDistanceInVoxels(float dist, Vec3i &ptVoxelDist) const
{
	// note this will screw up with zero voxel size
	if ((m_VoxelSize.x == 0.0f) || (m_VoxelSize.y == 0.0f) || (m_VoxelSize.z == 0.0f))
	{
		ptVoxelDist.Set(0, 0, 0);
		return;
	}

	ptVoxelDist.x = (dist / m_VoxelSize.x)+1;//+1;
	ptVoxelDist.y = (dist / m_VoxelSize.y)+1;//+1;
	ptVoxelDist.z = (dist / m_VoxelSize.z)+1;//+1;
}

// ray translated to voxel space
bool LightTracer::RayTrace_Start(Ray ray, Ray &voxel_ray, Vec3i &start_voxel)
{
	// if tracing from outside, try to trace to the edge of the world bound
	if (!m_SceneWorldBound_expanded.has_point(ray.o))
	{
		Vector3 clip;
		//if (!IntersectRayAABB(ray, m_SceneWorldBound_contracted, clip))
		if (!IntersectRayAABB(ray, m_SceneWorldBound_expanded, clip))
			return false;

		// does hit the world bound
		ray.o = clip;
	}

//	m_BFTrisHit.Blank();

	voxel_ray.o = ray.o - m_SceneWorldBound_expanded.position;
	voxel_ray.o.x /= m_VoxelSize.x;
	voxel_ray.o.y /= m_VoxelSize.y;
	voxel_ray.o.z /= m_VoxelSize.z;
	voxel_ray.d.x = ray.d.x / m_VoxelSize.x;
	voxel_ray.d.y = ray.d.y / m_VoxelSize.y;
	voxel_ray.d.z = ray.d.z / m_VoxelSize.z;
	voxel_ray.d.normalize();

	start_voxel.x = voxel_ray.o.x;
	start_voxel.y = voxel_ray.o.y;
	start_voxel.z = voxel_ray.o.z;

//	m_TriHits.clear();

	// out of bounds?
	bool within = VoxelWithinBounds(start_voxel);

	// instead of trying to calculate these on the fly with intersection
	// tests we can use simple linear addition to calculate them all quickly.
	if (within)
	{
		//PrecalcRayCuttingPoints(voxel_ray);
	}
//	if (within)
//		DebugCheckPointInVoxel(voxel_ray.o, start_voxel);

	// debug check the voxel number and bounding box are correct

//	if (m_bUseSDF)
//		print_line("Ray start");


	return within;
}

void LightTracer::DebugCheckWorldPointInVoxel(Vector3 pt, const Vec3i &ptVoxel)
{
	int iVoxelNum = GetVoxelNum(ptVoxel);
	AABB bb = m_VoxelBounds[iVoxelNum];
	bb.grow_by(0.01f);
	assert (bb.has_point(pt));
}


//bool LightTracer::RayTrace(const Ray &ray_orig, Ray &ray_out, Vec3i &ptVoxel)
const Voxel * LightTracer::RayTrace(const Ray &ray_orig, Ray &ray_out, Vec3i &ptVoxel)
{
	//m_TriHits.clear();

#ifdef LIGHTTRACER_IGNORE_VOXELS
	for (int n=0; n<m_iNumTris; n++)
	{
		m_TriHits.push_back(n);
	}
	return true;
#endif

	if (!VoxelWithinBounds(ptVoxel))
		return 0;



	// debug check
	DebugCheckLocalPointInVoxel(ray_orig.o, ptVoxel);

	// add the tris in this voxel
	int iVoxelNum = GetVoxelNum(ptVoxel);
	const Voxel &vox = m_Voxels[iVoxelNum];
	const Voxel * pCurrVoxel = &vox;

//	const AABB &bb = m_VoxelBounds[iVoxelNum];
//	print_line("Checking Voxel " + ptVoxel.ToString() + " bound " + String(bb));


	if (!m_bSIMD)
	{
		/*
		for (int n=0; n<vox.m_TriIDs.size(); n++)
		{
			unsigned int id = vox.m_TriIDs[n];

			// check bitfield, the tri may already have been added
			if (m_BFTrisHit.CheckAndSet(id))
				m_TriHits.push_back(id);
		}
		*/
	}

//#define LLIGHT_USE_SDF
#ifdef LLIGHT_USE_SDF
	int sdf = vox.m_SDF;

	sdf-= 1;
	if (sdf >= 0)
	{
		;
	}
	else
	{
		sdf = 0;
	}
//	sdf = 0;

	if (!m_bUseSDF)
	{
		sdf = 0;
	}

//		print_line("\tvoxel (" + itos(ptVoxel.x) + " " + itos(ptVoxel.y) + " " + itos(ptVoxel.z) +
//		"),	SDF " + itos(vox.m_SDF) + " num_tris " + itos (vox.m_iNumTriangles) + ", jump " + itos(sdf));



	int change_withsdf = sdf + 1;

//	sdf = 0;
	// SDF is now how many planes to cross.

	// PLANES are in INTEGER space (voxel space)
	// attempt to cross out of planes
//	const AABB &bb = m_VoxelBounds[iVoxelNum];
	Vector3 mins = Vector3(ptVoxel.x-sdf, ptVoxel.y-sdf, ptVoxel.z-sdf);
	Vector3 maxs = mins + Vector3(change_withsdf, change_withsdf, change_withsdf);

	// add a bit of bias to ensure we cross the boundary into another voxel and don't get float error
	// (note we may also have to catch triangles slightly outside the voxel to compensate for this)
	const float plane_bias = 0.001f;
	Vector3 plane_bias3 = Vector3(plane_bias, plane_bias, plane_bias);
	mins -= plane_bias3;
	maxs += plane_bias3;
#else
	Vector3 mins = Vector3(ptVoxel.x, ptVoxel.y, ptVoxel.z);
	Vector3 maxs = mins + Vector3(1, 1, 1);
#endif


//	const Vector3 &mins = bb.position;
//	Vector3 maxs = bb.position + bb.size;
	const Vector3 &dir = ray_orig.d;

	// the 3 intersection points
	Vector3 ptIntersect[3];
	float nearest_hit = FLT_MAX;
	int nearest_hit_plane = -1;

	//Vector3 ptBias;

	// planes from constants
	if (dir.x >= 0.0f)
	{
		ptIntersect[0].x = maxs.x;
		IntersectAAPlane(ray_orig, 0, ptIntersect[0], nearest_hit, 0, nearest_hit_plane);
		//ptBias.x = 0.5f;
	}
	else
	{
		ptIntersect[0].x = mins.x;
		IntersectAAPlane(ray_orig, 0, ptIntersect[0], nearest_hit, 1, nearest_hit_plane);
		//ptBias.x = -0.5f;
	}

	if (dir.y >= 0.0f)
	{
		ptIntersect[1].y = maxs.y;
		IntersectAAPlane(ray_orig, 1, ptIntersect[1], nearest_hit, 2, nearest_hit_plane);
		//ptBias.y = 0.5f;
	}
	else
	{
		ptIntersect[1].y = mins.y;
		IntersectAAPlane(ray_orig, 1, ptIntersect[1], nearest_hit, 3, nearest_hit_plane);
		//ptBias.y = -0.5f;
	}

	if (dir.z >= 0.0f)
	{
		ptIntersect[2].z = maxs.z;
		IntersectAAPlane(ray_orig, 2, ptIntersect[2], nearest_hit, 4, nearest_hit_plane);
		//ptBias.z = 0.5f;
	}
	else
	{
		ptIntersect[2].z = mins.z;
		IntersectAAPlane(ray_orig, 2, ptIntersect[2], nearest_hit, 5, nearest_hit_plane);
		//ptBias.z = -0.5f;
	}

	// ray out
	ray_out.d = ray_orig.d;
	ray_out.o = ptIntersect[nearest_hit_plane/2];

#ifdef LLIGHT_USE_SDF
	const Vector3 out = ray_out.o;// + ptBias;
	ptVoxel.x = floorf(out.x);
	ptVoxel.y = floorf(out.y);
	ptVoxel.z = floorf(out.z);
#else


	switch (nearest_hit_plane)
	{
	case 0:
		ptVoxel.x += 1;
		break;
	case 1:
		ptVoxel.x -= 1;
		break;
	case 2:
		ptVoxel.y += 1;
		break;
	case 3:
		ptVoxel.y -= 1;
		break;
	case 4:
		ptVoxel.z += 1;
		break;
	case 5:
		ptVoxel.z -= 1;
		break;
	default:
		assert (0 && "LightTracer::RayTrace");
		break;
	}

#endif
	return pCurrVoxel;
}


void LightTracer::Debug_SaveSDF()
{
	int width = m_Dims.x;
	int height = m_Dims.z;

	Ref<Image> image = memnew(Image(width, height, false, Image::FORMAT_RGBA8));
	image->lock();
	int y = m_Dims.y-1;

	for (int z=0; z<m_Dims.z; z++)
	{
//		for (int y=0; y<m_Dims.y; y++)
//		{
			for (int x=0; x<m_Dims.x; x++)
			{
				Vec3i pt;
				pt.Set(x, y, z);

				int sdf = GetVoxel(pt).m_SDF;
				float f = sdf * 0.25f;
				image->set_pixel(x, z, Color(f, f, f, 1.0f));
			}
//		}
	}

	image->unlock();
	image->save_png("sdf.png");

}

void LightTracer::CalculateSDF()
{
	return;

	print_line("Calculating SDF");

	// look at the surrounding neighbours. We should be at a minimum, the lowest neighbour +1
	int iters = m_Dims.x;
	if (m_Dims.y > iters) iters = m_Dims.y;
	if (m_Dims.z > iters) iters = m_Dims.z;

	for (int i=0; i<iters; i++)
	{

		// SDF is seeded with zero in the filled voxels
		for (int z=0; z<m_Dims.z; z++)
		{
			for (int y=0; y<m_Dims.y; y++)
			{
				for (int x=0; x<m_Dims.x; x++)
				{
					Vec3i pt;
					pt.Set(x, y, z);
					CalculateSDF_Voxel(pt);
				}
			}
		}

	} // for iteration

	// debug print sdf
//	Debug_SaveSDF();

}

void LightTracer::CalculateSDF_AssessNeighbour(const Vec3i &pt, unsigned int &min_SDF)
{
	// on map?
	if (!VoxelWithinBounds(pt))
		return;

	const Voxel &vox = GetVoxel(pt);
	if (vox.m_SDF < min_SDF)
		min_SDF = vox.m_SDF;
}

void LightTracer::CalculateSDF_Voxel(const Vec3i &ptCentre)
{
	Voxel &vox = GetVoxel(ptCentre);


	unsigned int lowest = UINT_MAX-1;

		for (int nz=-1; nz<=1; nz++)
		{
			for (int ny=-1; ny<=1; ny++)
			{
				for (int nx=-1; nx<=1; nx++)
				{
					if ((nx == 0) && (ny == 0) && (nz == 0))
					{
					}
					else
					{
						Vec3i pt;
						pt.Set(ptCentre.x + nx, ptCentre.y + ny, ptCentre.z + nz);
						CalculateSDF_AssessNeighbour(pt, lowest);
					}
				} // for nx
			} // for ny
		} // for nz

	lowest += 1;
	if (vox.m_SDF > lowest)
		vox.m_SDF = lowest;
}


void LightTracer::FillVoxels()
{
	print_line("FillVoxels : Num AABBs " + itos(m_pScene->m_TriPos_aabbs.size()));
	print_line("NumTris " + itos (m_iNumTris));

	int count = 0;
	for (int z=0; z<m_Dims.z; z++)
	{
		for (int y=0; y<m_Dims.y; y++)
		{
			for (int x=0; x<m_Dims.x; x++)
			{
				Voxel &vox = m_Voxels[count];
				vox.Reset();
				AABB aabb = m_VoxelBounds[count++];

				// expand the aabb just a little to account for float error
				aabb.grow_by(0.001f);


//				if ((z == 1) && (y == 1) && (x == 0))
//				{
//					print_line("AABB voxel is " + String(Variant(aabb)));
//				}

				// find all tris within
				for (int t=0; t<m_iNumTris; t++)
				{
					if (m_pScene->m_TriPos_aabbs[t].intersects(aabb))
					{
						// add tri to voxel
						//vox.m_TriIDs.push_back(t);
						vox.AddTriangle(m_pScene->m_Tris_EdgeForm[t], t);

//						if ((z == 1) && (y == 1) && (x == 0))
//						{
//							print_line("tri " + itos (t) + " AABB " + String(Variant(m_pScene->m_TriPos_aabbs[t])));
//						}
					}
				}
				vox.Finalize();

//				if ((z == 1) && (y == 1) && (x == 0))
//				{
//					print_line("\tvoxel line x " + itos(x) + ", " + itos(vox.m_TriIDs.size()) + " tris.");
//				}
			} // for x
		} // for y
	} // for z

	CalculateSDF();
}

void LightTracer::CalculateVoxelDims(int voxel_density)
{
	const AABB &aabb = m_SceneWorldBound_expanded;
	float max_length = aabb.get_longest_axis_size();

	m_Dims.x = ((aabb.size.x / max_length) * voxel_density) + 0.01f;
	m_Dims.y = ((aabb.size.y / max_length) * voxel_density) + 0.01f;
	m_Dims.z = ((aabb.size.z / max_length) * voxel_density) + 0.01f;

	// minimum of 1
	m_Dims.x = MAX(m_Dims.x, 1);
	m_Dims.y = MAX(m_Dims.y, 1);
	m_Dims.z = MAX(m_Dims.z, 1);

	print_line("voxels dims : " + itos(m_Dims.x) + ", " + itos(m_Dims.y) + ", " + itos (m_Dims.z));
}

void LightTracer::CalculateWorldBound()
{
	if (!m_iNumTris)
		return;

	AABB &aabb = m_SceneWorldBound_expanded;
	aabb.position = m_pScene->m_Tris[0].pos[0];
	aabb.size = Vector3(0, 0, 0);

	for (int n=0; n<m_iNumTris; n++)
	{
		const Tri &tri = m_pScene->m_Tris[n];
		aabb.expand_to(tri.pos[0]);
		aabb.expand_to(tri.pos[1]);
		aabb.expand_to(tri.pos[2]);
	}

	// exact
	m_SceneWorldBound_contracted = aabb;

	// expanded

	// it is CRUCIAL that the expansion here is more than the push in provided
	// by LightTracer::IntersectRayAABB
	// otherwise triangles at the very edges of the world will be missed by the ray tracing.

	aabb.grow_by(LIGHTTRACER_EXPANDED_BOUND);

	m_SceneWorldBound_mid = m_SceneWorldBound_contracted;
	m_SceneWorldBound_mid.grow_by(LIGHTTRACER_HALF_EXPANSION);
}

bool LightTracer::IntersectRayAABB(const Ray &ray, const AABB &aabb, Vector3 &ptInter)
{
	// the 3 intersection points
	const Vector3 &mins = aabb.position;
	Vector3 maxs = aabb.position + aabb.size;


	Vector3 ptIntersect[3];
	float nearest_hit = FLT_MAX;
	int nearest_hit_plane = -1;
	const Vector3 &dir = ray.d;

	// planes from constants
	if (dir.x <= 0.0f)
	{
		ptIntersect[0].x = maxs.x;
		IntersectAAPlane_OnlyWithinAABB(aabb, ray, 0, ptIntersect[0], nearest_hit, 0, nearest_hit_plane);
	}
	else
	{
		ptIntersect[0].x = mins.x;
		IntersectAAPlane_OnlyWithinAABB(aabb, ray, 0, ptIntersect[0], nearest_hit, 1, nearest_hit_plane);
	}

	if (dir.y <= 0.0f)
	{
		ptIntersect[1].y = maxs.y;
		IntersectAAPlane_OnlyWithinAABB(aabb, ray, 1, ptIntersect[1], nearest_hit, 2, nearest_hit_plane);
	}
	else
	{
		ptIntersect[1].y = mins.y;
		IntersectAAPlane_OnlyWithinAABB(aabb, ray, 1, ptIntersect[1], nearest_hit, 3, nearest_hit_plane);
	}

	if (dir.z <= 0.0f)
	{
		ptIntersect[2].z = maxs.z;
		IntersectAAPlane_OnlyWithinAABB(aabb, ray, 2, ptIntersect[2], nearest_hit, 4, nearest_hit_plane);
	}
	else
	{
		ptIntersect[2].z = mins.z;
		IntersectAAPlane_OnlyWithinAABB(aabb, ray, 2, ptIntersect[2], nearest_hit, 5, nearest_hit_plane);
	}

	ptInter = ptIntersect[nearest_hit_plane/2];

	if (nearest_hit_plane == -1)
		return false;

	// recalculate intersect using distance plus epsilon
	float nearest_length = sqrtf(nearest_hit);

	// this epsilon MUST be less than the world expansion in LightTracer::CalculateWorldBound
	ptInter = ray.o + (ray.d * (nearest_length + LIGHTTRACER_HALF_EXPANSION));

	if (aabb.has_point(ptInter))
		return true;

	return false;
}
