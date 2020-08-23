#pragma once

#include "lvector.h"
#include "lbitfield_dynamic.h"
#include "core/math/aabb.h"
#include "llighttypes.h"
#include "llighttests_simd.h"
#include <limits.h>

//#define LIGHTTRACER_IGNORE_VOXELS

#define LIGHTTRACER_EXPANDED_BOUND (0.2f)
#define LIGHTTRACER_HALF_EXPANSION (LIGHTTRACER_EXPANDED_BOUND / 2.0f)


namespace LM
{

class LightScene;

class Voxel
{
public:
	void Reset() {m_TriIDs.clear(); m_PackedTriangles.clear(); m_iNumTriangles = 0; m_SDF = UINT_MAX;}
	LVector<uint32_t> m_TriIDs;

	// a COPY of the triangles in SIMD format, edge form
	// contiguous in memory for faster testing
	LVector<PackedTriangles> m_PackedTriangles;
	int m_iNumTriangles;
	unsigned int m_SDF; // measured in voxels

	void AddTriangle(const Tri &tri, uint32_t tri_id)
	{
		m_TriIDs.push_back(tri_id);
		uint32_t packed = m_iNumTriangles / 4;
		uint32_t mod = m_iNumTriangles % 4;
		if (packed >= m_PackedTriangles.size())
		{
			PackedTriangles * pNew = m_PackedTriangles.request();
			pNew->Create();
		}

		// fill the relevant packed triangle
		PackedTriangles &tris = m_PackedTriangles[packed];
		tris.Set(mod, tri);
		m_iNumTriangles++;
	}
	void Finalize()
	{
		uint32_t packed = m_iNumTriangles / 4;
		uint32_t mod = m_iNumTriangles % 4;
		if (mod)
		{
			PackedTriangles &tris = m_PackedTriangles[packed];
			tris.Finalize(mod);
		}
		if (m_iNumTriangles)
			m_SDF = 0; // seed the SDF
	}
};

class LightTracer
{
public:
	friend class RayBank;

	void Reset();
//	void Create(const LightScene &scene, const Vec3i &voxel_dims);
	void Create(const LightScene &scene, int voxel_density);


	// translate a real world distance into a number of voxels in each direction
	// (this can be used to limit the distance in ray traces)
	void GetDistanceInVoxels(float dist, Vec3i &ptVoxelDist) const;

	bool RayTrace_Start(Ray ray, Ray &voxel_ray, Vec3i &start_voxel);
	const Voxel * RayTrace(const Ray &ray_orig, Ray &ray_out, Vec3i &ptVoxel);

	LVector<uint32_t> m_TriHits;

	bool m_bSIMD;
	bool m_bUseSDF;

	void FindNearestVoxel(const Vector3 &ptWorld, Vec3i &ptVoxel) const;
	const AABB &GetWorldBound_expanded() const {return m_SceneWorldBound_expanded;}
	const AABB &GetWorldBound_mid() const {return m_SceneWorldBound_mid;}

//	void ClampVoxelToBounds(Vec3i &v) const
//	{
//		v.x = CLAMP(v.x, 0, m_Dims.x-1);
//		v.y = CLAMP(v.y, 0, m_Dims.y-1);
//		v.z = CLAMP(v.z, 0, m_Dims.z-1);
//	}

private:
	void CalculateWorldBound();
	void CalculateVoxelDims(int voxel_density);
	void FillVoxels();
	void CalculateSDF();
	void Debug_SaveSDF();
	void CalculateSDF_Voxel(const Vec3i &ptCentre);
	void CalculateSDF_AssessNeighbour(const Vec3i &pt, unsigned int &min_SDF);
	bool VoxelWithinBounds(Vec3i v) const
	{
		if (v.x < 0) return false;
		if (v.y < 0) return false;
		if (v.z < 0) return false;
		if (v.x >= m_Dims.x) return false;
		if (v.y >= m_Dims.y) return false;
		if (v.z >= m_Dims.z) return false;
		return true;
	}
	void DebugCheckWorldPointInVoxel(Vector3 pt, const Vec3i &ptVoxel);
	void DebugCheckLocalPointInVoxel(Vector3 pt, const Vec3i &ptVoxel)
	{
//		assert ((int) (pt.x+0.5f) == ptVoxel.x);
//		assert ((int) (pt.y+0.5f) == ptVoxel.y);
//		assert ((int) (pt.z+0.5f) == ptVoxel.z);
	}

	LVector<Voxel> m_Voxels;
	LVector<AABB> m_VoxelBounds;
	LBitField_Dynamic m_BFTrisHit;
	int m_iNumTris;

	// slightly expanded
	AABB m_SceneWorldBound_expanded;
	AABB m_SceneWorldBound_mid;

	// exact
	AABB m_SceneWorldBound_contracted;

	const LightScene * m_pScene;

	Vec3i m_Dims;
	int m_DimsXTimesY;
	Vector3 m_VoxelSize;
	int m_iNumVoxels;
//	Plane m_VoxelPlanes[6];


	int GetVoxelNum(const Vec3i &pos) const
	{
		int v = pos.z * m_DimsXTimesY;
		v += pos.y * m_Dims.x;
		v += pos.x;
		return v;
	}

	const Voxel &GetVoxel(const Vec3i &pos) const
	{
		int v = GetVoxelNum(pos);
		assert (v < m_iNumVoxels);
		return m_Voxels[v];
	}
	Voxel &GetVoxel(const Vec3i &pos)
	{
		int v = GetVoxelNum(pos);
		assert (v < m_iNumVoxels);
		return m_Voxels[v];
	}

	// the pt is both the output intersection point, and the input plane constant (in the correct tuple)
	void IntersectAAPlane_OnlyWithinAABB(const AABB &aabb, const Ray &ray, int axis, Vector3 &pt, float &nearest_hit, int plane_id, int &nearest_plane_id, float aabb_epsilon = 0.0f) const
	{
		if (ray.IntersectAAPlane(axis, pt))
		{
			Vector3 offset = (pt - ray.o);
			float dist = offset.length_squared();
			if (dist < nearest_hit)
			{
				// must be within aabb
				if (aabb.has_point(pt))
				{
					nearest_hit = dist;
					nearest_plane_id = plane_id;
				}
			}
		}
	}


	// the pt is both the output intersection point, and the input plane constant (in the correct tuple)
	void IntersectAAPlane(const Ray &ray, int axis, Vector3 &pt, float &nearest_hit, int plane_id, int &nearest_plane_id) const
	{
		if (ray.IntersectAAPlane(axis, pt))
		{
			Vector3 offset = (pt - ray.o);
			float dist = offset.length_squared();
			if (dist < nearest_hit)
			{
				nearest_hit = dist;
				nearest_plane_id = plane_id;
			}
		}
	}

	bool IntersectRayAABB(const Ray &ray, const AABB &aabb, Vector3 &ptInter);
};

}
