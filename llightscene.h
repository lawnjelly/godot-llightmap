#pragma once

#include "lvector.h"
#include "llighttypes.h"
#include "llighttracer.h"
#include "llightimage.h"
#include "lmaterials.h"
#include "scene/3d/mesh_instance.h"



namespace LM
{

class LightMapper_Base;

class LightScene
{
	friend class LLightMapper;
	friend class LightTracer;
public:

	// each texel can have more than 1 triangle crossing it.
	// this is important for sub texel anti-aliasing, so we need a quick record
	// of all the triangles to test
	struct TexelTriList
	{
		enum {MAX_TRIS = 11};
		uint32_t tri_ids[MAX_TRIS];
		uint32_t num_tris;
	};

	void Reset();
	bool Create(Spatial * pMeshesRoot, int width, int height, int voxel_density, int max_material_size, float emission_density);

	// returns triangle ID (or -1) and barycentric coords
	int FindIntersect_Ray(const Ray &ray, float &u, float &v, float &w, float &nearest_t, const Vec3i * pVoxelRange = nullptr);//, int ignore_tri_p1 = 0);
	int FindIntersect_Ray(const Ray &ray, Vector3 &bary, float &nearest_t, const Vec3i * pVoxelRange = nullptr)
	{
		return FindIntersect_Ray(ray, bary.x, bary.y, bary.z, nearest_t, pVoxelRange);
	}

	// simple test returns true if any collision
	bool TestIntersect_Ray(const Ray &ray, float max_dist, const Vec3i &voxel_range, bool bCullBackFaces = false);
	bool TestIntersect_Ray(const Ray &ray, float max_dist, bool bCullBackFaces = false);
	bool TestIntersect_Line(const Vector3 &a, const Vector3 &b, bool bCullBackFaces = false);

	void FindUVsBarycentric(int tri, Vector2 &uvs, float u, float v, float w) const
	{
		m_UVTris[tri].FindUVBarycentric(uvs, u, v, w);
	}
	void FindUVsBarycentric(int tri, Vector2 &uvs, const Vector3 &bary) const
	{
		m_UVTris[tri].FindUVBarycentric(uvs, bary);
	}

	bool FindPrimaryTextureColors(int tri_id, const Vector3 &bary, Color &albedo, bool &bTransparent);
	bool FindEmissionColor(int tri_id, const Vector3 &bary, Color &texture_col, Color &col);

	// single function
	bool FindAllTextureColors(int tri_id, const Vector3 &bary, Color &albedo, Color &emission, bool &bTransparent, bool &bEmitter);


	// setup
	void RasterizeTriangleIDs(LightMapper_Base &base, LightImage<uint32_t> &im_p1, LightImage<Vector3> &im_bary);
	int GetNumTris() const {return m_UVTris.size();}

private:
	void FindMeshes(Spatial * pNode);

	bool Create_FromMesh(int mesh_id, int width, int height);
	bool Create_FromMeshSurface(int mesh_id, int surf_id, Ref<Mesh> rmesh, int width, int height);


	void CalculateTriTexelSize(int tri_id, int width, int height);

	void Transform_Verts(const PoolVector<Vector3> &ptsLocal, PoolVector<Vector3> &ptsWorld, const Transform &tr) const;
	void Transform_Norms(const PoolVector<Vector3> &normsLocal, PoolVector<Vector3> &normsWorld, const Transform &tr) const;
	void ProcessVoxelHits(const Ray &ray, const PackedRay &pray, const Voxel &voxel, float &r_nearest_t, int &r_nearest_tri); // int ignore_triangle_id_p1);
	void ProcessVoxelHits_Old(const Ray &ray, const Voxel &voxel, float &r_nearest_t, int &r_nearest_tri);

	bool TestVoxelHits(const Ray &ray, const PackedRay &pray, const Voxel &voxel, float max_dist, bool bCullBackFaces);


//	PoolVector<Vector3> m_ptPositions;
//	PoolVector<Vector3> m_ptNormals;
	//PoolVector<Vector2> m_UVs;
	//PoolVector<int> m_Inds;

	LVector<Rect2> m_TriUVaabbs;
	LVector<AABB> m_TriPos_aabbs;
	LVector<MeshInstance *> m_Meshes;

	// precalculate these .. useful for finding cutting tris.
	LVector<float> m_Tri_TexelSizeWorldSpace;

	LMaterials m_Materials;

protected:
public:
	LVector<UVTri> m_UVTris;

	LightTracer m_Tracer;

	LVector<Tri> m_Tris;
	LVector<Tri> m_TriNormals;

	LVector<Plane> m_TriPlanes;

	// we maintain a list of tris in the form of 2 edges plus a point. These
	// are precalculated as they are used in the intersection test.
	LVector<Tri> m_Tris_EdgeForm;

	// stuff for mapping to original meshes
//	LVector<int> m_Tri_MeshIDs;
//	LVector<uint16_t> m_Tri_SurfIDs;

	// these are plus 1
	LVector<uint16_t> m_Tri_LMaterialIDs;

	// a list of triangles that have emission materials
	LVector<EmissionTri> m_EmissionTris;


	// these are UVs in the first channel, if present, or 0.
	// This allows mapping back to the albedo etc texture.
	LVector<UVTri> m_UVTris_Primary;

	Vec3i m_VoxelRange;

	bool m_bUseSIMD;
};



}
