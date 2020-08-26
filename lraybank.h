#pragma once

#include "llighttypes.h"
#include "lvector.h"
#include "llightmapper_base.h"

#define RAYBANK_USE_THREADING

namespace LM
{

struct RB_Voxel
{
	LVector<FRay> m_Rays;
};

class RayBank : public LightMapper_Base
{
public:
	void RayBank_Reset();
	void RayBank_Create();

	// every time we want to queue a new ray for processing
	FRay * RayBank_RequestNewRay(Ray ray,  int num_rays_left, const FColor &col, const Vec3i * pStartVoxel = nullptr);

	// multithread accelerated .. do intersection tests on rays, calculate hit points and new rays
	void RayBank_Process();

	// flush ray results to the lightmap
	void RayBank_Flush();

	void RayBank_CheckVoxelsClear();
private:
	// used for below multithread routine
	RB_Voxel * m_pCurrentThreadVoxel;
	void RayBank_ProcessRay_MT(uint32_t ray_id, int start_ray);
	void RayBank_ProcessRay_MT_Old(uint32_t ray_id, int start_ray);

	void RayBank_FlushRay(RB_Voxel &vox, int ray_id);

	RB_Voxel &RayBank_GetVoxelWrite(const Vec3i &pt) {int n = GetTracer().GetVoxelNum(pt); return m_Data_RB.GetVoxels_Write()[n];}
	RB_Voxel &RayBank_GetVoxelRead(const Vec3i &pt) {int n = GetTracer().GetVoxelNum(pt); return m_Data_RB.GetVoxels_Read()[n];}

public:
	LightTracer &GetTracer() {return m_Scene.m_Tracer;}
	const LightTracer &GetTracer() const {return m_Scene.m_Tracer;}
private:

	struct RayBank_Data
	{
		LVector<RB_Voxel> &GetVoxels_Read() {return m_Voxels[m_MapRead];}
		LVector<RB_Voxel> &GetVoxels_Write() {return m_Voxels[m_MapWrite];}
		LVector<RB_Voxel> m_Voxels[2];
		void Swap();
		int m_MapRead;
		int m_MapWrite;
	} m_Data_RB;

protected:
	bool HitBackFace(const Ray &r, int tri_id, const Vector3 &bary, Vector3 &face_normal) const
	{
		const Tri &triangle_normal = m_Scene.m_TriNormals[tri_id];
		triangle_normal.InterpolateBarycentric(face_normal, bary);
		face_normal.normalize(); // is this necessary as we are just checking a dot product polarity?

		float dot = face_normal.dot(r.d);
		if (dot >= 0.0f)
		{
			return true;
		}

		return false;
	}

	void CalculateTransmittance(const Color &albedo, FColor &ray_color)
	{
		// rapidly converge to the surface color
		float surf_fraction = albedo.a * 2.0f;
		if (surf_fraction > 1.0f)
			surf_fraction = 1.0f;

		FColor mixed_color;
		mixed_color.Set(albedo);
		mixed_color = mixed_color * ray_color;

		ray_color.Lerp(mixed_color, surf_fraction);

		// darken
		float dark_fraction = albedo.a;
		dark_fraction -= 0.5f;
		dark_fraction *= 2.0f;

		if (dark_fraction < 0.0f)
			dark_fraction = 0.0f;

		ray_color *= 1.0f - dark_fraction;
	}

};


} // namespace
