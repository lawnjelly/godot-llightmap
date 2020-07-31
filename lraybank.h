#pragma once

#include "llighttypes.h"
#include "lvector.h"
#include "llightmapper_base.h"


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
	FRay * RayBank_RequestNewRay(Ray ray,  int num_rays_left, const FColor &col, const Vec3i * pStartVoxel);

	// can be used from several threads
	void RayBank_Process();

	// flush ray results to the lightmap
	void RayBank_Flush();

	void RayBank_CheckVoxelsClear();
private:
	// used for below multithread routine
	RB_Voxel * m_pCurrentThreadVoxel;
	void RayBank_ProcessRay_MT(uint32_t ray_id, int start_ray);
//	void RayBank_ProcessRay(uint32_t ray_id, RB_Voxel &vox);

	void RayBank_FlushRay(RB_Voxel &vox, int ray_id);

	RB_Voxel &RayBank_GetVoxelWrite(const Vec3i &pt) {int n = GetTracer().GetVoxelNum(pt); return m_Data_RB.GetVoxels_Write()[n];}
	RB_Voxel &RayBank_GetVoxelRead(const Vec3i &pt) {int n = GetTracer().GetVoxelNum(pt); return m_Data_RB.GetVoxels_Read()[n];}

	LightTracer &GetTracer() {return m_Scene.m_Tracer;}
	const LightTracer &GetTracer() const {return m_Scene.m_Tracer;}

	struct RayBank_Data
	{
		LVector<RB_Voxel> &GetVoxels_Read() {return m_Voxels[m_MapRead];}
		LVector<RB_Voxel> &GetVoxels_Write() {return m_Voxels[m_MapWrite];}
		LVector<RB_Voxel> m_Voxels[2];
		void Swap();
		int m_MapRead;
		int m_MapWrite;
	} m_Data_RB;

};


} // namespace
