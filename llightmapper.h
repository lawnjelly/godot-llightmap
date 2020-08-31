#pragma once

#include "lambient_occlusion.h"

namespace LM {

class LightMapper : public AmbientOcclusion
{
public:

	// main function called from the godot class
	bool lightmap_mesh(Spatial * pMeshesRoot, Spatial * pLR, Image * pIm_Lightmap, Image * pIm_AO, Image * pIm_Combined);
	bool uv_map_meshes(Spatial * pRoot);

private:
	bool LightmapMesh(Spatial * pMeshesRoot, const Spatial &light_root, Image &out_image_lightmap, Image &out_image_ao, Image &out_image_combined);
	void Reset();

private:
	// forward tracing
	void ProcessLights();
	void ProcessLight(int light_id, int num_rays);

	void ProcessEmissionTris();
	void ProcessEmissionTris_Section(float fraction_of_total);
	void ProcessEmissionTri(int etri_id, float fraction_of_total);

	// backward tracing
	void ProcessTexels();
	void ProcessTexel(int tx, int ty);
	void ProcessTexel_Light(int light_id, const Vector3 &ptSource, const Vector3 &ptNormal, FColor &color); //, uint32_t tri_ignore);

	void ProcessTexels_Bounce();
	FColor ProcessTexel_Bounce(int x, int y);
	bool ProcessTexel_Bounce_Sample(const Vector3 &plane_norm, const Vector3 &ray_origin, FColor &total_col);

	// light probes
	void ProcessLightProbes();

	void Refresh_Process_State();
//	const int m_iRaysPerSection = 1024 * 1024 * 4; // 64
	const int m_iRaysPerSection = 1024 * 64; // 64
	// 1024*1024 is 46 megs
};


} // namespace
