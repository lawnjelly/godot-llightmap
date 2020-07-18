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
	void ProcessLights();
	void ProcessLight(int light_id, int num_rays);
	void ProcessRay(LM::Ray r, int depth, float power, int dest_tri_id = 0, const Vector2i * pUV = 0);


	void ProcessTexels();
	void ProcessTexel(int tx, int ty);
	float ProcessTexel_Light(int light_id, const Vector3 &ptDest, const Vector3 &ptNormal, uint32_t tri_ignore);

	void ProcessTexels_Bounce();
	float ProcessTexel_Bounce(int x, int y);

	void Refresh_Process_State();

	const int m_iRaysPerSection = 1024 * 1024 * 4; // 64
	// 1024*1024 is 46 megs
};


} // namespace
