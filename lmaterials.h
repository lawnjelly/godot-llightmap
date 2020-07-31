#pragma once

#include "scene/3d/mesh_instance.h"
#include "lvector.h"

namespace LM
{

struct LTexture
{
	Vector<Color> colors;
	int width;
	int height;

	void Sample(const Vector2 &uv, Color &col) const;
};


struct LMaterial
{
	void Create() {pAlbedo = 0; pGodotMaterial = 0;}
	void Destroy();

	const Material * pGodotMaterial;
	LTexture * pAlbedo;
};

class LMaterials
{
public:
	LMaterials();
	~LMaterials();
	void Reset();

	int FindOrCreateMaterial(const MeshInstance &mi, Ref<Mesh> rmesh, int surf_id);
	bool FindColors(int mat_id, const Vector2 &uv, Color &albedo);


private:
	Variant FindShaderTex(Ref<Material> src_material);
	LTexture * _get_bake_texture(Ref<Image> p_image, const Color &p_color_mul, const Color &p_color_add);
	LTexture * _make_dummy_texture(LTexture * pLTexture, Color col);

	LVector<LMaterial> m_Materials;
};



} // namespace
