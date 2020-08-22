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
	void Create() {pAlbedo = 0; pGodotMaterial = 0; m_bEmitter = false; m_Power_Emission = 0.0f;}
	void Destroy();

	const Material * pGodotMaterial;
	LTexture * pAlbedo;

	bool m_bEmitter;

	float m_Power_Emission;
	Color m_Col_Emission; // color multiplied by emission power
};

class LMaterials
{
public:
	LMaterials();
	~LMaterials();
	void Reset();
	void Prepare(unsigned int max_material_size) {m_uiMaxMaterialSize = max_material_size;}

	int FindOrCreateMaterial(const MeshInstance &mi, Ref<Mesh> rmesh, int surf_id);
	bool FindColors(int mat_id, const Vector2 &uv, Color &albedo);

	const LMaterial &GetMaterial(int i) const {return m_Materials[i];}

	// in order to account for emission density we need to reduce
	// power to keep brightness the same
	void AdjustMaterials(float emission_density);

private:
	Variant FindCustom_AlbedoTex(Ref<Material> src_material);
	void FindCustom_ShaderParams(Ref<Material> src_material, float &emission, Color &emission_color);

	LTexture * _get_bake_texture(Ref<Image> p_image, const Color &p_color_mul, const Color &p_color_add);
	LTexture * _make_dummy_texture(LTexture * pLTexture, Color col);

	LVector<LMaterial> m_Materials;
	unsigned int m_uiMaxMaterialSize;
};



} // namespace
