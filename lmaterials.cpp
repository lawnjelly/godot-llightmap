#include "lmaterials.h"

#define LM_STRING_TRANSPARENT "_T_"

namespace LM
{

void LMaterial::Destroy()
{
	if (pAlbedo)
	{
		memdelete(pAlbedo);
		pAlbedo = 0;
	}
}

/////////////////////////////////

LMaterials::LMaterials()
{
	m_uiMaxMaterialSize = 256;
}

LMaterials::~LMaterials()
{
	Reset();
}

void LMaterials::Reset()
{
	for (int n=0; n<m_Materials.size(); n++)
	{
		m_Materials[n].Destroy();
	}

	m_Materials.clear(true);
}

void LMaterials::AdjustMaterials(float emission_density)
{
	if (emission_density == 0.0f)
		return;

	float emission_multiplier = 1.0f / emission_density;

	for (int n=0; n<m_Materials.size(); n++)
	{
		LMaterial &mat = m_Materials[n];

		if (mat.m_bEmitter)
		{
			mat.m_Power_Emission *= emission_multiplier;
			mat.m_Col_Emission *= emission_multiplier;
		}
	}

}

int LMaterials::FindOrCreateMaterial(const MeshInstance &mi, Ref<Mesh> rmesh, int surf_id)
{
	Ref<Material> src_material;

	// mesh instance has the material?
	src_material = mi.get_surface_material(surf_id);
	if (src_material.ptr())
	{
		//mi.set_surface_material(0, mat);
	}
	else
	{
		// mesh has the material?
		src_material = rmesh->surface_get_material(surf_id);
		//mi.set_surface_material(0, smat);
	}



//	Ref<Material> src_material = rmesh->surface_get_material(surf_id);
	const Material * pSrcMaterial = src_material.ptr();

	if (!pSrcMaterial)
		return 0;

	// already exists?
	for (int n=0; n<m_Materials.size(); n++)
	{
		if (m_Materials[n].pGodotMaterial == pSrcMaterial)
			return n+1;
	}

	// doesn't exist create a new material
	LMaterial * pMat = m_Materials.request();
	pMat->Create();
	pMat->pGodotMaterial = pSrcMaterial;

	// spatial material?
	Ref<SpatialMaterial> spatial_mat = src_material;
	Ref<Texture> albedo_tex;
	Color albedo = Color(1, 1, 1, 1);


	float emission = 0.0f;
	Color emission_color(1, 1, 1, 1);

	if (spatial_mat.is_valid())
	{
		albedo_tex = spatial_mat->get_texture(SpatialMaterial::TEXTURE_ALBEDO);
		albedo = spatial_mat->get_albedo();

		emission = spatial_mat->get_emission_energy();
		emission_color = spatial_mat->get_emission();
	}
	else
	{
		// shader material?
		Variant shader_tex = FindCustom_AlbedoTex(src_material);
		albedo_tex = shader_tex;

		// check the name of the material to allow emission
		String szMat = src_material->get_path();
		if (szMat.find(LM_STRING_TRANSPARENT) != -1)
		{
			pMat->m_bTransparent = true;
		}

		FindCustom_ShaderParams(src_material, emission, emission_color);

	} // not spatial mat

	// emission
	if (emission > 0.0f)
	{
		pMat->m_bEmitter = true;
		pMat->m_Power_Emission = emission / 1000.0f; // some constant to be comparable to lights
		pMat->m_Col_Emission = emission_color * pMat->m_Power_Emission;

		// apply a modifier for the emission density. As the number of samples go up, the power per sample
		// must reduce in order to prevent brightness changing.
	}

	Ref<Image> img_albedo;
	if (albedo_tex.is_valid())
	{
		img_albedo = albedo_tex->get_data();
		pMat->pAlbedo = _get_bake_texture(img_albedo, albedo, Color(0, 0, 0)); // albedo texture, color is multiplicative
		//albedo_texture = _get_bake_texture(img_albedo, size, mat->get_albedo(), Color(0, 0, 0)); // albedo texture, color is multiplicative
	} else
	{
		//albedo_texture = _get_bake_texture(img_albedo, size, Color(1, 1, 1), mat->get_albedo()); // no albedo texture, color is additive
	}

	// emission?


	// returns the new material ID plus 1
	return m_Materials.size();
}

void LMaterials::FindCustom_ShaderParams(Ref<Material> src_material, float &emission, Color &emission_color)
{
	// defaults
	emission = 0.0f;
	emission_color = Color(1, 1, 1, 1);

	Ref<ShaderMaterial> shader_mat = src_material;

	if (!shader_mat.is_valid())
		return;

	Variant p_emission = shader_mat->get_shader_param("emission");
	if (p_emission)
		emission = p_emission;

	Variant p_emission_color = shader_mat->get_shader_param("emission_color");
	if (p_emission_color)
		emission_color = p_emission_color;

}


Variant LMaterials::FindCustom_AlbedoTex(Ref<Material> src_material)
{
	Ref<ShaderMaterial> shader_mat = src_material;

	if (!shader_mat.is_valid())
		return Variant::NIL;

	// get the shader
	Ref<Shader> shader = shader_mat->get_shader();
	if (!shader.is_valid())
		return Variant::NIL;

	// first - is there a named albedo texture?
	Variant named_param = shader_mat->get_shader_param("texture_albedo");
//	if (named_param)
//		return named_param;
	return named_param;

	/*
	// find the most likely albedo texture
	List<PropertyInfo> plist;
	shader->get_param_list(&plist);

	String sz_first_obj_param;

	for (List<PropertyInfo>::Element *E = plist.front(); E; E = E->next()) {
		String szName = E->get().name;
		Variant::Type t = E->get().type;
//				print_line("shader param : " + szName);
//				print_line("shader type : " + String(Variant(t)));
		if (t == Variant::OBJECT)
		{
			sz_first_obj_param = szName;
			break;
		}

		//r_options->push_back(quote_style + E->get().name.replace_first("shader_param/", "") + quote_style);
	}

	if (sz_first_obj_param == "")
		return Variant::NIL;

	StringName pr = shader->remap_param(sz_first_obj_param);
	if (!pr) {
		String n = sz_first_obj_param;
		if (n.find("param/") == 0) { //backwards compatibility
			pr = n.substr(6, n.length());
		}
		if (n.find("shader_param/") == 0) { //backwards compatibility
			pr = n.replace_first("shader_param/", "");
		}
	}

	if (!pr)
		return Variant::NIL;

	Variant param = shader_mat->get_shader_param(pr);

	print_line("\tparam is " + String(param));
	return param;
	*/
}


LTexture * LMaterials::_make_dummy_texture(LTexture * pLTexture, Color col)
{
	pLTexture->colors.resize(1);
	pLTexture->colors.set(0, col);
	pLTexture->width = 1;
	pLTexture->height = 1;
	return pLTexture;
}

LTexture * LMaterials::_get_bake_texture(Ref<Image> p_image, const Color &p_color_mul, const Color &p_color_add)
{
	LTexture * lt = memnew(LTexture);

	// no image exists, use dummy texture
	if (p_image.is_null() || p_image->empty())
	{
		return _make_dummy_texture(lt, p_color_mul);
	}

	p_image = p_image->duplicate();

	if (p_image->is_compressed()) {
		Error err = p_image->decompress();
		if (err != OK)
		{
			// could not decompress
			WARN_PRINT("LMaterials::_get_bake_texture : could not decompress texture");
			return _make_dummy_texture(lt, p_color_mul);
		}
	}

	p_image->convert(Image::FORMAT_RGBA8);

	// downsize if necessary
	int w = p_image->get_width();
	int h = p_image->get_height();

	bool bResize = false;
	while (true)
	{
		if ((w > m_uiMaxMaterialSize) || (h > m_uiMaxMaterialSize))
		{
			w /= 2;
			h /= 2;
			bResize = true;
		}
		else
		{
			w = MAX(w, 1);
			h = MAX(h, 1);
			break;
		}
	}
	if (bResize)
		p_image->resize(w, h, Image::INTERPOLATE_CUBIC);


	w = p_image->get_width();
	h = p_image->get_height();
	int size = w * h;

	lt->width = w;
	lt->height = h;

	lt->colors.resize(size);

	PoolVector<uint8_t>::Read r = p_image->get_data().read();

	for (int i = 0; i<size; i++)
	{
		Color c;
		c.r = (r[i * 4 + 0] / 255.0) * p_color_mul.r + p_color_add.r;
		c.g = (r[i * 4 + 1] / 255.0) * p_color_mul.g + p_color_add.g;
		c.b = (r[i * 4 + 2] / 255.0) * p_color_mul.b + p_color_add.b;

		// srgb to linear?


		c.a = r[i * 4 + 3] / 255.0;

		lt->colors.set(i, c);
	}

	return lt;
}

void LTexture::Sample(const Vector2 &uv, Color &col) const
{
	// mod to surface (tiling)
	float x = fmodf(uv.x, 1.0f);
	float y = fmodf(uv.y, 1.0f);

	// we need these because fmod can produce negative results
	if (x < 0.0f)
		x = 1.0f + x;
	if (y < 0.0f)
		y = 1.0f + y;

	x *= width;
	y *= height;

	// no filtering as yet
	int tx = x;
	int ty = y;

	tx = MIN(tx, width);
	ty = MIN(ty, height);

	int i = (ty * width) + tx;

	col = colors[i];
}


bool LMaterials::FindColors(int mat_id, const Vector2 &uv, Color &albedo, bool &bTransparent)
{
	// mat_id is plus one
	if (!mat_id)
	{
		albedo = Color(1, 1, 1, 1);
		bTransparent = false;
		return false;
	}

	mat_id--;
	const LMaterial &mat = m_Materials[mat_id];

	// return whether transparent
	bTransparent = mat.m_bTransparent;

	if (!mat.pAlbedo)
	{
		albedo = Color(1, 1, 1, 1);
		return false;
	}

	const LTexture &tex = *mat.pAlbedo;
	tex.Sample(uv, albedo);
	return true;
}


} // namespace
