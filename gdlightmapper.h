#pragma once

#include "llightmapper.h"
#include "scene/3d/mesh_instance.h"
#include "scene/3d/spatial.h"

class LLightmap : public Spatial {
	GDCLASS(LLightmap, Spatial);

public:
	enum eMode {
		MODE_FORWARD = LM::LightMapper::LMMODE_FORWARD,
		MODE_BACKWARD = LM::LightMapper::LMMODE_BACKWARD,
	};

	enum eBakeMode {
		BAKEMODE_UVMAP = LM::LightMapper::LMBAKEMODE_UVMAP,
		BAKEMODE_LIGHTMAP = LM::LightMapper::LMBAKEMODE_LIGHTMAP,
		BAKEMODE_AO = LM::LightMapper::LMBAKEMODE_AO,
		BAKEMODE_MERGE = LM::LightMapper::LMBAKEMODE_MERGE,
		BAKEMODE_PROBES = LM::LightMapper::LMBAKEMODE_PROBES,
		BAKEMODE_COMBINED = LM::LightMapper::LMBAKEMODE_COMBINED,
	};

	enum eQuality {
		QUALITY_LOW = LM::LightMapper::LM_QUALITY_LOW,
		QUALITY_MEDIUM = LM::LightMapper::LM_QUALITY_MEDIUM,
		QUALITY_HIGH = LM::LightMapper::LM_QUALITY_HIGH,
		QUALITY_FINAL = LM::LightMapper::LM_QUALITY_FINAL,
	};

	bool lightmap_mesh(Node *pMeshRoot, Node *pLightRoot, Object *pOutputImage_Lightmap, Object *pOutputImage_AO, Object *pOutputImage_Combined);
	bool lightmap_bake();
	bool lightmap_bake_to_image(Object *pOutputLightmapImage, Object *pOutputAOImage, Object *pOutputCombinedImage);

	bool uvmap();

	void set_mode(LLightmap::eMode p_mode);
	LLightmap::eMode get_mode() const;

	void set_bake_mode(LLightmap::eBakeMode p_mode);
	LLightmap::eBakeMode get_bake_mode() const;

	void set_quality(LLightmap::eQuality p_quality);
	LLightmap::eQuality get_quality() const;

	void set_max_light_distance(int dist);
	int get_max_light_distance() const;

	void set_mesh_path(const NodePath &p_path);
	NodePath get_mesh_path() const;
	void set_lights_path(const NodePath &p_path);
	NodePath get_lights_path() const;

	////////////////////////////
	void set_num_samples(int num_samples);
	int get_num_samples() const;

	//	void set_forward_ray_power(float ray_power);
	//	float get_forward_ray_power() const;

	void set_emission_density(float density);
	float get_emission_density() const;

	void set_glow(float glow);
	float get_glow() const;

	////////////////////////////
	void set_backward_num_rays(int num_rays);
	int get_backward_num_rays() const;

	//	void set_backward_num_bounces(int num_bounces);
	//	int get_backward_num_bounces() const;

	//	void set_backward_ray_power(float ray_power);
	//	float get_backward_ray_power() const;

	////////////////////////////
	void set_num_bounces(int num_bounces);
	int get_num_bounces() const;

	void set_bounce_power(float bounce_power);
	float get_bounce_power() const;

	void set_roughness(float roughness);
	float get_roughness() const;

	void set_num_ambient_bounces(int num_bounces);
	int get_num_ambient_bounces() const;

	void set_num_ambient_bounce_samples(int num_samples);
	int get_num_ambient_bounce_samples() const;

	void set_ambient_bounce_power(float bounce_power);
	float get_ambient_bounce_power() const;

	////////////////////////////

	void set_ao_range(float ao_range);
	float get_ao_range() const;

	void set_ao_cut_range(float ao_cut_range);
	float get_ao_cut_range() const;

	void set_ao_num_samples(int ao_num_samples);
	int get_ao_num_samples() const;
	////////////////////////////

	void set_tex_width(int width);
	int get_tex_width() const;

	void set_tex_height(int height);
	int get_tex_height() const;

	void set_voxel_density(int density);
	int get_voxel_density() const;

	void set_surface_bias(float bias);
	float get_surface_bias() const;

	void set_material_size(int size);
	int get_material_size() const;

	void set_normalize(bool norm);
	bool get_normalize() const;

	void set_normalize_multiplier(float bias);
	float get_normalize_multiplier() const;

	void set_light_ao_ratio(float ratio);
	float get_light_ao_ratio() const;

	void set_gamma(float gamma);
	float get_gamma() const;

	void set_lightmap_filename(const String &p_filename);
	String get_lightmap_filename() const;

	void set_ao_filename(const String &p_filename);
	String get_ao_filename() const;

	void set_combined_filename(const String &p_filename);
	String get_combined_filename() const;

	// UV
	void set_uv_filename(const String &p_filename);
	String get_uv_filename() const;

	void set_uv_padding(int pad);
	int get_uv_padding() const;

	// Probes
	void set_probe_density(int density);
	int get_probe_density() const;

	void set_probe_samples(int samples);
	int get_probe_samples() const;

	//	void set_probe_filename(const String &p_filename);
	//	String get_probe_filename() const;

	// Noise reduction
	void set_noise_reduction(float nr);
	float get_noise_reduction() const;

	void set_noise_threshold(float threshold);
	float get_noise_threshold() const;

	void set_noise_reduction_method(int method);
	int get_noise_reduction_method() const;

	void set_seam_stitching(bool active);
	bool get_seam_stitching() const;

	void set_seam_distance_threshold(float threshold);
	float get_seam_distance_threshold() const;

	void set_seam_normal_threshold(float threshold);
	float get_seam_normal_threshold() const;

	void set_visualize_seams(bool active);
	bool get_visualize_seams() const;

	void set_dilate(bool active);
	bool get_dilate() const;

private:
	LM::LightMapper m_LM;

	void ShowWarning(String sz, bool bAlert = true) { m_LM.ShowWarning(sz, bAlert); }

protected:
	static void _bind_methods();
};

VARIANT_ENUM_CAST(LLightmap::eMode);
VARIANT_ENUM_CAST(LLightmap::eBakeMode);
VARIANT_ENUM_CAST(LLightmap::eQuality);
