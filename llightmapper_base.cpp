#include "llightmapper.h"
#include "ldilate.h"
#include "core/os/os.h"
#include "scene/3d/light.h"
#include "core/math/plane.h"

using namespace LM;

LightMapper_Base::BakeBeginFunc LightMapper_Base::bake_begin_function = NULL;
LightMapper_Base::BakeStepFunc LightMapper_Base::bake_step_function = NULL;
LightMapper_Base::BakeEndFunc LightMapper_Base::bake_end_function = NULL;

LightMapper_Base::LightMapper_Base()
{
	m_iNumRays = 1;
	m_Settings_Forward_NumRays = 16;
	m_Settings_Forward_NumBounces = 0;
	m_Settings_Forward_RayPower = 0.01f;
	m_Settings_Forward_BouncePower = 1.0f;
	m_Settings_Forward_BounceDirectionality = 0.5f;
	m_Settings_Forward_Emission_Density = 1.0f;

	m_Settings_Backward_NumRays = 128;
	m_Settings_Backward_NumBounceRays = 128;
	m_Settings_Backward_NumBounces = 0;
	m_Settings_Backward_RayPower = 1.0f;
	m_Settings_Backward_BouncePower = 0.5f;

	m_Settings_AO_Range = 2.0f;
	m_Settings_AO_Samples = 256;
	m_Settings_AO_CutRange = 1.5f;
	m_Settings_AO_ReverseBias = 0.005f;

	m_Settings_Mode = LMMODE_FORWARD;
	m_Settings_BakeMode = LMBAKEMODE_LIGHTMAP;
	m_Settings_Quality = LM_QUALITY_MEDIUM;

	m_Settings_TexWidth = 512;
	m_Settings_TexHeight = 512;
	m_Settings_VoxelDensity = 20;
	m_Settings_SurfaceBias = 0.005f;

	m_Settings_Max_Material_Size = 256;

	m_Settings_Normalize = true;
	m_Settings_NormalizeBias = 4.0f;
	m_Settings_Light_AO_Ratio = 0.5f;
	m_Settings_Gamma = 2.2f;

	m_Settings_LightmapIsHDR = false;
	m_Settings_AmbientIsHDR = false;
	m_Settings_CombinedIsHDR = false;

	m_Settings_Process_Lightmap = true;
	m_Settings_Process_AO = true;

	m_Settings_UVPadding = 4;
}

void LightMapper_Base::CalculateQualityAdjustedSettings()
{
	// set them initially to the same
	AdjustedSettings &as = m_AdjustedSettings;

	as.m_Forward_NumRays = m_Settings_Forward_NumRays;
	as.m_Forward_NumBounces = m_Settings_Forward_NumBounces;
	as.m_Forward_Emission_Density = m_Settings_Forward_Emission_Density;

	as.m_Backward_NumRays= m_Settings_Backward_NumRays;
	as.m_Backward_NumBounceRays = m_Settings_Backward_NumBounceRays;
	as.m_Backward_NumBounces = m_Settings_Backward_NumBounces;

	as.m_AO_Samples = m_Settings_AO_Samples;

	as.m_Max_Material_Size = m_Settings_Max_Material_Size;

	// overrides
	switch (m_Settings_Quality)
	{
	case LM_QUALITY_LOW:
		{
			as.m_Forward_NumRays = 1;
			as.m_Forward_NumBounces = 0;
			as.m_Backward_NumRays = 4;
			as.m_Backward_NumBounces = 0;
			as.m_AO_Samples = 1;
			as.m_Max_Material_Size = 32;
		}
		break;
	case LM_QUALITY_MEDIUM:
		{
			as.m_Forward_NumRays /= 2;
			as.m_Backward_NumRays /= 2;
			as.m_Backward_NumBounceRays /= 2;
			as.m_AO_Samples /= 2;
			as.m_Max_Material_Size /= 4;
		}
		break;
	default:
		// high is default
		break;
	case LM_QUALITY_FINAL:
		as.m_Forward_NumRays *= 2;
		as.m_Backward_NumRays *= 2;
		as.m_Backward_NumBounceRays *= 2;
		as.m_AO_Samples *= 2;
		break;
	}

	// minimums
	as.m_Forward_NumRays = MAX(as.m_Forward_NumRays, 1);
	as.m_Backward_NumRays = MAX(as.m_Backward_NumRays, 1);
	as.m_Backward_NumBounceRays = MAX(as.m_Backward_NumBounceRays, 1);
	as.m_AO_Samples = MAX(as.m_AO_Samples, 1);
	as.m_Max_Material_Size = MAX(as.m_Max_Material_Size, 32);

}

void LightMapper_Base::FindLight(const Node * pNode)
{
	const Light * pLight = Object::cast_to<const Light>(pNode);
	if (!pLight)
		return;

	// visibility or bake mode?
	if (pLight->get_bake_mode() == Light::BAKE_DISABLED)
		return;

	// is it visible?
//	if (!pLight->is_visible_in_tree())
//		return;

	LLight * l = m_Lights.request();

	// blank
	memset (l, 0, sizeof (LLight));

	l->m_pLight = pLight;
	// get global transform only works if glight is in the tree
	Transform trans = pLight->get_global_transform();
	l->pos = trans.origin;
	l->dir = -trans.basis.get_axis(2); // or possibly get_axis .. z is what we want
	l->dir.normalize();

	trans = pLight->get_transform();
	l->scale = trans.basis.get_scale();

	l->energy = pLight->get_param(Light::PARAM_ENERGY);
	l->indirect_energy = pLight->get_param(Light::PARAM_INDIRECT_ENERGY);
	l->range = pLight->get_param(Light::PARAM_RANGE);
	l->spot_angle_radians = Math::deg2rad(pLight->get_param(Light::PARAM_SPOT_ANGLE));
	l->spot_dot_max = Math::cos(l->spot_angle_radians);
	l->spot_dot_max = MIN(l->spot_dot_max, 0.9999f); // just to prevent divide by zero in cone of spotlight

	// the spot emanation point is used for spotlight cone culling.
	// if we used the dot from the pos to cull, we would miss cases where
	// the sample origin is offset by scale from pos. So we push back the pos
	// in order to account for the scale 'cloud' of origins.
	float radius = MAX(l->scale.x, MAX(l->scale.y, l->scale.z));
	l->spot_emanation_point =l->pos - (l->dir * radius);

	// pre apply intensity
	l->color.Set(pLight->get_color() * l->energy);

	const DirectionalLight * pDLight = Object::cast_to<DirectionalLight>(pLight);
	if (pDLight)
		l->type = LLight::LT_DIRECTIONAL;

	const SpotLight * pSLight = Object::cast_to<SpotLight>(pLight);
	if (pSLight)
		l->type = LLight::LT_SPOT;

	const OmniLight * pOLight = Object::cast_to<OmniLight>(pLight);
	if (pOLight)
		l->type = LLight::LT_OMNI;

}

void LightMapper_Base::PrepareLights()
{
	for (int n=0; n<m_Lights.size(); n++)
	{
		LLight &light = m_Lights[n];

		if (light.type == LLight::LT_DIRECTIONAL)
			LightToPlane(light);
	}
}


void LightMapper_Base::FindLights_Recursive(const Node * pNode)
{
	FindLight(pNode);

	int nChildren = pNode->get_child_count();

	for (int n=0; n<nChildren; n++)
	{
		Node * pChild = pNode->get_child(n);
		FindLights_Recursive(pChild);
	}
}

Plane LightMapper_Base::FindContainmentPlane(const Vector3 &dir, Vector3 pts[8], float &range, float padding)
{
	// construct a plane with one of the points
	Plane pl(pts[0], dir);

	float furthest_dist = 0.0f;
	int furthest = 0;

	// find the furthest point
	for (int n=0; n<8; n++)
	{
		float d = pl.distance_to(pts[n]);

		if (d < furthest_dist)
		{
			furthest_dist = d;
			furthest = n;
		}
	}

	// reconstruct the plane based on the furthest point
	pl = Plane(pts[furthest], dir);

	// find the range
	range = 0.0f;

	for (int n=0; n<8; n++)
	{
		float d = pl.distance_to(pts[n]);

		if (d > range)
		{
			range = d;
		}
	}

//	const float padding = 8.0f;

	// move plane backward a bit for luck
	pl.d -= padding;

	// add a boost to the range
	range += padding * 2.0f;

	return pl;
}

void LightMapper_Base::LightToPlane(LLight &light)
{
	AABB bb = m_Scene.m_Tracer.GetWorldBound_expanded();
	Vector3 minmax[2];
	minmax[0] = bb.position;
	minmax[1] = bb.position + bb.size;

	if (light.dir.y == 0.0f)
		return;

	// find the shift in x and z caused by y offset to top of scene
	float units = bb.size.y / light.dir.y;
	Vector3 offset = light.dir * -units;

	// add to which side of scene
	if (offset.x >= 0.0f)
		minmax[1].x += offset.x;
	else
		minmax[0].x += offset.x;

	if (offset.z >= 0.0f)
		minmax[1].z += offset.z;
	else
		minmax[0].z += offset.z;


	light.dl_plane_pt = minmax[0];
	light.dl_tangent = Vector3(1, 0, 0);
	light.dl_bitangent = Vector3(0, 0, 1);
	light.dl_tangent_range = minmax[1].x - minmax[0].x;
	light.dl_bitangent_range = minmax[1].z - minmax[0].z;

	print_line("plane mins : " + String(Variant(minmax[0])));
	print_line("plane maxs : " + String(Variant(minmax[1])));


	return;




	Vector3 pts[8];
	for (int n=0; n<8; n++)
	{
		// which x etc. either 0 or 1 for each axis
		int wx = MIN(n & 1, 1);
		int wy = MIN(n & 2, 1);
		int wz = MIN(n & 4, 1);

		pts[n].x = minmax[wx].x;
		pts[n].y = minmax[wy].y;
		pts[n].z = minmax[wz].z;
	}

	// new .. don't use light direction as plane normal, always use
	// ceiling type plane (for sky) or from below.
	// This will deal with most common cases .. for side lights,
	// area light is better.
	Vector3 plane_normal = light.dir;
	if (light.dir.y < 0.0f)
	{
		plane_normal = Vector3(0, -1, 0);
	}
	else
	{
		plane_normal = Vector3(0, 1, 0);
	}

	const float PLANE_PUSH = 2.0f;

	float main_range;
	Plane pl = FindContainmentPlane(plane_normal, pts, main_range, PLANE_PUSH);

	// push it back even further for safety
	//pl.d -= 2.0f;

	// now create a bound on this plane

	// find a good tangent
	Vector3 cross[3];
	cross[0] = plane_normal.cross(Vector3(1, 0, 0));
	cross[1] = plane_normal.cross(Vector3(0, 1, 0));
	cross[2] = plane_normal.cross(Vector3(0, 0, 1));

	float lx = cross[0].length();
	float ly = cross[1].length();
	float lz = cross[2].length();

	int best_cross = 0;
	if (ly > lx)
	{
		best_cross = 1;
		if (lz > ly)
			best_cross = 2;
	}
	if (lz > lx)
		best_cross = 2;

	Vector3 tangent = cross[best_cross];
	tangent.normalize();

	Vector3 bitangent = plane_normal.cross(tangent);
	bitangent.normalize();

	float tangent_range;
	Plane pl_tangent = FindContainmentPlane(tangent, pts, tangent_range, 0.0f);

	float bitangent_range;
	Plane pl_bitangent = FindContainmentPlane(bitangent, pts, bitangent_range, 0.0f);

	// find point at mins of the planes
	Vector3 ptPlaneMins;
	bool res = pl.intersect_3(pl_tangent, pl_bitangent, &ptPlaneMins);
	assert (res);

	// for flat sky, adjust the point to account for the incoming light direction
	// so as not to have part of the mesh in shadow
//	Vector3 offset = light.dir * -PLANE_PUSH;
//	ptPlaneMins.x += offset.x;
//	ptPlaneMins.z += offset.z;

	// we now have a point, 2 vectors (tangent and bitangent) and ranges,
	// all that is needed for a random distribution!

//	Vector3 dl_plane_pt;
//	Vector3 dl_plane_tangent;
//	Vector3 dl_plane_bitangent;
//	float dl_plane_tangent_range;
//	float dl_plane_bitangent_range;
	light.dl_plane_pt = ptPlaneMins;
	light.dl_tangent = tangent;
	light.dl_bitangent = bitangent;
	light.dl_tangent_range = tangent_range;
	light.dl_bitangent_range = bitangent_range;


	// debug output the positions
	Vector3 pA = ptPlaneMins;
	Vector3 pB = ptPlaneMins + (tangent * tangent_range);
	Vector3 pC = ptPlaneMins + (tangent * tangent_range) + (bitangent * bitangent_range);
	Vector3 pD = ptPlaneMins + (bitangent * bitangent_range);

	print_line("dir light A : " + String(Variant(pA)));
	print_line("dir light B : " + String(Variant(pB)));
	print_line("dir light C : " + String(Variant(pC)));
	print_line("dir light D : " + String(Variant(pD)));

}


void LightMapper_Base::PrepareImageMaps()
{
	m_Image_ID_p1.Blank();
	m_Image_ID2_p1.Blank();

	// rasterize each triangle in turn
	m_Scene.RasterizeTriangleIDs(*this, m_Image_ID_p1, m_Image_ID2_p1, m_Image_Barycentric);

	/*
	// go through each texel
	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			// use the texel centre!
			// find the triangle at this UV
			float u = (x + 0.5f) / (float) m_iWidth;
			float v = (y + 0.5f) / (float) m_iHeight;


			Vector3 &bary = *m_Image_Barycentric.Get(x, y);
			*m_Image_ID_p1.Get(x, y) = m_Scene.FindTriAtUV(u, v, bary.x, bary.y, bary.z);
		}
	}
	*/
}

void LightMapper_Base::Normalize_AO()
{
	int nPixels = m_Image_AO.GetNumPixels();
	float fmax = 0.0f;

	// first find the max
	for (int n=0; n<nPixels; n++)
	{
		float f = *m_Image_AO.Get(n);
		if (f > fmax)
			fmax = f;
	}

	if (fmax < 0.001f)
	{
		WARN_PRINT_ONCE("LightMapper_Base::Normalize_AO : values too small to normalize");
		return;
	}

	// multiplier to normal is 1.0f / fmax
	float mult = 1.0f / fmax;

	// apply bias
	//mult *= m_Settings_NormalizeBias;

	// apply multiplier
	for (int n=0; n<nPixels; n++)
	{
		float &f = *m_Image_AO.Get(n);
		f *= mult;

		// negate AO
		f = 1.0f - f;
		if (f < 0.0f)
			f = 0.0f;
	}

}

void LightMapper_Base::Normalize()
{
	if (!m_Settings_Normalize)
		return;

	int nPixels = m_Image_L.GetNumPixels();
	float fmax = 0.0f;

	// first find the max
	for (int n=0; n<nPixels; n++)
	{
		float f = m_Image_L.Get(n)->Max();
		if (f > fmax)
			fmax = f;
	}

	if (fmax < 0.001f)
	{
		WARN_PRINT_ONCE("LightMapper_Base::Normalize : values too small to normalize");
		return;
	}

	// multiplier to normal is 1.0f / fmax
	float mult = 1.0f / fmax;

	// apply bias
	mult *= m_Settings_NormalizeBias;

	// apply multiplier
	for (int n=0; n<nPixels; n++)
	{
		FColor &col = *m_Image_L.Get(n);
		col = col * mult;
	}
}


void LightMapper_Base::LoadLightmap(Image &image)
{
	assert (image.get_width() == m_iWidth);
	assert (image.get_height() == m_iHeight);

	Error res = image.load(m_Settings_LightmapFilename);
	if (res != OK)
	{
		WARN_PRINT_ONCE("LoadLightmap failed")
		return;
	}

	image.lock();
	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			m_Image_L.GetItem(x, y).Set(image.get_pixel(x, y));
		}
	}
	image.unlock();
}

void LightMapper_Base::LoadAO(Image &image)
{
	assert (image.get_width() == m_iWidth);
	assert (image.get_height() == m_iHeight);

	Error res = image.load(m_Settings_AmbientFilename);
	if (res != OK)
	{
		WARN_PRINT_ONCE("LoadAO failed")
		return;
	}

	image.lock();
	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			m_Image_AO.GetItem(x, y) = image.get_pixel(x, y).r;
		}
	}
	image.unlock();
}


void LightMapper_Base::Merge_AndWriteOutputImage_Combined(Image &image)
{
	// normalize lightmap on combine
	Normalize();

	// assuming both lightmap and AO are already dilated
	// final version
	image.lock();


	float gamma = 1.0f / m_Settings_Gamma;

	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			float ao = m_Image_AO.GetItem(x, y);
			FColor lum = m_Image_L.GetItem(x, y);

			// combined
			FColor f;
			switch (m_Settings_BakeMode)
			{
			case LMBAKEMODE_LIGHTMAP:
				{
					f = lum;
				}
				break;
			case LMBAKEMODE_AO:
				{
					f.Set(ao);
				}
				break;
			default:
				{
					FColor mid = lum * ao;

					if (m_Settings_Light_AO_Ratio < 0.5f)
					{
						float r = m_Settings_Light_AO_Ratio / 0.5f;
						f.Set((1.0f - r) * ao);
						f += mid * r;
					}
					else
					{
						float r = (m_Settings_Light_AO_Ratio-0.5f) / 0.5f;
						f =  mid * (1.0f - r);
						f += lum * r;
					}
				}
				break;
			}

			// gamma correction
			if (!m_Settings_CombinedIsHDR)
			{
				f.r = powf(f.r, gamma);
				f.g = powf(f.g, gamma);
				f.b = powf(f.b, gamma);
			}

			Color col;
			col = Color(f.r, f.g, f.b, 1);

			// new... RGBM .. use a multiplier in the alpha to get increased dynamic range!
			//ColorToRGBM(col);

			image.set_pixel(x, y, col);
		}
	}

	image.unlock();
}


void LightMapper_Base::WriteOutputImage_AO(Image &image)
{
	Dilate<float> dilate;
	dilate.DilateImage(m_Image_AO, m_Image_ID_p1, 256);

	// final version
	image.lock();

	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			const float * pf = m_Image_AO.Get(x, y);
			assert (pf);
			float f = *pf;

			// gamma correction
			if (!m_Settings_AmbientIsHDR)
			{
				float gamma = 1.0f / 2.2f;
				f = powf(f, gamma);
			}

			Color col;
			col = Color(f, f, f, 1);


			// debug mark the dilated pixels
//#define MARK_AO_DILATED
#ifdef MARK_AO_DILATED
			if (!m_Image_ID_p1.GetItem(x, y))
			{
				col = Color(1.0f, 0.33f, 0.66f, 1);
			}
#endif
			image.set_pixel(x, y, col);
		}
	}

	image.unlock();
}



void LightMapper_Base::WriteOutputImage_Lightmap(Image &image)
{
	Dilate<FColor> dilate;
	dilate.DilateImage(m_Image_L, m_Image_ID_p1, 256);

	// test
//	int test_size = 7;
//	LightImage<float> imf;
//	imf.Create(test_size, test_size);
//	LightImage<uint32_t> imi;
//	imi.Create(test_size, test_size);
//	imi.GetItem(3, 3) = 255;
//	dilate.DilateImage(imf, imi);

//	Normalize();

	////
	// write some debug
//#define LLIGHTMAPPER_OUTPUT_TRIIDS
#ifdef LLIGHTMAPPER_OUTPUT_TRIIDS
	output_image.lock();
	Color cols[1024];
	for (int n=0; n<m_Scene.GetNumTris(); n++)
	{
		if (n == 1024)
			break;

		cols[n] = Color(Math::randf(), Math::randf(), Math::randf(), 1.0f);
	}
	cols[0] = Color(0, 0, 0, 1.0f);

	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			int coln = m_Image_ID_p1.GetItem(x, y) % 1024;

			output_image.set_pixel(x, y, cols[coln]);
		}
	}

	output_image.unlock();
	output_image.save_png("tri_ids.png");
#endif

	// final version
	image.lock();

	for (int y=0; y<m_iHeight; y++)
	{
		for (int x=0; x<m_iWidth; x++)
		{
			FColor f = *m_Image_L.Get(x, y);

			// gamma correction
			if (!m_Settings_LightmapIsHDR)
			{
				float gamma = 1.0f / 2.2f;
				f.r = powf(f.r, gamma);
				f.g = powf(f.g, gamma);
				f.b = powf(f.b, gamma);
			}

			Color col;
			col = Color(f.r, f.g, f.b, 1);


			// debug mark the dilated pixels
//#define MARK_DILATED
#ifdef MARK_DILATED
			if (!m_Image_ID_p1.GetItem(x, y))
			{
				col = Color(1.0f, 0.33f, 0.66f, 1);
			}
#endif
			//			if (m_Image_ID_p1.GetItem(x, y))
//			{
//				output_image.set_pixel(x, y, Color(f, f, f, 255));
//			}
//			else
//			{
//				output_image.set_pixel(x, y, Color(0, 0, 0, 255));
//			}

			// visual cuts
//			if (m_Image_Cuts.GetItem(x, y).num)
//			{
//				col = Color(1.0f, 0.33f, 0.66f, 1);
//			}
//			else
//			{
//				col = Color(0, 0, 0, 1);
//			}

			// visualize concave
//			const MiniList_Cuts &cuts = m_Image_Cuts.GetItem(x, y);
//			if (cuts.num == 2)
//			{
//				if (cuts.convex)
//				{
//					col = Color(1.0f, 0, 0, 1);
//				}
//				else
//				{
//					col = Color(0, 0, 1, 1);
//				}
//			}

			image.set_pixel(x, y, col);
		}
	}

	image.unlock();
}


