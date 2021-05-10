#pragma once

#include "core/math/vector3.h"
#include "llightimage.h"
#include "llighttypes.h"

namespace LM {

class LSky {
public:
	bool load_sky(String p_filename, float p_blur, int p_tex_size);
	void unload_sky();

	bool is_active() const { return _active; }
	void read_sky(const Vector3 &ptDir, FColor &color);

private:
	FColor _bilinear_sample(const Vector2 &p_uv, bool p_clamp_x, bool p_clamp_y);
	void _blur();
	void _blur_horz(FColor *p_scan, const float *p_curve, int p_num_scan, LightImage<FColor> &p_output, int p_y);
	void _blur_vert(FColor *p_scan, const float *p_curve, int p_num_scan, LightImage<FColor> &p_output, int p_x);
	void _debug_save();
	void _create_curve(float *p_curve, int p_num_scan);
	void _adjust_brightness();

	const FColor &_blur_read_x(int x, int y) const {
		// wrap x
		int w = _texture.GetWidth();
		x += w;
		x %= w;
		return _texture.GetItem(x, y);
	}

	const FColor &_blur_read_y(int x, int y) const {
		// wrap y
		int h = _texture.GetHeight();
		y = CLAMP(y, 0, h - 1);
		return _texture.GetItem(x, y);
	}

	FColor _gaussian_blur(const FColor *p_scan, int p_num_scan, int p_scan_pointer, const float *p_curve) const {
		FColor total;
		total.Set(0.0);

		for (int n = 0; n < p_num_scan; n++) {
			total += p_scan[p_scan_pointer] * p_curve[n];
			p_scan_pointer++;
			p_scan_pointer %= p_num_scan;
		}

		return total;
	}

	// gaussian distribution function
	float _normal_pdf(float x, float m, float s) const {
		static const float inv_sqrt_2pi = 0.3989422804014327;
		float a = (x - m) / s;

		return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
	}

	LightImage<FColor> _texture;
	bool _active = false;
	int _blur_pixels = 0;
};

} // namespace LM
