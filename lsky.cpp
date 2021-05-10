#include "lsky.h"

using namespace LM;

bool LSky::load_sky(String p_filename, float p_blur, int p_tex_size) {
	if (p_filename == "") {
		_active = false;
		return false;
	}

	Ref<Image> im;
	im.instance();
	if (im->load(p_filename) != OK) {
		_active = false;
		WARN_PRINT("ERROR loading " + p_filename + ", ignoring sky.");
		return false;
	}

	// shrink
	p_tex_size *= p_tex_size;
	while ((im->get_width() * im->get_height()) > p_tex_size) {
		im->shrink_x2();
	}

	// convert to our image format
	int w = im->get_width();
	int h = im->get_height();

	_texture.Create(w, h);

	im->lock();

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			Color c = im->get_pixel(x, y);
			_texture.GetItem(x, y).Set(c.r, c.g, c.b);
		}
	}

	im->unlock();

	// just using the x dimension here, we could use longest axis etc...
	_blur_pixels = p_blur * w;

	// blur
	_blur();

	// normalize and brightness
	_adjust_brightness();

	// debug save
	//_debug_save();

	_active = true;

	return true;
}

void LSky::unload_sky() {
	_texture.Reset();
	_active = false;
}

void LSky::_blur_horz(FColor *p_scan, const float *p_curve, int p_num_scan, LightImage<FColor> &p_output, int p_y) {
	int w = _texture.GetWidth();

	// horizontal
	// prime
	int count = 0;
	FColor average;
	average.Set(0.0);

	for (int x = -(_blur_pixels); x < _blur_pixels; x++) {
		p_scan[count] = _blur_read_x(x, p_y);

		// average the whole scan except the last
		average += p_scan[count];
		count++;
	}
	p_scan[p_num_scan - 1].Set(0.0);

	// next write
	int scan_pointer = p_num_scan - 1;

	for (int x = 0; x < w; x++) {
		// remove one from scan
		average -= p_scan[scan_pointer];

		// add one to scan
		p_scan[scan_pointer] = _blur_read_x(x + _blur_pixels, p_y);
		average += p_scan[scan_pointer];

		// move on scan pointer
		scan_pointer++;
		scan_pointer %= p_num_scan; // wraparound

		// write output
		//p_output.GetItem(x, p_y) = average / p_num_scan;
		p_output.GetItem(x, p_y) = _gaussian_blur(p_scan, p_num_scan, scan_pointer, p_curve);
	}
}

void LSky::_blur_vert(FColor *p_scan, const float *p_curve, int p_num_scan, LightImage<FColor> &p_output, int p_x) {
	int h = _texture.GetHeight();

	// horizontal
	// prime
	int count = 0;
	FColor average;
	average.Set(0.0);

	for (int y = -(_blur_pixels); y < _blur_pixels; y++) {
		p_scan[count] = _blur_read_y(p_x, y);

		// average the whole scan except the last
		average += p_scan[count];
		count++;
	}
	p_scan[p_num_scan - 1].Set(0.0);

	// next write
	int scan_pointer = p_num_scan - 1;

	for (int y = 0; y < h; y++) {
		// remove one from scan
		average -= p_scan[scan_pointer];

		// add one to scan
		p_scan[scan_pointer] = _blur_read_y(p_x, y + _blur_pixels);
		average += p_scan[scan_pointer];

		// move on scan pointer
		scan_pointer++;
		scan_pointer %= p_num_scan; // wraparound

		// write output
		//		p_output.GetItem(p_x, y) = average / p_num_scan;
		p_output.GetItem(p_x, y) = _gaussian_blur(p_scan, p_num_scan, scan_pointer, p_curve);
	}
}

void LSky::_blur() {
	int w = _texture.GetWidth();
	int h = _texture.GetHeight();

	LightImage<FColor> output;
	output.Create(w, h);

	int num_scan = (_blur_pixels * 2) + 1;

	FColor *scan = (FColor *)alloca(num_scan * sizeof(FColor));
	float *curve = (float *)alloca(num_scan * sizeof(float));

	_create_curve(curve, num_scan);

	for (int y = 0; y < h; y++) {
		_blur_horz(scan, curve, num_scan, output, y);
	}

	// copy the intermediate output to the texture
	output.CopyTo(_texture);

	for (int x = 0; x < w; x++) {
		_blur_vert(scan, curve, num_scan, output, x);
	}

	// copy the intermediate output to the texture
	output.CopyTo(_texture);
}

void LSky::_create_curve(float *p_curve, int p_num_scan) {
	int centre = p_num_scan / 2;
	int range = centre;

	float standard_deviation = range * 0.15f;

	float total = 0.0f;

	for (int n = 0; n < p_num_scan; n++) {
		// use gaussian function
		float val = _normal_pdf(n, centre, standard_deviation);
		p_curve[n] = val;
		total += val;
	}

	// normalize curve
	for (int n = 0; n < p_num_scan; n++) {
		p_curve[n] /= total;
	}
}

void LSky::_adjust_brightness() {
	int w = _texture.GetWidth();
	int h = _texture.GetHeight();

	// normalize at the same time, to an average per pixel value
	FColor total;
	total.Set(0.0);

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			total += _texture.GetItem(x, y);
		}
	}

	// lightness per pixel
	float lightness = MAX(MAX(total.r, total.g), total.b);
	lightness /= w * h;

	// we want to bump this up to an average of say .. 0.5?
	float brightness = (0.25f) / lightness;

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			_texture.GetItem(x, y) *= brightness;
		}
	}
}

void LSky::_debug_save() {
	int w = _texture.GetWidth();
	int h = _texture.GetHeight();

	Ref<Image> im;
	im.instance();

	im->create(w, h, false, Image::FORMAT_RGBF);

	im->lock();

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			const FColor &c = _texture.GetItem(x, y);
			im->set_pixel(x, y, Color(c.r, c.g, c.b, 1.0));
		}
	}

	im->unlock();

	im->save_exr("blurred.exr", false);
}

void LSky::read_sky(const Vector3 &ptDir, FColor &color) {
	//direction = parameters.environment_transform.xform_inv(direction);
	Vector3 dir = ptDir;
	real_t dir_y = CLAMP(dir.y, 0.0, 1.0);
	Vector2 st = Vector2(Math::atan2(dir.z, dir.x), Math::acos(dir_y));

	// already clamped
	//	if (Math::is_nan(st.y)) {
	//		st.y = direction.y > 0.0 ? 0.0 : Math_PI;
	//	}

	st.x += Math_PI;
	st /= Vector2(Math_TAU, Math_PI);
	st.x = Math::fmod(st.x + 0.75, 1.0);
	FColor c = _bilinear_sample(st, false, true);
	//	color += throughput * Vector3(c.r, c.g, c.b) * c.a;

	color += c;
	//	color.r += c.r;
	//	color.g += c.g;
	//	color.b += c.b;
}

FColor LSky::_bilinear_sample(const Vector2 &p_uv, bool p_clamp_x, bool p_clamp_y) {
	int width = _texture.GetWidth();
	int height = _texture.GetHeight();

	Vector2 uv;
	uv.x = p_clamp_x ? p_uv.x : Math::fposmod(p_uv.x, 1.0f);
	uv.y = p_clamp_y ? p_uv.y : Math::fposmod(p_uv.y, 1.0f);

	float xf = uv.x * width;
	float yf = uv.y * height;

	int xi = (int)xf;
	int yi = (int)yf;

	FColor texels[4];
	for (int i = 0; i < 4; i++) {
		int sample_x = xi + i % 2;
		int sample_y = yi + i / 2;

		sample_x = CLAMP(sample_x, 0, width - 1);
		sample_y = CLAMP(sample_y, 0, height - 1);

		texels[i] = _texture.GetItem(sample_x, sample_y);
	}

	float tx = xf - xi;
	float ty = yf - yi;

	FColor c;
	c.Set(0);
	//	for (int i = 0; i < 4; i++) {
	FColor a = texels[0];
	a.Lerp(texels[1], tx);

	FColor b = texels[2];
	b.Lerp(texels[3], tx);

	c = a;
	c.Lerp(b, ty);

	//		c[i] = Math::lerp(Math::lerp(texels[0][i], texels[1][i], tx), Math::lerp(texels[2][i], texels[3][i], tx), ty);
	//	}
	return c;
}
