#pragma once

#include "lvector.h"

namespace LM
{

template <class T> class LightImage
{
public:
	LightImage()
	{
		m_uiWidth = 0;
		m_uiHeight = 0;
		m_uiNumPixels = 0;
	}

	void Create(uint32_t width, uint32_t height, bool blank = true)
	{
		m_uiWidth = width;
		m_uiHeight = height;
		m_uiNumPixels = width * height;
		m_Pixels.resize(m_uiNumPixels);
		if (blank)
			Blank();
	}

	T * Get(uint32_t p)
	{
		if (p < m_uiNumPixels)
			return &m_Pixels[p];
		return 0;
	}

	const T * Get(uint32_t p) const
	{
		if (p < m_uiNumPixels)
			return &m_Pixels[p];
		return 0;
	}

	T * Get(uint32_t x, uint32_t y)
	{
		uint32_t p = GetPixelNum(x, y);
		if (p < m_uiNumPixels)
			return &m_Pixels[p];
		return 0;
	}
	const T * Get(uint32_t x, uint32_t y) const
	{
		uint32_t p = GetPixelNum(x, y);
		if (p < m_uiNumPixels)
			return &m_Pixels[p];
		return 0;
	}

	const T &GetItem(uint32_t x, uint32_t y) const {return *Get(x, y);}
	T &GetItem(uint32_t x, uint32_t y) {return *Get(x, y);}

	void Fill(const T &val)
	{
		for (uint32_t n=0; n<m_uiNumPixels; n++)
		{
			m_Pixels[n] = val;
		}
	}
	void Blank()
	{
		T val;
		memset(&val, 0, sizeof (T));
		Fill(val);
	}

	uint32_t GetWidth() const {return m_uiWidth;}
	uint32_t GetHeight() const {return m_uiHeight;}
	uint32_t GetNumPixels() const {return m_uiNumPixels;}
	bool IsWithin(int x, int y) const
	{
		if (x < 0) return false;
		if (y < 0) return false;
		if (x >= (int) m_uiWidth) return false;
		if (y >= (int) m_uiHeight) return false;
		return true;
	}

	void CopyTo(LightImage<T> &dest) const
	{
		assert (dest.GetNumPixels() == GetNumPixels());

		for (unsigned int n=0; n<m_uiNumPixels; n++)
		{
			dest.m_Pixels[n] = m_Pixels[n];
		}
	}

private:
	uint32_t GetPixelNum(uint32_t x, uint32_t y) const
	{
		return (y * m_uiWidth) + x;
	}

	uint32_t m_uiWidth;
	uint32_t m_uiHeight;
	uint32_t m_uiNumPixels;
	LVector<T> m_Pixels;
};

}
