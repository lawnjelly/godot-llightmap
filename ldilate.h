#pragma once

#include "llightimage.h"
#include "lvector.h"
#include "llighttypes.h"
#include <limits.h>

//#define DILATE_VERBOSE

namespace LM
{

template <class T> class Dilate
{
public:
	bool DilateImage(LightImage<T> &im, const LightImage <uint32_t> &orig_mask, unsigned int uiMaxDist)
	{
		print_line("Dilating image");

		int w = im.GetWidth();
		int h = im.GetHeight();

		// construct mask
		LightImage<uint8_t> mask;
		mask.Create(w, h);
		for (int y=0; y<h; y++)
		{
			for (int x=0; x<w; x++)
			{
				if (*orig_mask.Get(x, y))
					*mask.Get(x, y) = 255;
			}
		}

		LightImage<Vec2_i16> source_pts;
		Run(mask, &source_pts, uiMaxDist);

		ConvertFinalImage(im);
		return true;
	}

	// mem vars
	LightImage<uint8_t> m_Done;
	LightImage<Vec2_i16> * m_pSource_pts;
	LVector<Vec2_i16> m_Actives;
	int m_iWidth;
	int m_iHeight;

	void Run(LightImage<uint8_t> &mask, LightImage<Vec2_i16> * pSource_pts, unsigned int uiMaxDist)
	{
		m_iWidth = mask.GetWidth();
		m_iHeight = mask.GetHeight();
		int w = m_iWidth;
		int h = m_iHeight;

		if ((w * h) == 0)
			return;

		// note which pixels have been done (activated)
		m_Done.Create(w, h);
		m_Done.Blank();

		// source coordinates for best source
		m_pSource_pts = pSource_pts;
		m_pSource_pts->Create(w, h);

		// the active points
		m_Actives.reserve((w*h) / 8); // just a rough number to start with, the vector can grow

		// first find the edges
		FindEdges(mask);

		// required number of iterations
		for (unsigned int i=0; i<uiMaxDist; i++)
		{
			// this must be recorded before starting, as we will be adding to the list
			unsigned int uiNumActives = m_Actives.size();

			for (unsigned int n=0; n<uiNumActives; n++)
			{
				ProcessPixel(n, mask);
			}

			// delete the previous active items
			m_Actives.delete_items_first(uiNumActives);

			PrintActives();
		}

	}

	void ConvertFinalImage(LightImage<T> &image)
	{
		for (unsigned int n=0; n<m_Done.GetNumPixels(); n++)
		{
			if (!m_Done.Get(n))
				continue;

			// get the x y
			int y = n / m_iWidth;
			int x = n % m_iWidth;

			const Vec2_i16 &scoords = *m_pSource_pts->Get(x, y);
			if ((scoords.x != 0) || (scoords.y != 0))
				image.GetItem(x, y) = image.GetItem(scoords.x, scoords.y);
		}
	}

	unsigned int SquareLength(const Vector2i &v) const {return (v.x * v.x) + (v.y * v.y);}

	void ProcessPixel(unsigned int uiPixel, LightImage<uint8_t> &mask)
	{
		const Vec2_i16 &loc = m_Actives[uiPixel];
		int iRange = 3;

		Rect2i rect = Rect2i(loc.x, loc.y, 0, 0);
		rect = rect.grow(iRange);
		rect.size.x += 1;
		rect.size.y += 1;

		// clip
		rect = rect.clip(Rect2i(0, 0, m_iWidth, m_iHeight));
		if ((rect.size.x == 0) && (rect.size.y == 0))
			return;

		// get loc quick format
		Vector2i iloc = Vector2i(loc.x, loc.y);

		Vec2_i16 bestcoord;
		unsigned int sl_best = UINT_MAX;

		int bottom = rect.position.y + rect.size.y;
		int right = rect.position.x + rect.size.x;

		for (int y=rect.position.y; y<bottom; y++)
		{
			for (int x=rect.position.x; x<right; x++)
			{
				//if ((x == iloc.x) && (y == iloc.y))
				//	continue;

				Vec2_i16 coord = *m_pSource_pts->Get(x, y);

				// ignore
				if (coord.IsZero())
					continue;

				// calculate distance from HERE to the referenced point
				Vector2i offset = Vector2i(coord.x, coord.y);
				offset -= iloc;
				unsigned int sl = SquareLength(offset);

				// is this the best fit?
				if (sl < sl_best)
				{
					sl_best = sl;
					bestcoord = coord;
				}
			}
		}

		// assert that we have found a coord
		assert (sl_best < UINT_MAX);
		*m_pSource_pts->Get(loc.x, loc.y) = bestcoord;

		// add suitable neighbours to active list
		AddActiveNeighs(loc.x, loc.y, mask);
	}

	void FindEdges(LightImage<uint8_t> &mask)
	{
		int w = mask.GetWidth();
		int h = mask.GetHeight();

		for (int y=0; y<h; y++)
		{
			for (int x=0; x<w; x++)
			{
				// assess this point
				if (mask.GetItem(x, y) || m_Done.GetItem(x, y))
//				if (mask.GetItem(x, y))
					continue;

				// do we have a used neighbour?
				//bool bSuitable = false;

//#define TEST_NEIGH(a, b) if (image.GetItem(a, b)) {m_Actives.Add(Prim::CoPoint2S(a, b)); bSuitable = true; }
#define TEST_NEIGH(a, b) if (mask.IsWithin(a, b) && mask.GetItem(a, b)) {m_pSource_pts->GetItem(x, y).Set(a, b);\
goto suitable; }
//image.GetItem(x, y) = image.GetItem(a, b);

				// closest 4
				TEST_NEIGH(x, y-1)
				TEST_NEIGH(x-1, y)
				TEST_NEIGH(x+1, y)
				TEST_NEIGH(x, y+1)

				// then the diagnals
				TEST_NEIGH(x-1, y-1)
				TEST_NEIGH(x+1, y-1)
				TEST_NEIGH(x-1, y+1)
				TEST_NEIGH(x+1, y+1)
#undef TEST_NEIGH

				// not suitable, just continue to next x
				continue;

suitable:
#ifdef DILATE_VERBOSE
				print_line("suitable : " + itos(x) + ", " + itos(y));
#endif
				m_Done.GetItem(x, y) = 255;
				//AddActiveNeighs(x, y, mask);
				m_Actives.push_back(Vec2_i16(x, y));
			} // for x
		} // for y

		PrintActives();
	}

	void CheckNeigh(int a, int b, LightImage<uint8_t> &mask)
	{
		bool bWithin = mask.IsWithin(a, b);
		if (!bWithin)
			return;

		bool bSourcePtSet = m_pSource_pts->GetItem(a, b).IsNonZero();

		if (!bSourcePtSet)
		{
			bool bDone = m_Done.GetItem(a, b);
			bool bMaskSet = mask.GetItem(a, b);
			if ((!bDone)
			&& (!bMaskSet))
			{
#ifdef DILATE_VERBOSE
				print_line("\tadding active " + itos(a) + ", " + itos(b));
#endif
				m_Actives.push_back(Vec2_i16(a, b));
				m_Done.GetItem(a, b) = 255;
			}

		}

	}


	void AddActiveNeighs(int x, int y, LightImage<uint8_t> &mask)
	{
//#define CHECK_NEIGH(a, b) if ((image.IsWithin(a, b)) && (!m_pSource_pts->GetItem(a, b).IsNonZero()))\
//{\
//if ((!m_Done.GetItem(a, b)) && (!image.GetItem(a, b)))\
//m_Actives.push_back(Vec2_i16(a, b));\
//m_Done.GetItem(a, b) = 255;\
//}

				CheckNeigh(x-1, y-1, mask);
				CheckNeigh(x, y-1, mask);
				CheckNeigh(x+1, y-1, mask);
				CheckNeigh(x-1, y, mask);
				CheckNeigh(x+1, y, mask);
				CheckNeigh(x-1, y+1, mask);
				CheckNeigh(x, y+1, mask);
				CheckNeigh(x+1, y+1, mask);

//#undef CHECK_NEIGH
	}

	bool OnActiveList(int x, int y) const
	{
		for (int n=0; n<m_Actives.size(); n++)
		{
			Vec2_i16 a = m_Actives[n];
			if ((a.x == x) && (a.y == y))
				return true;
		}
		return false;
	}

	void PrintActives()
	{
#ifdef DILATE_VERBOSE
		for (int y=0; y<m_iHeight; y++)
		{
			String sz;
			for (int x=0; x<m_iWidth; x++)
			{
				if (m_Done.GetItem(x, y))
					sz += "*";
				else
					sz += ".";

				if (OnActiveList(x, y))
					sz += "a";
				else
					sz += " ";

				Vec2_i16 cd = m_pSource_pts->GetItem(x, y);
				sz += "(" + itos(cd.x) + ", " + itos(cd.y) + ")";
			}
			print_line(sz);
		}

//		print_line (itos(m_Actives.size()) + " Actives:");
//		for (int n=0; n<m_Actives.size(); n++)
//		{
//			Vec2_i16 a = m_Actives[n];
//			print_line("\t\tactive " + itos(a.x) + ", " + itos(a.y));
//		}
#endif
	}

};


} // namespace

