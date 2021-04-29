#pragma once

#include "core/os/memory.h"
#include "llightimage.h"
#include "llighttypes.h"

namespace LM {

template <class T>
class Convolution {
public:
	enum {
		NUM_LINES = 3
	};

	Convolution() {
		for (int n = 0; n < NUM_LINES; n++) {
			m_pLines_Source[n] = nullptr;
			m_pLines_Dest[n] = nullptr;
		}
	}

	~Convolution() {
		Destroy();
	}

	void Destroy() {
		for (int n = 0; n < NUM_LINES; n++) {
			if (m_pLines_Source[n]) {
				memdelete_arr(m_pLines_Source[n]);
				m_pLines_Source[n] = nullptr;
			}
			if (m_pLines_Dest[n]) {
				memdelete_arr(m_pLines_Dest[n]);
				m_pLines_Dest[n] = nullptr;
			}
		}
	}

	void CreateLines() {
		CRASH_COND(m_pLines_Source[0]);
		CRASH_COND(!m_iLineLength_Expanded);

		for (int n = 0; n < NUM_LINES; n++) {
			m_pLines_Source[n] = memnew_arr(T, m_iLineLength_Expanded);
			m_pLines_Dest[n] = memnew_arr(T, m_iLineLength_Expanded);
		}
	}

	void CopySourceLine(int source_y, int dest_y) {
		const T *pSource = m_pLines_Source[source_y];
		T *pDest = m_pLines_Source[dest_y];

		for (int x = 0; x < m_iLineLength_Expanded; x++) {
			pDest[x] = pSource[x];
		}
	}

	void LoadLine(int source_y, int dest_y) {
		const T *pSource = &m_pImage->GetItem(0, source_y);
		T *pDest = m_pLines_Source[dest_y];

		// leave a pixel gap at start and end for filtering
		for (int x = 0; x < m_iLineLength_Orig; x++) {
			pDest[x + 1] = pSource[x];
		}

		// wrap start and end
		pDest[0] = pDest[1];
		pDest[m_iLineLength_Expanded - 1] = pDest[m_iLineLength_Expanded - 2];
	}

	void SaveLine(int source_y, int dest_y) {
		const T *pSource = m_pLines_Dest[source_y];
		pSource++;
		T *pDest = &m_pImage->GetItem(0, dest_y);

		for (int x = 0; x < m_iLineLength_Orig; x++) {
			pDest[x] = pSource[x];
		}
	}

	void ZeroPixel(float &f) { f = 0.0f; }
	void ZeroPixel(FColor &c) { c.Set(0.0f); }
	//	void ZeroPixel(Color &c) { c = Color(0.0f, 0.0f, 0.0f, 1.0f);}

	void AdjustCentre(float &centre, float average) {
		float diff = average - centre;
		if (fabsf(diff) < m_fTolerance) {
			centre = m_fAmount * average + (centre * (1.0f - m_fAmount));
		}
	}

	void AdjustCentre(FColor &centre, FColor average) {
		FColor diff = average;
		diff -= centre;

		if ((fabsf(diff.r) < m_fTolerance) &&
				(fabsf(diff.g) < m_fTolerance) &&
				(fabsf(diff.b) < m_fTolerance)) {
			centre *= 1.0f - m_fAmount;
			average *= m_fAmount;
			centre += average;
		}
	}

	//	void AdjustCentre(Color &centre, Color average) {
	//		Color diff = average;
	//		diff -= centre;

	//		if ((fabsf(diff.r) < m_fTolerance) &&
	//				(fabsf(diff.g) < m_fTolerance) &&
	//				(fabsf(diff.b) < m_fTolerance)) {
	//			centre *= 1.0f - m_fAmount;
	//			average *= m_fAmount;
	//			centre += average;
	//		}

	//		centre.a = 1.0f;
	//	}

	void ConvolvePixel(int x, int top, int mid, int bot) {
		T total;
		ZeroPixel(total);

		total += m_pLines_Source[top][x - 1];
		total += m_pLines_Source[top][x];
		total += m_pLines_Source[top][x + 1];

		total += m_pLines_Source[mid][x - 1];
		total += m_pLines_Source[mid][x + 1];

		total += m_pLines_Source[bot][x - 1];
		total += m_pLines_Source[bot][x];
		total += m_pLines_Source[bot][x + 1];

		T centre = m_pLines_Source[mid][x];

		total += centre;

		total *= 1.0f / 9.0f;

		AdjustCentre(centre, total);

		// save
		m_pLines_Dest[mid][x] = centre;
	}

	void Run(LightImage<T> &image, float tolerance = 0.2f, float amount = 1.0f) {
		m_pImage = &image;
		m_fTolerance = tolerance;
		m_fAmount = amount;

		// noop?
		if (amount <= 0.001f)
			return;

		if (image.GetWidth() < 2) return;
		if (image.GetHeight() < 2) return;

		m_iLineLength_Orig = image.GetWidth();
		m_iLineLength_Expanded = m_iLineLength_Orig + 2;

		CreateLines();

		int height = image.GetHeight();

		// preload for first line
		LoadLine(0, 1);
		LoadLine(1, 2);
		CopySourceLine(1, 0);

		int top = 0;
		int mid = 1;
		int bot = 2;

		for (int y = 1; y < height + 1; y++) {
			// convolve
			for (int x = 1; x < m_iLineLength_Expanded - 1; x++) {
				ConvolvePixel(x, top, mid, bot);
			}

			// save out top
			if (y > 1)
				SaveLine(top, y - 2);

			// shuffle
			int temp = top;
			top = mid;
			mid = bot;
			bot = temp;

			// load bot
			if ((y + 1) < height) {
				LoadLine(y + 1, bot);
			} else {
				// copy last line
				CopySourceLine(mid, bot);
			}
		}

		Destroy();
	}

private:
	LightImage<T> *m_pImage;
	T *m_pLines_Source[NUM_LINES];
	T *m_pLines_Dest[NUM_LINES];

	int m_iLineLength_Orig;
	int m_iLineLength_Expanded;
	float m_fTolerance;
	float m_fAmount;
};

} // namespace LM
