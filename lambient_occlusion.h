#pragma once

#include "lraybank.h"

namespace LM {

class AmbientOcclusion : public RayBank
{
public:

protected:
	struct AOSample
	{
		Vector2 uv;
		uint32_t tri_id;
		Vector3 pos;
		Vector3 normal;
	};


	enum {MAX_COMPLEX_AO_TEXEL_SAMPLES = 16};

	void ProcessAO();
	void ProcessAO_LineMT(uint32_t y_offset, int y_section_start);
	void ProcessAO_Texel(int tx, int ty, int qmc_variation);
	float CalculateAO(int tx, int ty, int qmc_variation, const MiniList &ml);
	float CalculateAO_Complex(int tx, int ty, int qmc_variation, const MiniList &ml);

private:
	int AO_FindSamplePoints(int tx, int ty, const MiniList &ml, AOSample samples[MAX_COMPLEX_AO_TEXEL_SAMPLES]);


	void AO_RandomTexelSample(Vector2 &st, int tx, int ty, int n) const
	{
		if (n)
		{
			st.x = Math::randf() + tx;
			st.y = Math::randf() + ty;
		}
		else
		{
			// fix first to centre of texel
			st.x = 0.5f + tx;
			st.y = 0.5f + ty;
		}

		// has to be ranged 0 to 1
		st.x /= m_iWidth;
		st.y /= m_iHeight;
	}

	bool AO_FindTexelTriangle(const MiniList &ml, const Vector2 &st, uint32_t &tri_inside, Vector3 &bary) const
	{
		// barycentric coords.
		const UVTri * pUVTri;

		for (uint32_t i=0; i<ml.num; i++)
		{
			tri_inside = m_TriIDs[ml.first + i];
			pUVTri = &m_Scene.m_UVTris[tri_inside];

			// within?
			pUVTri->FindBarycentricCoords(st, bary);
			if (BarycentricInside(bary))
			{
				return true;
			}
		}

		// not inside any triangles
		return false;
	}

	void AO_RandomQMCDirection(Vector3 &dir, const Vector3 &ptNormal, int n, int qmc_variation) const
	{
		Vector3 push = ptNormal * m_Settings_SurfaceBias; //0.005f;

		m_QMC.QMCRandomUnitDir(dir, n, qmc_variation);

		// clip?
		float dot = dir.dot(ptNormal);
		if (dot < 0.0f)
			dir = -dir;

		// prevent parallel lines
		dir += push;

		// collision detect
		dir.normalize();
	}

	void AO_RandomQMCRay(Ray &r, const Vector3 &ptNormal, int n, int qmc_variation) const
	{
		Vector3 push = ptNormal * m_Settings_SurfaceBias; // 0.005f;

		// push ray origin
		r.o += push;

//		RandomUnitDir(r.d);
		m_QMC.QMCRandomUnitDir(r.d, n, qmc_variation);

		// clip?
		float dot = r.d.dot(ptNormal);
		if (dot < 0.0f)
		{
			// make dot always positive for calculations as to the weight given to this hit
			//dot = -dot;
			r.d = -r.d;
		}

		// prevent parallel lines
		r.d += push;
	//	r.d = ptNormal;

		// collision detect
		r.d.normalize();
	}

	void AO_FindSamplesRandomRay(Ray &r, const Vector3 &ptNormal) const
	{
		//Vector3 push = ptNormal * m_Settings_AO_ReverseBias; //0.005f;
		Vector3 push = ptNormal * m_Settings_SurfaceBias; //0.005f;

		// push ray origin
		r.o -= push;

		RandomUnitDir(r.d);

		// clip?
		float dot = r.d.dot(ptNormal);
		if (dot < 0.0f)
		{
			// make dot always positive for calculations as to the weight given to this hit
			//dot = -dot;
			r.d = -r.d;
		}

		// prevent parallel lines
		r.d += push;
	//	r.d = ptNormal;

		// collision detect
		r.d.normalize();
	}

};

} // namespace
