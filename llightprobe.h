#pragma once

#include "lvector.h"


namespace LM
{

class LightMapper;

class LightProbe
{
public:
	void Save(FileAccess *f);
	void Save_Secondary(FileAccess *f);
	struct Contribution
	{
		uint32_t light_id;
		float power; // linear
	};

//	LightProbe()
//	{
//		m_fLight = 0.0f;
//	}
//	float m_fLight;
	LVector<Contribution> m_Contributions;
};


class LightProbes
{
public:
	void Create(LightMapper &lm);



private:
	void Save(const char * pszFilename);
	void CalculateProbe(const Vec3i &pt);
	float CalculatePower(const LightMapper_Base::LLight &light, float dist, const Vector3 &pos) const;

	void Save_Vector3(FileAccess *f, const Vector3 &v);

	LightProbe * GetProbe(const Vec3i &pt)
	{
		unsigned int ui = GetProbeNum(pt);
		if (ui < m_Probes.size())
		{
			return &m_Probes[ui];
		}
		return nullptr;
	}

	void GetProbePosition(const Vec3i &pt, Vector3 &pos) const
	{
		pos = m_ptMin;
		pos.x += pt.x * m_VoxelSize.x;
		pos.y += pt.y * m_VoxelSize.y;
		pos.z += pt.z * m_VoxelSize.z;
	}

	unsigned int GetProbeNum(const Vec3i &pt) const
	{
		return (pt.z * m_DimXTimesY) + (pt.y * m_Dims.x) + pt.x;
	}

	// number of probes on each axis
	Vec3i m_Dims;
	int m_DimXTimesY;

	// world size of each voxel separating probes
	Vector3 m_VoxelSize;

	// bottom left probe point, can be used to derive the other points using the voxel size
	Vector3 m_ptMin;

	LVector<LightProbe> m_Probes;

	LightMapper * m_pLightMapper;
};


} // namespace
