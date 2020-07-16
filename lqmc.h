#pragma once

#include "core/math/vector3.h"
#include "lvector.h"

namespace LM
{

class QMC
{
	enum
	{
		NUM_VARIATIONS = 1024
	};

	struct Sample
	{
		Vector3 dir[NUM_VARIATIONS];
	};

	struct Group
	{
		LVector<Sample> m_Samples;
	};

public:
	void Create(int num_samples);
	void QMCRandomUnitDir(Vector3 &dir, int count, int variation) const
	{
		dir = m_Group.m_Samples[count].dir[variation];
	}
	int GetNextVariation(int previous) const;
	int RandomVariation() const;

private:
	void GenerateVariation(Group &group, int var);
	void RandomUnitDir(Vector3 &dir) const;

	int m_CurrentVariation;
	Group m_Group;
};


} // namespace
