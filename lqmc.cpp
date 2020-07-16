#include "lqmc.h"

namespace LM
{

void QMC::Create(int num_samples)
{
	m_CurrentVariation = 0;
	m_Group.m_Samples.resize(num_samples);

	for (int n=0; n<NUM_VARIATIONS; n++)
	{
		GenerateVariation(m_Group, n);
	}

}


void QMC::RandomUnitDir(Vector3 &dir) const
{
	while (true)
	{
		dir.x = Math::random(-1.0f, 1.0f);
		dir.y = Math::random(-1.0f, 1.0f);
		dir.z = Math::random(-1.0f, 1.0f);

		float l = dir.length();
		if (l > 0.001f)
		{
			dir /= l;
			return;
		}
	}
}

void QMC::GenerateVariation(Group &group, int var)
{
	int nSamples = m_Group.m_Samples.size();

	const int nTestSamples = 8;
	Vector3 tests[nTestSamples];
	float dots[nTestSamples];

	for (int s=0; s<nSamples; s++)
	{
		// choose some random test directions
		for (int n=0; n<nTestSamples; n++)
		{
			RandomUnitDir(tests[n]);
		}

		// find the best
		int best_t = 0;

		// blank dots
		for (int d=0; d<nTestSamples; d++)
		{
			dots[d] = -1.0f;
		}

		// compare to all the rest
		for (int n=0; n<s-1; n++)
		{
			const Vector3 &existing_dir = group.m_Samples[n].dir[var];

			for (int t=0; t<nTestSamples; t++)
			{
				float dot = existing_dir.dot(tests[t]);

				if (dot > dots[t])
				{
					dots[t] = dot;
				}
			}
		}

		// find the best test
		for (int t=0; t<nTestSamples; t++)
		{
			if (dots[t] < dots[best_t])
			{
				best_t = t;
			}
		}

		// now we know the best, add it to the list
		group.m_Samples[s].dir[var] = tests[best_t];
	} // for variation
}

int QMC::RandomVariation() const
{
	return Math::rand() % NUM_VARIATIONS;
}

int QMC::GetNextVariation(int previous) const
{
	int next = previous +1;
	if (next >= NUM_VARIATIONS)
		next = 0;

	return next;
}


//void QMC::QMCRandomUnitDir(Vector3 &dir, int count, int variation)
//{
//	if (count == 0)
//	{
//		m_CurrentVariation++;

//		// wraparound
//		if (m_CurrentVariation >= NUM_VARIATIONS)
//			m_CurrentVariation = 0;
//	}

//	dir = m_Group.m_Samples[count].dir[m_CurrentVariation];
//}


} // namespace
