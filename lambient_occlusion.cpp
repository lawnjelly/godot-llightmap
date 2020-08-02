#include "lambient_occlusion.h"

namespace LM {

void AmbientOcclusion::ProcessAO_Texel(int tx, int ty, int qmc_variation)
{
//	if ((tx == 77) && (ty == 221))
//		print_line("test");

	const MiniList &ml = m_Image_TriIDs.GetItem(tx, ty);
	if (!ml.num)
		return; // no triangles in this UV

	// could be off the image
	float * pfTexel = m_Image_AO.Get(tx, ty);

	float power = 0.0f;

	// simple or complex? (are there cuts on this texel)
	if (1)
//	if (!m_Image_Cuts.GetItem(tx, ty))
	{
		power = CalculateAO(tx, ty, qmc_variation, ml);
	}
	else
	{
		power = CalculateAO_Complex(tx, ty, qmc_variation, ml);
	}

	*pfTexel = power;
}

void AmbientOcclusion::ProcessAO()
{
	if (bake_begin_function) {
		bake_begin_function(m_iHeight);
	}

	// find the max range in voxels. This can be used to speed up the ray trace
	const float range = m_Settings_AO_Range;
	m_Scene.m_Tracer.GetDistanceInVoxels(range, m_Scene.m_VoxelRange);

#define LAMBIENT_OCCLUSION_USE_THREADS
#ifdef LAMBIENT_OCCLUSION_USE_THREADS
//	int nCores = OS::get_singleton()->get_processor_count();
	int nSections = m_iHeight / 64;
	if (!nSections) nSections = 1;
	int y_section_size = m_iHeight / nSections;
	int leftover_start = y_section_size * nSections;

	for (int s=0; s<nSections; s++)
	{
		int y_section_start = s * y_section_size;

		if (bake_step_function) {
			m_bCancel = bake_step_function(y_section_start, String("Process Texels: ") + " (" + itos(y_section_start) + ")");
			if (m_bCancel)
			{
				if (bake_end_function) {
					bake_end_function();
				}
				return;
			}
		}

		thread_process_array(y_section_size, this, &AmbientOcclusion::ProcessAO_LineMT, y_section_start);
	}

	// leftovers
	int nLeftover = m_iHeight - leftover_start;

	if (nLeftover)
		thread_process_array(nLeftover, this, &AmbientOcclusion::ProcessAO_LineMT, leftover_start);
#else

	for (int y=0; y<m_iHeight; y++)
	{
		if ((y % 10) == 0)
		{
			if (bake_step_function) {
				m_bCancel = bake_step_function(y, String("Process Texels: ") + " (" + itos(y) + ")");
				if (m_bCancel)
				{
					if (bake_end_function) {
						bake_end_function();
					}
					return;
				}
			}
		}

		ProcessAO_LineMT(0, y);
	}
#endif

	if (bake_end_function) {
		bake_end_function();
	}
}

void AmbientOcclusion::ProcessAO_LineMT(uint32_t y_offset, int y_section_start)
{
	int ty = y_section_start + y_offset;

	// seed based on the line
	int qmc_variation = m_QMC.RandomVariation();

	for (int x=0; x<m_iWidth; x++)
	{
		qmc_variation = m_QMC.GetNextVariation(qmc_variation);
		ProcessAO_Texel(x, ty, qmc_variation);
	}
}



float AmbientOcclusion::CalculateAO(int tx, int ty, int qmc_variation, const MiniList &ml)
{
	Ray r;
	int nSamples = m_AdjustedSettings.m_AO_Samples;

	// find the max range in voxels. This can be used to speed up the ray trace
	Vec3i voxel_range = m_Scene.m_VoxelRange;

	int nHits = 0;
	int nSamplesInside = 0;

	for (int n=0; n<nSamples; n++)
	{
		// pick a float position within the texel
		Vector2 st;
		AO_RandomTexelSample(st, tx, ty, n);

		// find which triangle on the minilist we are inside (if any)
		uint32_t tri_inside_id;
		Vector3 bary;
		if (!AO_FindTexelTriangle(ml, st, tri_inside_id, bary))
			continue;

		nSamplesInside++;

		// calculate world position ray origin from barycentric
		m_Scene.m_Tris[tri_inside_id].InterpolateBarycentric(r.o, bary);

		// calculate surface normal (should be use plane?)
		Vector3 ptNormal;
		m_Scene.m_TriNormals[tri_inside_id].InterpolateBarycentric(ptNormal, bary);
		//const Vector3 &ptNormal = m_Scene.m_TriPlanes[tri_inside_id].normal;

		// construct a random ray to test
		AO_RandomQMCRay(r, ptNormal, n, qmc_variation);

		// test ray
		if (m_Scene.TestIntersect_Ray(r, m_Settings_AO_Range, voxel_range))
		{
			nHits++;
		}

	} // for samples

	float fTotal = (float) nHits / nSamplesInside;
	fTotal = 1.0f - (fTotal * 1.0f);

	if (fTotal < 0.0f)
		fTotal = 0.0f;

	return fTotal;
}


float AmbientOcclusion::CalculateAO_Complex(int tx, int ty, int qmc_variation, const MiniList &ml)
{
	// first we need to identify some sample points
	AOSample samples[MAX_COMPLEX_AO_TEXEL_SAMPLES];
	int nSampleLocs = AO_FindSamplePoints(tx, ty, ml, samples);

	// if no good samples locs found
	if (!nSampleLocs)
	{
		// set to dilate
		m_Image_ID_p1.GetItem(tx, ty) = 0;
		return 0.5f; // 0.5 could use an intermediate value for texture filtering to look better?
	}

	int sample_counter = 0;
	int nHits = 0;
	Ray r;

	for (int n=0; n<m_Settings_AO_Samples; n++)
	{
		// get the sample to look from
		const AOSample &sample = samples[sample_counter++];

		// wraparound
		if (sample_counter == nSampleLocs)
			sample_counter = 0;

		r.o = sample.pos;

		// construct a random ray to test
		AO_RandomQMCDirection(r.d, sample.normal, n, qmc_variation);

		// test ray
		if (m_Scene.TestIntersect_Ray(r, m_Settings_AO_Range, m_Scene.m_VoxelRange))
		{
			nHits++;
		}
	} // for n

	float fTotal = (float) nHits / m_Settings_AO_Samples;
	fTotal = 1.0f - (fTotal * 1.0f);

	if (fTotal < 0.0f)
		fTotal = 0.0f;

	return fTotal;
}


int AmbientOcclusion::AO_FindSamplePoints(int tx, int ty, const MiniList &ml, AOSample samples[MAX_COMPLEX_AO_TEXEL_SAMPLES])
{
	int samples_found = 0;
	int attempts = m_Settings_AO_Samples + 64;

	// scale number of sample points roughly to the user interface quality
	int num_desired_samples = m_Settings_AO_Samples;
	if (num_desired_samples > MAX_COMPLEX_AO_TEXEL_SAMPLES)
		num_desired_samples = MAX_COMPLEX_AO_TEXEL_SAMPLES;


	Ray r;

	for (int n=0; n<attempts; n++)
	{
		Vector2 st;
		AO_RandomTexelSample(st, tx, ty, 10);

		// find which triangle on the minilist we are inside (if any)
		uint32_t tri_inside_id;
		Vector3 bary;

		if (!AO_FindTexelTriangle(ml, st, tri_inside_id, bary))
			continue; // not inside a triangle, failed this attempt

		// calculate world position ray origin from barycentric
		Vector3 pos;
		m_Scene.m_Tris[tri_inside_id].InterpolateBarycentric(pos, bary);
		r.o = pos;

		// calculate surface normal (should be use plane?)
		const Vector3 &ptNormal = m_Scene.m_TriPlanes[tri_inside_id].normal;

		// construct a random ray to test
		AO_FindSamplesRandomRay(r, ptNormal);

		// test ray
		if (m_Scene.TestIntersect_Ray(r, m_Settings_AO_Range, m_Scene.m_VoxelRange, true))
			continue; // hit, not interested

		// no hit .. we can use this sample!
		AOSample &sample = samples[samples_found++];

		// be super careful with this offset.
		// the find test offsets BACKWARDS to find tris on the floor
		// the actual ambient test offsets FORWARDS to avoid self intersection.
		// the ambient test could also use the backface culling test, but this is slower.
		sample.pos = pos + (ptNormal * m_Settings_SurfaceBias); //0.005f);
		sample.uv = st;
		sample.normal = ptNormal;
		sample.tri_id = tri_inside_id;

		// finished?
		if (samples_found == num_desired_samples)
			return samples_found;
	}

	return samples_found;
}



} // namespace
