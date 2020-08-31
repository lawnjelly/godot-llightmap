#include "llightprobe.h"
#include "llightmapper.h"
#include "llighttypes.h"

namespace LM {


void LightProbes::Create(LightMapper &lm)
{
	return;

	// save locally
	m_pLightMapper = &lm;

	// get the aabb (smallest) of the world bound.
	AABB bb = lm.GetTracer().GetWorldBound_contracted();

	// can't make probes for a non existant world
	if (bb.get_shortest_axis_size() <= 0.0f)
		return;

	// start hard coded
	m_Dims.Set(64, 8, 64);
	m_DimXTimesY = (m_Dims.x * m_Dims.y);

	Vector3 voxel_dims;
	m_Dims.To(voxel_dims);

	// voxel dimensions and start point
	m_VoxelSize = bb.size / voxel_dims;

	// start point is offset by half a voxel
	m_ptMin = bb.position + (m_VoxelSize * 0.5f);

	m_Probes.resize(m_Dims.x * m_Dims.y * m_Dims.z);

	for (int z=0; z<m_Dims.z; z++)
	{
		for (int y=0; y<m_Dims.y; y++)
		{
			for (int x=0; x<m_Dims.x; x++)
			{
				CalculateProbe(Vec3i(x, y, z));
			}
		}
	}

	Save("lightprobes.probe");
}

void LightProbes::CalculateProbe(const Vec3i &pt)
{
	Vector3 pos;
	GetProbePosition(pt, pos);

	LightProbe *pProbe = GetProbe(pt);
	assert (pProbe);

	for (int l=0; l<m_pLightMapper->m_Lights.size(); l++)
	{
		const LightMapper_Base::LLight &light = m_pLightMapper->m_Lights[l];

		// ray from probe to light
		Ray r;
		r.o = pos;
		Vector3 offset = light.pos - pos;
		float dist_to_light = offset.length();

		// NYI
		if (dist_to_light == 0.0f)
			continue;

		r.d = offset / dist_to_light;

		// to start with just trace a ray to each light
		Vector3 bary;
		float nearest_t;
		int tri_id = m_pLightMapper->m_Scene.FindIntersect_Ray(r, bary, nearest_t);

		// light blocked
		if ((tri_id != -1) && (nearest_t < dist_to_light))
			continue;

		// light got through..
		// pProbe->m_fLight += 1.0f;
		LightProbe::Contribution * pCont = pProbe->m_Contributions.request();
		pCont->light_id = l;

		// calculate power based on distance
		pCont->power = CalculatePower(light, dist_to_light, pos);

		print_line("\tprobe " + pt.ToString() + " light " + itos (l) + " power " + String(Variant(pCont->power)));
	}

}

float LightProbes::CalculatePower(const LightMapper_Base::LLight &light, float dist, const Vector3 &pos) const
{
	// ignore light energy for getting more detail, we can apply the light energy at runtime
	float power = m_pLightMapper->LightDistanceDropoff(dist, light, 1.0f);
	return power;
}


void LightProbes::Save_Vector3(FileAccess *f, const Vector3 &v)
{
	f->store_float(v.x);
	f->store_float(v.y);
	f->store_float(v.z);
}

void LightProbes::Save(const char * pszFilename)
{
	// try to open file for writing
	FileAccess *f = FileAccess::open(pszFilename, FileAccess::WRITE);
	ERR_FAIL_COND_MSG(!f, "Cannot create file at path '" + String(pszFilename) + "'.");

	// fourcc
	f->store_string("Prob");

	// first save dimensions
	f->store_16(m_Dims.x);
	f->store_16(m_Dims.y);
	f->store_16(m_Dims.z);

	// minimum pos
	Save_Vector3(f, m_ptMin);

	// voxel size
	Save_Vector3(f, m_VoxelSize);

	// num lights
	f->store_16(m_pLightMapper->m_Lights.size());

	for (int n=0; n<m_pLightMapper->m_Lights.size(); n++)
	{
		// light info
		const LightMapper_Base::LLight &light = m_pLightMapper->m_Lights[n];

		// type
		f->store_8(light.type);

		// pos, direction
		Save_Vector3(f, light.pos);
		Save_Vector3(f, light.dir);

		// energy
		f->store_float(light.energy);
		f->store_float(light.range);

		// color
		f->store_float(light.color.r);
		f->store_float(light.color.g);
		f->store_float(light.color.b);

		f->store_float(light.spot_angle_radians);
	}


	// now save voxels .. 2 passes to make easier to compress with zip.
	// the first pass will be easily compressable
	for (int z=0; z<m_Dims.z; z++)
	{
		for (int y=0; y<m_Dims.y; y++)
		{
			for (int x=0; x<m_Dims.x; x++)
			{
				LightProbe * probe = GetProbe(Vec3i(x, y, z));
				probe->Save(f);
			}
		}
	}

	for (int z=0; z<m_Dims.z; z++)
	{
		for (int y=0; y<m_Dims.y; y++)
		{
			for (int x=0; x<m_Dims.x; x++)
			{
				LightProbe * probe = GetProbe(Vec3i(x, y, z));
				probe->Save_Secondary(f);
			}
		}
	}


	if (f) {
		f->close();
		memdelete(f);
	}
}

void LightProbe::Save_Secondary(FileAccess *f)
{
	int nContribs = m_Contributions.size();
	if (nContribs > 255) nContribs = 255;

	for (int n=0; n<nContribs; n++)
	{
		const Contribution &c = m_Contributions[n];

		//f->store_8(c.light_id);

		// convert power to linear 0 to 255
		float p = c.power;
		p *= 0.5f;
		p *= 255.0f;

		int i = p;
		i = CLAMP(i, 0, 255);
		f->store_8(i);
	}
}

void LightProbe::Save(FileAccess *f)
{
	// number of contributions
	int nContribs = m_Contributions.size();
	if (nContribs > 255) nContribs = 255;

	// store number of contribs
	f->store_8(nContribs);

	for (int n=0; n<nContribs; n++)
	{
		const Contribution &c = m_Contributions[n];

		f->store_8(c.light_id);

		// convert power to linear 0 to 255
//		float p = c.power;
//		p *= 0.5f;
//		p *= 255.0f;

//		int i = p;
//		i = CLAMP(i, 0, 255);
//		f->store_8(i);
	}
}


} // namespace
