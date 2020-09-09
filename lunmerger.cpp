#include "lunmerger.h"
#include "lmerger.h"

namespace LM {


UnMerger::UnMerger()
{
	SetUnMergeParams(0.1f, 0.95f);
}

void UnMerger::SetUnMergeParams(float thresh_dist, float thresh_dot)
{
	m_MergeParams.m_fThresholdDist = thresh_dist;
	m_MergeParams.m_fThresholdDist_Squared = thresh_dist * thresh_dist;
	m_MergeParams.m_fThresholdDot = thresh_dot;
}


String UnMerger::LFace::ToString() const
{
	String sz;

	for (int c=0; c<3; c++)
	{
		sz += String(Variant(m_Pos[c]));
		sz += ", ";
	}

	sz += "norm : ";
	for (int c=0; c<3; c++)
	{
		sz += String(Variant(m_Norm[c]));
		sz += ", ";
	}

	return sz;
}


bool UnMerger::LVert::ApproxEqual(const LVert &o) const
{
	if (m_Pos != o.m_Pos)
		return false;
	if (m_Norm != o.m_Norm)
		return false;
	if (m_UV != o.m_UV)
		return false;
	if (m_UV2 != o.m_UV2)
		return false;

	return true;
}

bool UnMerger::UnMerge(const Merger &merger, const MeshInstance &source_mi)
{
	LMerged merged;
	if (!FillMergedFromMesh(merged, source_mi))
		return false;

	for (int n=0; n<merger.m_Meshes.size(); n++)
	{
		if (!UnMerge_Mesh(*merger.m_Meshes[n], merged))
			return false;
	}

	return true;
}

bool UnMerger::UnMerge_Mesh(MeshInstance &mi, LMerged &merged)
{
	//LPRINT(2, "UnMerge_SOB " + mi.get_name());

	Ref<Mesh> rmesh = mi.get_mesh();

	if (rmesh->get_surface_count() == 0)
	{
		WARN_PRINT_ONCE("UnMerger::UnMerge_Mesh mesh has no surfaces");
		return false;
	}

	Array arrays = rmesh->surface_get_arrays(0);
	PoolVector<Vector3> verts = arrays[VS::ARRAY_VERTEX];
	PoolVector<Vector3> norms = arrays[VS::ARRAY_NORMAL];
	PoolVector<Vector2> uv1s = arrays[VS::ARRAY_TEX_UV];
	PoolVector<int> inds = arrays[VS::ARRAY_INDEX];

	// if the mesh contains no indices, create some
	EnsureIndicesValid(inds, verts);

	// we need to get the vert positions and normals from local space to world space to match up with the
	// world space coords in the merged mesh
	PoolVector<Vector3> world_verts;
	PoolVector<Vector3> world_norms;
	Transform trans = mi.get_global_transform();
	Transform_Verts(verts, world_verts, trans);
	Transform_Norms(norms, world_norms, trans);


	// these are the uvs to be filled in the sob
	PoolVector<Vector2> uv2s;
	uv2s.resize(verts.size());

	// for each face in the SOB, attempt to find matching face in the merged mesh
	int nFaces = inds.size() / 3;
	int iCount = 0;

	int nMergedFaces = merged.m_nFaces;


	iCount = 0;

	// the number of unique verts in the UV2 mapped mesh may be HIGHER
	// than the original mesh, because verts with same pos / norm / uv may now have
	// different UV2. So we will be recreating the entire
	// mesh data with a new set of unique verts
	LVector<LVert> UniqueVerts;
	PoolVector<int> UniqueIndices;

	for (int f=0; f<nFaces; f++)
	{
		// create the original face in world space and standardize order
		LFace lf;
		// create a copy original face with model space coords
		LFace lf_modelspace;

		for (int c=0; c<3; c++)
		{
			int ind = inds[iCount++];
			lf.m_Pos[c] = world_verts[ind];
			lf.m_Norm[c] = world_norms[ind];
			lf.m_UVs[c] = uv1s[ind];

			lf.m_index[c] = ind;

			// model space .. needed for final reconstruction
			lf_modelspace.m_Pos[c] = verts[ind];
			lf_modelspace.m_Norm[c] = norms[ind];
		}

		lf.StandardizeOrder();
		// DONT STANDARDIZE THE MODEL SPACE VERSION, it may need
		// to be rejigged when finally adding the final face

#ifdef LDEBUG_UNMERGE
		LPRINT(5, "lface : " + lf.ToString());
#endif

		// find matching face
		float best_goodness_of_fit = FLT_MAX;
		int best_matching_merged_face = -1;

		for (int mf=0; mf<nMergedFaces; mf++)
		{
			// construct merged lface
			// no need to standardize, this is done as a pre-process now.
			LFace mlf = merged.m_LFaces[mf];

			// determine how good a match these 2 faces are
			float goodness_of_fit = DoFacesMatch(lf, mlf);

			// is it a better fit than the best match so far?
			if (goodness_of_fit < best_goodness_of_fit)
			{
				best_goodness_of_fit = goodness_of_fit;
				best_matching_merged_face = mf;

				// special case .. perfect goodness of fit
				if (goodness_of_fit == 0.0f)
					break;
			}

		} // for through merged faces

		// special case
		if (best_matching_merged_face == -1)
		{
			String sz = "\tface " + itos(f);
			sz += " no match";
			print_line (sz);
		}
		else
		{
			// standardized order is preprocessed
			LFace mlf = merged.m_LFaces[best_matching_merged_face];

			AddMergedTriangle(lf, lf_modelspace, mlf, UniqueVerts, UniqueIndices);
		}

	} // for through original faces

	print_line("UnMerge MI : " + mi.get_name());// + "\tFirstVert : " + itos(vert_count) + "\tNumUVs : " + itos(verts.size()));


	// rebuild the sob mesh, but now using the new uv2s

	//	LPRINT(2, "\t\tOrig Verts: " + itos(verts.size()) + ", Unique Verts: " + itos(UniqueVerts.size()) + ", Num Tris: " + itos(inds.size() / 3));

	// construct unique pool vectors to pass to mesh construction
	PoolVector<Vector3> unique_poss;
	PoolVector<Vector3> unique_norms;
	PoolVector<Vector2> unique_uv1s;
	PoolVector<Vector2> unique_uv2s;

	for (int n=0; n<UniqueVerts.size(); n++)
	{
		const LVert &v = UniqueVerts[n];
		unique_poss.push_back(v.m_Pos);
		unique_norms.push_back(v.m_Norm);
		unique_uv1s.push_back(v.m_UV);
		unique_uv2s.push_back(v.m_UV2);
	}



	Ref<ArrayMesh> am_new;
	am_new.instance();
	Array arr;
	arr.resize(Mesh::ARRAY_MAX);
	arr[Mesh::ARRAY_VERTEX] = unique_poss;
	arr[Mesh::ARRAY_NORMAL] = unique_norms;
	arr[Mesh::ARRAY_INDEX] = UniqueIndices;
	arr[Mesh::ARRAY_TEX_UV] = unique_uv1s;
	arr[Mesh::ARRAY_TEX_UV2] = unique_uv2s;

	am_new->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arr);

	//	am_new->surface_set_material(0, mat);

	// get the old material before setting the new mesh
	Ref<Material> mat = mi.get_surface_material(0);
	bool mat_on_mi = true;
	if (!mat.ptr())
	{
		mat_on_mi = false;
		mat = rmesh->surface_get_material(0);
	}


	// hopefully the old mesh will be ref count freed? ???
	mi.set_mesh(am_new);

	// mesh instance has the material?
	//	Ref<Material> mat = mi.get_surface_material(0);
	//	if (mat_on_mi)
	//	{
	if (!mat.ptr())
	{
		print_line("\t\tno material found");
	}

	mi.set_surface_material(0, mat);
	//	}
	//	else
	//	{
	// mesh has the material?
	//	Ref<Material> smat = rmesh->surface_get_material(0);
	//mi.set_surface_material(0, smat);
	//	}

	return true;
}

// we now have all the info needed in the orig face and the merged face
void UnMerger::AddMergedTriangle(const LFace &forig, const LFace &forig_modelspace, const LFace &fmerged, LVector<LVert> &UniqueVerts, PoolVector<int> &UniqueIndices)
{
	for (int c=0; c<3; c++)
	{
		// the model space tuple may be offset
		int t = c + (forig.offset);
		t %= 3;

		// construct the unique vert
		LVert uvert;
		//int orig_ind = forig.m_index[c];
		uvert.m_Pos = forig_modelspace.m_Pos[t];
		uvert.m_Norm = forig_modelspace.m_Norm[t];

		// first UV from the original
		uvert.m_UV = forig.m_UVs[c];

		// 2nd UV from the merged
		uvert.m_UV2 = fmerged.m_UVs[c];

		// find it or add to the list
		int ind_uni_vert = FindOrAddVert(UniqueVerts, uvert);

		// add the index to form the triangle list
		UniqueIndices.push_back(ind_uni_vert);
	}
}

bool UnMerger::FillMergedFromMesh(LMerged &merged, const MeshInstance &mesh)
{
	Ref<Mesh> rmesh = mesh.get_mesh();

	if (rmesh->get_surface_count() == 0)
	{
		WARN_PRINT("UnMerger::FillMergedFromMesh no surfaces");
		return false;
	}

	Array arrays = rmesh->surface_get_arrays(0);

	if (arrays.size() <= 1)
	{
		WARN_PRINT("UnMerger::FillMergedFromMesh no arrays");
		return false;
	}

	merged.m_Verts = arrays[VS::ARRAY_VERTEX];
	merged.m_Norms = arrays[VS::ARRAY_NORMAL];
	merged.m_UV2s = arrays[VS::ARRAY_TEX_UV2];
	merged.m_Inds = arrays[VS::ARRAY_INDEX];
	//	PoolVector<Vector2> p_UV1s = arrays[VS::ARRAY_TEX_UV];

	// special case, if no indices, create some
	if (!merged.m_Inds.size())
	{
		WARN_PRINT("UnMerger::FillMergedFromMesh no indices");
	}
	//	if (!EnsureIndicesValid(merged.m_Inds, merged.m_Verts))
	//	{
	//		print_line("invalid triangles due to duplicate indices detected in " + mesh.get_name());
	//	}

	merged.m_nFaces = merged.m_Inds.size() / 3;

	if (merged.m_UV2s.size() == 0)
	{
		print_line( "Merged mesh has no secondary UVs, using primary UVs");

		merged.m_UV2s = arrays[VS::ARRAY_TEX_UV];

		if (merged.m_UV2s.size() == 0)
		{
			print_line( "Merged mesh has no UVs");
			return false;
		}
	}

	int miCount = 0;
	for (int mf=0; mf<merged.m_nFaces; mf++)
	{
		// construct merged lface
		LFace mlf;
		for (int c=0; c<3; c++)
		{
			int ind = merged.m_Inds[miCount++];
			mlf.m_Pos[c] = merged.m_Verts[ind];
			mlf.m_Norm[c] = merged.m_Norms[ind].normalized();
			mlf.m_index[c] = ind;

			mlf.m_UVs[c] = merged.m_UV2s[ind];
		}

		// we precalculate the order of the merged faces .. these will be gone through lots of times
		mlf.StandardizeOrder();

		merged.m_LFaces.push_back(mlf);
	}

	return true;
}

void UnMerger::Transform_Verts(const PoolVector<Vector3> &ptsLocal, PoolVector<Vector3> &ptsWorld, const Transform &tr) const
{
	for (int n=0; n<ptsLocal.size(); n++)
	{
		Vector3 ptWorld = tr.xform(ptsLocal[n]);
		ptsWorld.push_back(ptWorld);
	}
}

void UnMerger::Transform_Norms(const PoolVector<Vector3> &normsLocal, PoolVector<Vector3> &normsWorld, const Transform &tr) const
{
	for (int n=0; n<normsLocal.size(); n++)
	{
		// hacky way for normals, we should use transpose of inverse matrix, dunno if godot supports this
		Vector3 ptNormA = Vector3(0, 0, 0);
		Vector3 ptNormWorldA = tr.xform(ptNormA);
		Vector3 ptNormWorldB = tr.xform(normsLocal[n]);

		Vector3 ptNorm = ptNormWorldB - ptNormWorldA;

		ptNorm = ptNorm.normalized();

		normsWorld.push_back(ptNorm);
	}
}



// returns -1 if no match, or the offset 0 1 or 2
float UnMerger::DoFacesMatch(const LFace &a, const LFace &b) const
{
	//	const float thresh_diff = m_MergeParams.m_fThresholdDist;
	const float thresh_diff_squared = m_MergeParams.m_fThresholdDist_Squared;

	// now we can compare by tuple
	float goodness_of_fit = 0.0f;

	for (int n=0; n<3; n++)
	{
		// positions
		Vector3 pos_diff = b.m_Pos[n] - a.m_Pos[n];
		float sl = pos_diff.length_squared();

		// quick reject
		if (sl > thresh_diff_squared)
			return FLT_MAX;

		goodness_of_fit += sl;

		// normals
		// CHECKING THE NORMALS was resulting in some incorrect results. xatlas must change them.
		// But in fact now the pos checks are accurate, with vertex ordering, back faces will not be selected
		// so normals don't need checking.

		//		Vector3 norm_a = a.m_Norm[n];
		//		Vector3 norm_b = b.m_Norm[n];

		//		norm_a.normalize();
		//		norm_b.normalize();

		// quick reject
		//float dot = norm_a.dot(norm_b);
		//if (dot < m_MergeParams.m_fThresholdDot)
		//	return FLT_MAX;
	}

	// if we got to here they matched to some degree
	return goodness_of_fit;
}


int UnMerger::FindOrAddVert(LVector<LVert> &uni_verts, const LVert &vert) const
{
	for (int n=0; n<uni_verts.size(); n++)
	{
		if (uni_verts[n].ApproxEqual(vert))
			return n;
	}

	// not found .. add to list
	uni_verts.push_back(vert);

	return uni_verts.size() - 1;
}


} // namespace
