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
		LFace lf;
		for (int c=0; c<3; c++)
		{
			int ind = inds[iCount++];
			lf.m_Pos[c] = world_verts[ind];
			lf.m_Norm[c] = world_norms[ind];

			lf.m_index[c] = ind;
		}

#ifdef LDEBUG_UNMERGE
		LPRINT(5, "lface : " + lf.ToString());
#endif

		// find matching face
//		int miCount = 0;
		bool bMatchFound = false;

		for (int mf=0; mf<nMergedFaces; mf++)
		{
			// construct merged lface
			const LFace &mlf = merged.m_LFaces[mf];
//			for (int c=0; c<3; c++)
//			{
//				int ind = merged.m_Inds[miCount++];
//				mlf.m_Pos[c] = merged.m_Verts[ind];
//				mlf.m_Norm[c] = merged.m_Norms[ind];

//				mlf.m_index[c] = ind;
//			}

			int match = DoFacesMatch(lf, mlf);
			if (match != -1)
			{
				// we found a match in the merged mesh! transfer the UV2s
				bMatchFound = true;

				//sz += " match found offset " + itos(match);

				// find the corresponding uv2s in the merged mesh face and add them (taking into account offset)
				Vector2 found_uvs[3];
				for (int c=0; c<3; c++)
				{
					found_uvs[c] = merged.m_UV2s[mlf.m_index[c]];
				}

				// add them
				for (int c=0; c<3; c++)
				{
					int which = (c + match)	% 3;

					// index for the uv2 should be the same as the vertex index in the 'to mesh' face
					int ind = lf.m_index[c];
					uv2s.set(ind, found_uvs[which]);

					{
						// construct the unique vert
						LVert uvert;
						int orig_ind = lf.m_index[c];
						uvert.m_Pos = verts[orig_ind];
						uvert.m_Norm = norms[orig_ind];
						uvert.m_UV = uv1s[orig_ind];
						uvert.m_UV2 = found_uvs[which];

						// find it or add to the list
						int ind_uni_vert = FindOrAddVert(UniqueVerts, uvert);

						// add the index to form the triangle list
						UniqueIndices.push_back(ind_uni_vert);
					}
				}

				break;
			}
		}

		// special case
		if (!bMatchFound)
		{
			// add some dummy uv2s
//			uv2s.push_back(Vector2(0, 0));
//			uv2s.push_back(Vector2(0, 0));
//			uv2s.push_back(Vector2(0, 0));

			String sz = "\tface " + itos(f);
			sz += " no match";
			print_line (sz);
		}

	}

//	LPRINT(2, "UnMerge MI : " + mi.get_name() + "\tFirstVert : " + itos(vert_count) + "\tNumUVs : " + itos(verts.size()));


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

	// hopefully the old mesh will be ref count freed? ???
	mi.set_mesh(am_new);


	// mesh instance has the material?
	Ref<Material> mat = mi.get_surface_material(0);
	if (mat.ptr())
	{
		mi.set_surface_material(0, mat);
	}
	else
	{
		// mesh has the material?
		Ref<Material> smat = rmesh->surface_get_material(0);
		mi.set_surface_material(0, smat);
	}

	return true;
}


bool UnMerger::FillMergedFromMesh(LMerged &merged, const MeshInstance &mesh)
{
	Ref<Mesh> rmesh = mesh.get_mesh();
	Array arrays = rmesh->surface_get_arrays(0);


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
		}
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

// -1 for no match, or 0 for 0 offset match, 1 for +1, 2 for +2 offset match...
int UnMerger::DoFacesMatch(const LFace& sob_f, const LFace &m_face) const
{
	// match one
	for (int offset = 0; offset < 3; offset++)
	{
		if (DoFaceVertsApproxMatch(sob_f, m_face, 0, offset, false))
		{
			int res = DoFacesMatch_Offset(sob_f, m_face, offset);

			if (res != -1)
				return res;
		}
	}

	return -1; // no match
}

int UnMerger::DoFacesMatch_Offset(const LFace& sob_f, const LFace &m_face, int offset) const
{
#ifdef LDEBUG_UNMERGE
	// debug
	String sz = "\t\tPOSS match sob : ";
	sz += sob_f.ToString();
	sz += "\n\t\t\tmerged : ";
	sz += m_face.ToString();
	LPRINT(2, sz);
#endif


	// does 2nd and third match?
	int offset1 = (offset + 1) % 3;
	if (!DoFaceVertsApproxMatch(sob_f, m_face, 1, offset1, true))
		return -1;

	int offset2 = (offset + 2) % 3;
	if (!DoFaceVertsApproxMatch(sob_f, m_face, 2, offset2, true))
		return -1;

	return offset;
}


bool UnMerger::DoFaceVertsApproxMatch(const LFace& sob_f, const LFace &m_face, int c0, int c1, bool bDebug) const
{
	return DoPosNormsApproxMatch(sob_f.m_Pos[c0], sob_f.m_Norm[c0], m_face.m_Pos[c1], m_face.m_Norm[c1], bDebug);
}

bool UnMerger::DoPosNormsApproxMatch(const Vector3 &a_pos, const Vector3 &a_norm, const Vector3 &b_pos, const Vector3 &b_norm, bool bDebug) const
{
	bDebug = false;

	float thresh_diff = m_MergeParams.m_fThresholdDist;
	float thresh_diff_squared = m_MergeParams.m_fThresholdDist_Squared;

	float x_diff = fabs (b_pos.x - a_pos.x);
	if (x_diff > thresh_diff)
	{
#ifdef LDEBUG_UNMERGE
		if (bDebug)
			LPRINT(5, "\t\t\t\tRejecting x_diff " + ftos(x_diff));
#endif
		return false;
	}

	float z_diff = fabs (b_pos.z - a_pos.z);
	if (z_diff > thresh_diff)
	{
#ifdef LDEBUG_UNMERGE
		if (bDebug)
			LPRINT(5, "\t\t\t\tRejecting z_diff " + ftos(z_diff));
#endif
		return false;
	}

	float y_diff = fabs (b_pos.y - a_pos.y);
	if (y_diff > thresh_diff)
	{
#ifdef LDEBUG_UNMERGE
		if (bDebug)
			LPRINT(5, "\t\t\t\tRejecting y_diff " + ftos(y_diff));
#endif
		return false;
	}



	Vector3 pos_diff = b_pos - a_pos;
	if (pos_diff.length_squared() > thresh_diff_squared) // 0.1
	{
#ifdef LDEBUG_UNMERGE
		if (bDebug)
			LPRINT(5, "\t\t\t\tRejecting length squared " + ftos(pos_diff.length_squared()));
#endif

		return false;
	}

	// make sure both are normalized
	Vector3 na = a_norm;//.normalized();
	Vector3 nb = b_norm;//.normalized();

	float norm_dot = na.dot(nb);
	if (norm_dot < m_MergeParams.m_fThresholdDot)
	{
#ifdef LDEBUG_UNMERGE
		if (bDebug)
			LPRINT(5, "\t\t\t\tRejecting normal " + ftos(norm_dot) + " na : " + String(na) + ", nb : " + String(nb));
#endif
		return false;
	}

	return true;
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
