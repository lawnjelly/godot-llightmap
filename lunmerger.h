#pragma once

#include "scene/3d/mesh_instance.h"
#include "lvector.h"

namespace LM {

class UnMerger
{
	struct LFace
	{
		Vector3 m_Pos[3];
		Vector3 m_Norm[3];
		int m_index[3];

		String ToString() const;
	};

	// unique vert
	struct LVert
	{
		Vector3 m_Pos;
		Vector3 m_Norm;
		Vector2 m_UV;
		Vector2 m_UV2;

		bool ApproxEqual(const LVert &o) const;
	};

	// the merged mesh data to be passed for unmerging meshes
	struct LMerged
	{
		PoolVector<Vector3> m_Verts;
		PoolVector<Vector3> m_Norms;
		PoolVector<Vector2> m_UV2s;
		PoolVector<int> m_Inds;

		int m_nFaces;

		// precreate LFaces
		LVector<LFace> m_LFaces;
	};


	// internal data
	// merging params
	struct LMergeParams
	{
		float m_fThresholdDist;
		float m_fThresholdDist_Squared;
		float m_fThresholdDot;
	} m_MergeParams;

public:
	UnMerger();
	bool UnMerge(const Merger &merger, const MeshInstance &source_mi);
	void SetUnMergeParams(float thresh_dist, float thresh_dot);


private:
	bool FillMergedFromMesh(LMerged &merged, const MeshInstance &mesh);
	bool UnMerge_Mesh(MeshInstance &mi, LMerged &merged);

	int FindOrAddVert(LVector<LVert> &uni_verts, const LVert &vert) const;
	int DoFacesMatch(const LFace& sob_f, const LFace &m_face) const;
	int DoFacesMatch_Offset(const LFace& sob_f, const LFace &m_face, int offset) const;
	bool DoFaceVertsApproxMatch(const LFace& sob_f, const LFace &m_face, int c0, int c1, bool bDebug) const;
	bool DoPosNormsApproxMatch(const Vector3 &a_pos, const Vector3 &a_norm, const Vector3 &b_pos, const Vector3 &b_norm, bool bDebug) const;

	void Transform_Verts(const PoolVector<Vector3> &ptsLocal, PoolVector<Vector3> &ptsWorld, const Transform &tr) const;
	void Transform_Norms(const PoolVector<Vector3> &normsLocal, PoolVector<Vector3> &normsWorld, const Transform &tr) const;
};


} // namespace
