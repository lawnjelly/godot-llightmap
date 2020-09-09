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

		// todo .. we only need to store one set of UVs
		// (the orig face and the merged face store UV and UV2 respectively)
		Vector2 m_UVs[3];

		// when standardized, a face has an offset
		// to get back to the original format
		int offset;

		bool Vector3Less(const Vector3 &a, const Vector3 &b) const
		{
			if (a.x < b.x) return true;
			if (a.y < b.y) return true;
			if (a.z < b.z) return true;
			return false;
		}

		template <class T> void Shift(int shift, T * p)
		{
			switch (shift)
			{
			case 1:
				{
					T temp = p[0];
					p[0] = p[1];
					p[1] = p[2];
					p[2] = temp;
				}
				break;
			case 2:
				{
					T temp = p[0];
					p[0] = p[2];
					p[2] = p[1];
					p[1] = temp;
				}
				break;
			}
		}

		// to make faces comparable
		int StandardizeOrder()
		{
			int first = 0;
			Vector3 biggest = m_Pos[0];
			if (Vector3Less(biggest, m_Pos[1]))
			{
				biggest = m_Pos[1];
				first = 1;
			}
			if (Vector3Less(biggest, m_Pos[2]))
			{
				biggest = m_Pos[2];
				first = 2;
			}

			// resync
			Shift(first, m_Pos);
			Shift(first, m_Norm);
			Shift(first, m_index);
			Shift(first, m_UVs);

			// save the standardize offset
			offset = first;

			return first;
		}

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

	// these are now not crucial because the algorithm finds the best fit.
	// These are only used for quick reject.
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

	// returns goodness of fit
	float DoFacesMatch(const LFace &a, const LFace &b) const;
	void AddMergedTriangle(const LFace &forig, const LFace &forig_modelspace, const LFace &fmerged, LVector<LVert> &UniqueVerts, PoolVector<int> &UniqueIndices);

	void Transform_Verts(const PoolVector<Vector3> &ptsLocal, PoolVector<Vector3> &ptsWorld, const Transform &tr) const;
	void Transform_Norms(const PoolVector<Vector3> &normsLocal, PoolVector<Vector3> &normsWorld, const Transform &tr) const;
};


} // namespace
