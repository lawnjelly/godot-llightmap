#pragma once

#include "llightimage.h"
#include "llighttypes.h"
#include "scene/3d/mesh_instance.h"

namespace LM {

class Stitcher {
public:
	struct UVSeam {
		Vector2 edge0[2];
		Vector2 edge1[2];
	};

	struct SeamEdge {
		Vector3 pos[2];
		Vector3 normal[2];
		Vector2 uv[2];

		_FORCE_INLINE_ bool operator<(const SeamEdge &p_edge) const {
			return pos[0].x < p_edge.pos[0].x;
		}
	};

	//	struct StitchMesh {
	//		MeshData data;
	//		int slice = 0;
	//		Vector2i offset;
	//		Vector2i size;
	//		bool cast_shadows;
	//		bool generate_lightmap;
	//		String node_name;
	//	};

	void StitchObjectSeams(const MeshInstance &p_mi, LightImage<FColor> &r_image, float distance_threshold, float normal_threshold, bool p_visualize_seams);

private:
	void _compute_seams(const Vector<Vector3> &points, const Vector<Vector2> &uv2s, const Vector<Vector3> &normals, Vector2i lm_size, LocalVector<UVSeam> &r_seams);
	void _fix_seams(const LocalVector<UVSeam> &p_seams, Vector3 *r_lightmap, Vector2i p_size);
	void _fix_seam(const Vector2 &p_pos0, const Vector2 &p_pos1, const Vector2 &p_uv0, const Vector2 &p_uv1, const Vector3 *p_read_buffer, Vector3 *r_write_buffer, const Vector2i &p_size);
};

float m_fDistanceThreshold;
float m_fNormalThreshold;
float m_fPositionEpsilon;
bool m_bVisualizeSeams;

} // namespace LM
