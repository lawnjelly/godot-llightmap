#pragma once

#if defined(__x86_64__) || defined(_M_X64)
#define LLIGHTMAPPER_USE_SIMD
#endif


#ifdef LLIGHTMAPPER_USE_SIMD
#include <emmintrin.h>
#endif

#include "llighttypes.h"

namespace LM
{

union u_m128
{
#ifdef LLIGHTMAPPER_USE_SIMD
	__m128 mm128;
#endif
	struct
	{
		float v[4];
	};
	struct
	{
		uint32_t ui[4];
	};

	static u_m128 set1_ps(float f) {u_m128 r; r.v[0] = f; r.v[1] = f; r.v[2] = f; r.v[3] = f; return r;}

	static u_m128 mul_ps(const u_m128 &a, const u_m128 &b) {u_m128 r; for (int n=0; n<4; n++) {r.v[n] = a.v[n] * b.v[n];} return r;}
	static u_m128 div_ps(const u_m128 &a, const u_m128 &b) {u_m128 r; for (int n=0; n<4; n++) {r.v[n] = a.v[n] / b.v[n];} return r;}
	static u_m128 add_ps(const u_m128 &a, const u_m128 &b) {u_m128 r; for (int n=0; n<4; n++) {r.v[n] = a.v[n] + b.v[n];} return r;}
	static u_m128 sub_ps(const u_m128 &a, const u_m128 &b) {u_m128 r; for (int n=0; n<4; n++) {r.v[n] = a.v[n] - b.v[n];} return r;}

	static u_m128 ref_cmple_ps(const u_m128 &a, const u_m128 &b)
	{
		u_m128 res;
		for (int n=0; n<4; n++)
		{
			if (a.v[n] <= b.v[n])
				res.ui[n] = 0xFFFFFFFF;
			else
				res.ui[n] = 0;
		}
		return res;
	}
	static u_m128 ref_cmpge_ps(const u_m128 &a, const u_m128 &b)
	{
		u_m128 res;
		for (int n=0; n<4; n++)
		{
			if (a.v[n] >= b.v[n])
				res.ui[n] = 0xFFFFFFFF;
			else
				res.ui[n] = 0;
		}
		return res;
	}

	int movemask_ps() const
	{ // possibly needs this order reversing
		int m = 0;
		if (ui[0]) m |= 1;
		if (ui[1]) m |= 2;
		if (ui[2]) m |= 4;
		if (ui[3]) m |= 8;
		return m;
	}

	static u_m128 mm_or(const u_m128 &a, const u_m128 &b) {u_m128 r; for (int n=0; n<4; n++) {r.ui[n] = a.ui[n] | b.ui[n];} return r;}
	static u_m128 mm_and(const u_m128 &a, const u_m128 &b) {u_m128 r; for (int n=0; n<4; n++) {r.ui[n] = a.ui[n] & b.ui[n];} return r;}

	static void multi_sub(u_m128 result[3], const u_m128 a[3], const u_m128 b[3])
	{
		for (int n=0; n<3; n++)
			result[n] = sub_ps(a[n], b[n]);
	}

	static void multi_cross(u_m128 result[3], const u_m128 a[3], const u_m128 b[3])
	{
		u_m128 tmp;
		u_m128 tmp2;

		tmp = mul_ps(a[1], b[2]);
		tmp2 = mul_ps(b[1], a[2]);
		result[0] = sub_ps(tmp, tmp2);

		tmp = mul_ps(a[2], b[0]);
		tmp2 = mul_ps(b[2], a[0]);
		result[1] = sub_ps(tmp, tmp2);

		tmp = mul_ps(a[0], b[1]);
		tmp2 = mul_ps(b[0], a[1]);
		result[2] = sub_ps(tmp, tmp2);
	}

	static u_m128 multi_dot(const u_m128 a[3], const u_m128 b[3])
	{
		u_m128 r;
		u_m128 tmp, tmp2, tmp3, tmp4;

		tmp = mul_ps(a[0], b[0]);
		tmp2 = mul_ps(a[1], b[1]);
		tmp3 = mul_ps(a[2], b[2]);

		tmp4 = add_ps(tmp, tmp2);
		r = add_ps(tmp4, tmp3); // this may be dodgy, as argument is same as result
		return r;
	}
};



// 4 triangles in edge form
// i.e. edge1, edge2 and a vertex
struct PackedTriangles
{
	u_m128 e1[3];
	u_m128 e2[3];
	u_m128 v0[3];
	// 0 to 4 filled
//	uint32_t num_tris;
	u_m128 inactiveMask;

	void Create() {memset (this, 0, sizeof (PackedTriangles));} // num_tris = 0;
	void Set(int which_tri, const Tri &tri) {Set_e1(which_tri, tri); Set_e2(which_tri, tri); Set_v0(which_tri, tri);
//		num_tris = which_tri+1;
	}
	void ExtractTriangle(int w, Tri &tri) const // debug
	{
		for (int m=0; m<3; m++)
		{
			tri.pos[0].coord[m] = e1[m].v[w];
			tri.pos[1].coord[m] = e2[m].v[w];
			tri.pos[2].coord[m] = v0[m].v[w];
		}
	}
	void Finalize(int num_tris)
	{
		for (int n=num_tris; n<4; n++)
		{
			inactiveMask.v[n] = 1.0f;
		}
	}
	void MakeInactive(int which)
	{
		inactiveMask.v[which] = 1.0f;
	}

	void Set_e1(int which_tri, const Tri &tri) {for (int m=0; m<3; m++) {e1[m].v[which_tri] = tri.pos[0].coord[m];} }
	void Set_e2(int which_tri, const Tri &tri) {for (int m=0; m<3; m++) {e2[m].v[which_tri] = tri.pos[1].coord[m];} }
	void Set_v0(int which_tri, const Tri &tri) {for (int m=0; m<3; m++) {v0[m].v[which_tri] = tri.pos[2].coord[m];} }
};


struct PackedRay
{
	u_m128 m_origin[3];
	u_m128 m_direction[3];
	u_m128 m_length;

	void Create(const Ray &ray);
	int Intersect(const PackedTriangles& packedTris, float &nearest_dist) const;
//	int Intersect_ORIG(const PackedTriangles& packedTris, float &nearest_dist) const;
	bool IntersectTest(const PackedTriangles& packedTris, float max_dist) const;
	bool IntersectTest_CullBackFaces(const PackedTriangles& packedTris, float max_dist) const;
};


class LightTests_SIMD
{
public:
	bool TestIntersect4(const Tri *tris[4], const Ray &ray, float &r_nearest_t, int &r_winner) const;
	bool TestIntersect4_Packed(const PackedTriangles &ptris, const Ray &ray, float &r_nearest_t, int &r_winner) const;
};


} // namespace
