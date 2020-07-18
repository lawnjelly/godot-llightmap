#include "llighttests_simd.h"

namespace LM {

#ifdef LLIGHTMAPPER_USE_SIMD

const __m128 oneM128 = _mm_set1_ps(1.0f);
const __m128 minusOneM128 = _mm_set1_ps(-1.0f);
const __m128 positiveEpsilonM128 = _mm_set1_ps(0.000001f);
const __m128 negativeEpsilonM128 = _mm_set1_ps(-0.000001f);
const __m128 zeroM128 = _mm_set1_ps(0.0f);

//const u_m128 ref_oneM128 = {1.0f, 1.0f, 1.0f, 1.0f};
//const u_m128 ref_minusOneM128 = {-1.0f, -1.0f, -1.0f, -1.0f};
//const u_m128 ref_positiveEpsilonM128 = {0.000001f, 0.000001f, 0.000001f, 0.000001f};
//const u_m128 ref_negativeEpsilonM128 = {-0.000001f, -0.000001f, -0.000001f, -0.000001f};
//const u_m128 ref_zeroM128 = {0.0f, 0.0f, 0.0f, 0.0f};



void multi_cross(__m128 result[3], const __m128 a[3], const __m128 b[3])
{
	__m128 tmp;
	__m128 tmp2;

	tmp = _mm_mul_ps(a[1], b[2]);
	tmp2 = _mm_mul_ps(b[1], a[2]);
	result[0] = _mm_sub_ps(tmp, tmp2);

	tmp = _mm_mul_ps(a[2], b[0]);
	tmp2 = _mm_mul_ps(b[2], a[0]);
	result[1] = _mm_sub_ps(tmp, tmp2);

	tmp = _mm_mul_ps(a[0], b[1]);
	tmp2 = _mm_mul_ps(b[0], a[1]);
	result[2] = _mm_sub_ps(tmp, tmp2);
}

//real_t Vector3::dot(const Vector3 &p_b) const {

//	return x * p_b.x + y * p_b.y + z * p_b.z;
//}


__m128 multi_dot(const __m128 a[3], const __m128 b[3])
{
	__m128 tmp = _mm_mul_ps(a[0], b[0]);
	__m128 tmp2 = _mm_mul_ps(a[1], b[1]);
	__m128 tmp3 = _mm_mul_ps(a[2], b[2]);

	return _mm_add_ps(tmp, _mm_add_ps(tmp2, tmp3));
	//    return fmadd(a[2], b[2], fmadd(a[1], b[1], mul(a[0], b[0])));
}

void multi_sub(__m128 result[3], const __m128 a[3], const __m128 b[3])
{
	result[0] = _mm_sub_ps(a[0], b[0]);
	result[1] = _mm_sub_ps(a[1], b[1]);
	result[2] = _mm_sub_ps(a[2], b[2]);
}


void PackedRay::Create(const Ray &ray)
{
	m_origin[0].mm128 = _mm_set1_ps(ray.o.x);
	m_origin[1].mm128 = _mm_set1_ps(ray.o.y);
	m_origin[2].mm128 = _mm_set1_ps(ray.o.z);

	m_direction[0].mm128 = _mm_set1_ps(ray.d.x);
	m_direction[1].mm128 = _mm_set1_ps(ray.d.y);
	m_direction[2].mm128 = _mm_set1_ps(ray.d.z);
}


// returns true if a hit
bool PackedRay::IntersectTest_CullBackFaces(const PackedTriangles& packedTris, float max_dist) const
{
	//Begin calculating determinant - also used to calculate u parameter
	// P
	__m128 q[3];
	multi_cross(q, &m_direction[0].mm128, &packedTris.e2[0].mm128);

	//if determinant is near zero, ray lies in plane of triangle
	// det
	__m128 a = multi_dot(&packedTris.e1[0].mm128, q);

	// reject based on det being close to zero
	// NYI

	// inv_det
	__m128 f = _mm_div_ps(oneM128, a);

	// distance from v1 to ray origin
	// T
	__m128 s[3];
	multi_sub(s, &m_origin[0].mm128, &packedTris.v0[0].mm128);

	// Calculate u parameter and test bound
	__m128 u = _mm_mul_ps(f, multi_dot(s, q));

	// the intersection lies outside triangle
	// NYI

	// Prepare to test v parameter
	// Q
	__m128 r[3];
	multi_cross(r, s, &packedTris.e1[0].mm128);

	// calculate V parameter and test bound
	// v
	__m128 v = _mm_mul_ps(f, multi_dot(&m_direction[0].mm128, r));

	// intersection outside of triangles?
	// NYI

	// t
	__m128 t = _mm_mul_ps(f, multi_dot(&packedTris.e2[0].mm128, r));

	// if t > epsilon, hit

	/////////////////////////////////////////
	// back face culling.
	// calculate face normal (not normalized)
	__m128 face_normals[3];
	multi_cross(face_normals, &packedTris.e2[0].mm128, &packedTris.e1[0].mm128);
	// dot ray direction
	__m128 ray_dot_normal = multi_dot(&m_direction[0].mm128, face_normals);
	/////////////////////////////////////////


	// Failure conditions
	// determinant close to zero?
	__m128 nohit = _mm_and_ps(_mm_cmpge_ps(a, negativeEpsilonM128), _mm_cmple_ps(a, positiveEpsilonM128));
	nohit = _mm_or_ps(nohit, _mm_cmple_ps(u, zeroM128));
	nohit = _mm_or_ps(nohit, _mm_cmple_ps(v, zeroM128));
	nohit = _mm_or_ps(nohit, _mm_cmpge_ps(_mm_add_ps(u, v), oneM128));
	nohit = _mm_or_ps(nohit, _mm_cmple_ps(t, zeroM128));
	//    failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_length));
	nohit = _mm_or_ps(nohit, packedTris.inactiveMask.mm128);

	// test against minimum t .. not sure if this will speed up or not.
	__m128 m_max_dist = _mm_set1_ps(max_dist);
	nohit = _mm_or_ps(nohit, _mm_cmpge_ps(t, m_max_dist));


	// backface culling - THE BIT OPERATION IS OPPOSITE.
	// We want to PREVENT failing if the backface is hit.
	__m128 backface_mask = _mm_cmpge_ps(ray_dot_normal, zeroM128);

	// always let through the bits set in backface mask
	nohit = _mm_or_ps(nohit, backface_mask);

	int mask = _mm_movemask_ps(nohit);
    if (mask != 15) // first 4 bits set
	{
		return true;
	} // mask

	return false;
}


// returns true if a hit
bool PackedRay::IntersectTest(const PackedTriangles& packedTris, float max_dist) const
{
	//Begin calculating determinant - also used to calculate u parameter
	// P
	__m128 q[3];
	multi_cross(q, &m_direction[0].mm128, &packedTris.e2[0].mm128);

	//if determinant is near zero, ray lies in plane of triangle
	// det
	__m128 a = multi_dot(&packedTris.e1[0].mm128, q);

	// reject based on det being close to zero
	// NYI

	// inv_det
	__m128 f = _mm_div_ps(oneM128, a);

	// distance from v1 to ray origin
	// T
	__m128 s[3];
	multi_sub(s, &m_origin[0].mm128, &packedTris.v0[0].mm128);

	// Calculate u parameter and test bound
	__m128 u = _mm_mul_ps(f, multi_dot(s, q));

	// the intersection lies outside triangle
	// NYI

	// Prepare to test v parameter
	// Q
	__m128 r[3];
	multi_cross(r, s, &packedTris.e1[0].mm128);

	// calculate V parameter and test bound
	// v
	__m128 v = _mm_mul_ps(f, multi_dot(&m_direction[0].mm128, r));

	// intersection outside of triangles?
	// NYI

	// t
	__m128 t = _mm_mul_ps(f, multi_dot(&packedTris.e2[0].mm128, r));

	// if t > epsilon, hit

	// Failure conditions
	// determinant close to zero?
	__m128 failed = _mm_and_ps(_mm_cmpge_ps(a, negativeEpsilonM128), _mm_cmple_ps(a, positiveEpsilonM128));
	failed = _mm_or_ps(failed, _mm_cmple_ps(u, zeroM128));
	failed = _mm_or_ps(failed, _mm_cmple_ps(v, zeroM128));

	failed = _mm_or_ps(failed, _mm_cmpge_ps(_mm_add_ps(u, v), oneM128));

	failed = _mm_or_ps(failed, _mm_cmple_ps(t, zeroM128));
	//    failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_length));

	failed = _mm_or_ps(failed, packedTris.inactiveMask.mm128);

	// test against minimum t .. not sure if this will speed up or not.
	__m128 m_max_dist = _mm_set1_ps(max_dist);
	failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_max_dist));

	//bool bHit  = false;
	int mask = _mm_movemask_ps(failed);
    if (mask != 15) // first 4 bits set
	{
		return true;
	} // mask

	return false;
}

/*
// returns winner index +1, or zero if no hit
int PackedRay::Intersect_TESTREF(const PackedTriangles& packedTris, float &nearest_dist) const
{
	//Begin calculating determinant - also used to calculate u parameter
	// P
	u_m128 q[3];
	u_m128::multi_cross(q, &m_direction[0], &packedTris.e2[0]);

	__m128 SIMD_q[3];
	multi_cross(SIMD_q, &m_direction[0].mm128, &packedTris.e2[0].mm128);


	//if determinant is near zero, ray lies in plane of triangle
	// det
	u_m128 a = u_m128::multi_dot(&packedTris.e1[0], q);


	__m128 SIMD_a = multi_dot(&packedTris.e1[0].mm128, SIMD_q);

	// reject based on det being close to zero
	// NYI

	// inv_det
	u_m128 f = u_m128::div_ps(ref_oneM128, a);

	__m128 SIMD_f = _mm_div_ps(oneM128, SIMD_a);


	// distance from v1 to ray origin
	// T
	u_m128 s[3];
	u_m128::multi_sub(s, &m_origin[0], &packedTris.v0[0]);


	__m128 SIMD_s[3];
	multi_sub(SIMD_s, &m_origin[0].mm128, &packedTris.v0[0].mm128);

	// Calculate u parameter and test bound
	u_m128 u;
	u_m128 dot_sq;
	dot_sq = u_m128::multi_dot(s, q);
	u = u_m128::mul_ps(f, dot_sq);


	__m128 SIMD_dot_sq = multi_dot(SIMD_s, SIMD_q);
	__m128 SIMD_u = _mm_mul_ps(SIMD_f, SIMD_dot_sq);

	// the intersection lies outside triangle
	// NYI

	// Prepare to test v parameter
	// Q
	u_m128 r[3];
	u_m128::multi_cross(r, s, &packedTris.e1[0]);

	__m128 SIMD_r[3];
	multi_cross(SIMD_r, SIMD_s, &packedTris.e1[0].mm128);


	// calculate V parameter and test bound
	// v
	u_m128 dot_dir_r;
	u_m128 v;
	dot_dir_r = u_m128::multi_dot(&m_direction[0], r);
	v = u_m128::mul_ps(f, dot_dir_r);


	__m128 SIMD_v = _mm_mul_ps(SIMD_f, multi_dot(&m_direction[0].mm128, SIMD_r));

	// intersection outside of triangles?
	// NYI

	// t
	u_m128 t, dot_edge_r;
	dot_edge_r = u_m128::multi_dot(&packedTris.e2[0], r);
	t = u_m128::mul_ps(f, dot_edge_r);


	__m128 SIMD_t = _mm_mul_ps(SIMD_f, multi_dot(&packedTris.e2[0].mm128, SIMD_r));

	// if t > epsilon, hit

	// Failure conditions
	// determinant close to zero?
	u_m128 failed = u_m128::mm_and(u_m128::ref_cmpge_ps(a, ref_negativeEpsilonM128), u_m128::ref_cmple_ps(a, ref_positiveEpsilonM128));

	__m128 SIMD_failed = _mm_and_ps(_mm_cmpge_ps(SIMD_a, negativeEpsilonM128), _mm_cmple_ps(SIMD_a, positiveEpsilonM128));



	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(u, ref_zeroM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(v, ref_zeroM128));
	//failed = u_m128::mm_or(failed, u_m128::ref_cmpge_ps(v, ref_zeroM128));

	SIMD_failed = _mm_or_ps(SIMD_failed, _mm_cmple_ps(SIMD_u, zeroM128));
	SIMD_failed = _mm_or_ps(SIMD_failed, _mm_cmple_ps(SIMD_v, zeroM128));


	failed = u_m128::mm_or(failed, u_m128::ref_cmpge_ps(u_m128::add_ps(u, v), ref_oneM128));
	SIMD_failed = _mm_or_ps(SIMD_failed, _mm_cmpge_ps(_mm_add_ps(SIMD_u, SIMD_v), oneM128));



	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(t, ref_zeroM128));
	//    failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_length));
	failed = u_m128::mm_or(failed, packedTris.inactiveMask);

	// test against minimum t .. not sure if this will speed up or not.
//	__m128 m_nearest_dist = _mm_set1_ps(nearest_dist);
//	failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_nearest_dist));

	//bool bHit  = false;
	int mask = failed.movemask_ps();
    if (mask != 15) // first 4 bits set
	{
		// there is at least one winner, so no need to set this
		int winner = 0;

		float * pT = (float*) &t;

		// use simplified way of getting results
		uint32_t * pResults = (uint32_t*) &failed;

		for (int n=0; n<4; n++)
		{
			// hit!
	//		if (!pResults[n] && n<packedTris.num_tris)
			if (!pResults[n])
			{
				if (pT[n] < nearest_dist)
				{
					nearest_dist = pT[n];
					//r_winner = n;
					winner = n+1;
					//bHit = true;
				}
			}
		} // for n

		return winner;
	} // mask

	return 0;
}
*/

// returns winner index +1, or zero if no hit
int PackedRay::Intersect(const PackedTriangles& packedTris, float &nearest_dist) const
{
	//Begin calculating determinant - also used to calculate u parameter
	// P
	__m128 q[3];
	multi_cross(q, &m_direction[0].mm128, &packedTris.e2[0].mm128);

	//if determinant is near zero, ray lies in plane of triangle
	// det
	__m128 a = multi_dot(&packedTris.e1[0].mm128, q);

	// reject based on det being close to zero
	// NYI

	// inv_det
	__m128 f = _mm_div_ps(oneM128, a);

	// distance from v1 to ray origin
	// T
	__m128 s[3];
	multi_sub(s, &m_origin[0].mm128, &packedTris.v0[0].mm128);

	// Calculate u parameter and test bound
	__m128 u = _mm_mul_ps(f, multi_dot(s, q));

	// the intersection lies outside triangle
	// NYI

	// Prepare to test v parameter
	// Q
	__m128 r[3];
	multi_cross(r, s, &packedTris.e1[0].mm128);

	// calculate V parameter and test bound
	// v
	__m128 v = _mm_mul_ps(f, multi_dot(&m_direction[0].mm128, r));

	// intersection outside of triangles?
	// NYI

	// t
	__m128 t = _mm_mul_ps(f, multi_dot(&packedTris.e2[0].mm128, r));

	// if t > epsilon, hit

	// Failure conditions
	// determinant close to zero?
	__m128 failed = _mm_and_ps(_mm_cmpge_ps(a, negativeEpsilonM128), _mm_cmple_ps(a, positiveEpsilonM128));

	//    __m128 failed = _mm_and_ps(
	//        _mm_cmp(a, negativeEpsilonM256, _CMP_GT_OQ),
	//        cmp(a, positiveEpsilonM256, _CMP_LT_OQ)
	//    );

	failed = _mm_or_ps(failed, _mm_cmple_ps(u, zeroM128));
	failed = _mm_or_ps(failed, _mm_cmple_ps(v, zeroM128));

	failed = _mm_or_ps(failed, _mm_cmpge_ps(_mm_add_ps(u, v), oneM128));

	failed = _mm_or_ps(failed, _mm_cmple_ps(t, zeroM128));
	//    failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_length));
	failed = _mm_or_ps(failed, packedTris.inactiveMask.mm128);

	// test against minimum t .. not sure if this will speed up or not.
//	__m128 m_nearest_dist = _mm_set1_ps(nearest_dist);
//	failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_nearest_dist));

	//bool bHit  = false;
	int mask = _mm_movemask_ps(failed);
    if (mask != 15) // first 4 bits set
	{
		// there is at least one winner, so no need to set this
		int winner = 0;

		float * pT = (float*) &t;

		// use simplified way of getting results
		uint32_t * pResults = (uint32_t*) &failed;

		for (int n=0; n<4; n++)
		{
			// hit!
	//		if (!pResults[n] && n<packedTris.num_tris)
			if (!pResults[n])
			{
				if (pT[n] < nearest_dist)
				{
					nearest_dist = pT[n];
					//r_winner = n;
					winner = n+1;
					//bHit = true;
				}
			}
		} // for n

		return winner;
	} // mask

	return 0;
}

bool LightTests_SIMD::TestIntersect4_Packed(const PackedTriangles &ptris, const Ray &ray, float &r_nearest_t, int &r_winner) const
{
	// simd ray
	PackedRay pray;
	pray.m_origin[0].mm128 = _mm_set1_ps(ray.o.x);
	pray.m_origin[1].mm128 = _mm_set1_ps(ray.o.y);
	pray.m_origin[2].mm128 = _mm_set1_ps(ray.o.z);

	pray.m_direction[0].mm128 = _mm_set1_ps(ray.d.x);
	pray.m_direction[1].mm128 = _mm_set1_ps(ray.d.y);
	pray.m_direction[2].mm128 = _mm_set1_ps(ray.d.z);

	r_winner = pray.Intersect(ptris, r_nearest_t);
	if (r_winner)
	{
		r_winner--;
		return true;
	}

	return false;
}


bool LightTests_SIMD::TestIntersect4(const Tri *tris[4], const Ray &ray, float &r_nearest_t, int &r_winner) const
{
	return true;
}

#else // not SIMD

// REFERENCE IMPLEMENTATION (for non SSE2 CPUs)

const u_m128 ref_oneM128 = {{1.0f, 1.0f, 1.0f, 1.0f}};
const u_m128 ref_minusOneM128 = {{-1.0f, -1.0f, -1.0f, -1.0f}};
const u_m128 ref_positiveEpsilonM128 = {{0.000001f, 0.000001f, 0.000001f, 0.000001f}};
const u_m128 ref_negativeEpsilonM128 = {{-0.000001f, -0.000001f, -0.000001f, -0.000001f}};
const u_m128 ref_zeroM128 = {{0.0f, 0.0f, 0.0f, 0.0f}};


void PackedRay::Create(const Ray &ray)
{
	m_origin[0] = u_m128::set1_ps(ray.o.x);
	m_origin[1] = u_m128::set1_ps(ray.o.y);
	m_origin[2] = u_m128::set1_ps(ray.o.z);

	m_direction[0] = u_m128::set1_ps(ray.d.x);
	m_direction[1] = u_m128::set1_ps(ray.d.y);
	m_direction[2] = u_m128::set1_ps(ray.d.z);
}


// returns true if a hit
bool PackedRay::IntersectTest_CullBackFaces(const PackedTriangles& packedTris, float max_dist) const
{
	//Begin calculating determinant - also used to calculate u parameter
	// P
	u_m128 q[3];
	u_m128::multi_cross(q, &m_direction[0], &packedTris.e2[0]);

	//if determinant is near zero, ray lies in plane of triangle
	// det
	u_m128 a = u_m128::multi_dot(&packedTris.e1[0], q);

	// reject based on det being close to zero
	// NYI

	// inv_det
	u_m128 f = u_m128::div_ps(ref_oneM128, a);

	// distance from v1 to ray origin
	// T
	u_m128 s[3];
	u_m128::multi_sub(s, &m_origin[0], &packedTris.v0[0]);

	// Calculate u parameter and test bound
	u_m128 u;
	u_m128 dot_sq;
	dot_sq = u_m128::multi_dot(s, q);
	u = u_m128::mul_ps(f, dot_sq);

	// the intersection lies outside triangle
	// NYI

	// Prepare to test v parameter
	// Q
	u_m128 r[3];
	u_m128::multi_cross(r, s, &packedTris.e1[0]);

	// calculate V parameter and test bound
	// v
	u_m128 dot_dir_r;
	u_m128 v;
	dot_dir_r = u_m128::multi_dot(&m_direction[0], r);
	v = u_m128::mul_ps(f, dot_dir_r);

	// intersection outside of triangles?
	// NYI

	// t
	u_m128 t, dot_edge_r;
	dot_edge_r = u_m128::multi_dot(&packedTris.e2[0], r);
	t = u_m128::mul_ps(f, dot_edge_r);

	// if t > epsilon, hit

	/////////////////////////////////////////
	// back face culling.
	// calculate face normal (not normalized)
	u_m128 face_normals[3];
	u_m128::multi_cross(face_normals, &packedTris.e2[0], &packedTris.e1[0]);
	// dot ray direction
	u_m128 ray_dot_normal = u_m128::multi_dot(&m_direction[0], face_normals);
	/////////////////////////////////////////


	// Failure conditions
	// determinant close to zero?
	u_m128 failed = u_m128::mm_and(u_m128::ref_cmpge_ps(a, ref_negativeEpsilonM128), u_m128::ref_cmple_ps(a, ref_positiveEpsilonM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(u, ref_zeroM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(v, ref_zeroM128));


	failed = u_m128::mm_or(failed, u_m128::ref_cmpge_ps(u_m128::add_ps(u, v), ref_oneM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(t, ref_zeroM128));
	//    failed = _mm_or_ps(failed, _mm_cmpge_ps(t, m_length));
	failed = u_m128::mm_or(failed, packedTris.inactiveMask);

	// test against max dist.
	u_m128 m_max_dist = u_m128::set1_ps(max_dist);
	failed = u_m128::mm_or(failed, u_m128::ref_cmpge_ps(t, m_max_dist));

	// backface culling - THE BIT OPERATION IS OPPOSITE.
	// We want to PREVENT failing if the backface is hit.
	u_m128 backface_mask = u_m128::ref_cmpge_ps(ray_dot_normal, ref_zeroM128);

	// always let through the bits set in backface mask
	failed = u_m128::mm_or(failed, backface_mask);

	//bool bHit  = false;
	int mask = failed.movemask_ps();
    if (mask != 15) // first 4 bits set
	{
		// there is at least one winner
		return true;
	}

	return false;
}

// returns true if a hit
bool PackedRay::IntersectTest(const PackedTriangles& packedTris, float max_dist) const
{
	//Begin calculating determinant - also used to calculate u parameter
	// P
	u_m128 q[3];
	u_m128::multi_cross(q, &m_direction[0], &packedTris.e2[0]);

	//if determinant is near zero, ray lies in plane of triangle
	// det
	u_m128 a = u_m128::multi_dot(&packedTris.e1[0], q);

	// reject based on det being close to zero
	// NYI

	// inv_det
	u_m128 f = u_m128::div_ps(ref_oneM128, a);

	// distance from v1 to ray origin
	// T
	u_m128 s[3];
	u_m128::multi_sub(s, &m_origin[0], &packedTris.v0[0]);

	// Calculate u parameter and test bound
	u_m128 u;
	u_m128 dot_sq;
	dot_sq = u_m128::multi_dot(s, q);
	u = u_m128::mul_ps(f, dot_sq);

	// the intersection lies outside triangle
	// NYI

	// Prepare to test v parameter
	// Q
	u_m128 r[3];
	u_m128::multi_cross(r, s, &packedTris.e1[0]);

	// calculate V parameter and test bound
	// v
	u_m128 dot_dir_r;
	u_m128 v;
	dot_dir_r = u_m128::multi_dot(&m_direction[0], r);
	v = u_m128::mul_ps(f, dot_dir_r);

	// intersection outside of triangles?
	// NYI

	// t
	u_m128 t, dot_edge_r;
	dot_edge_r = u_m128::multi_dot(&packedTris.e2[0], r);
	t = u_m128::mul_ps(f, dot_edge_r);

	// if t > epsilon, hit

	// Failure conditions
	// determinant close to zero?
	u_m128 failed = u_m128::mm_and(u_m128::ref_cmpge_ps(a, ref_negativeEpsilonM128), u_m128::ref_cmple_ps(a, ref_positiveEpsilonM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(u, ref_zeroM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(v, ref_zeroM128));

	failed = u_m128::mm_or(failed, u_m128::ref_cmpge_ps(u_m128::add_ps(u, v), ref_oneM128));

	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(t, ref_zeroM128));

	failed = u_m128::mm_or(failed, packedTris.inactiveMask);

	// test against max dist.
	u_m128 m_max_dist = u_m128::set1_ps(max_dist);
	failed = u_m128::mm_or(failed, u_m128::ref_cmpge_ps(t, m_max_dist));

	//bool bHit  = false;
	int mask = failed.movemask_ps();
    if (mask != 15) // first 4 bits set
	{
		// there is at least one winner
		return true;
	}

	return false;
}

// returns winner index +1, or zero if no hit
int PackedRay::Intersect(const PackedTriangles& packedTris, float &nearest_dist) const
{
	//Begin calculating determinant - also used to calculate u parameter
	// P
	u_m128 q[3];
	u_m128::multi_cross(q, &m_direction[0], &packedTris.e2[0]);

	//if determinant is near zero, ray lies in plane of triangle
	// det
	u_m128 a = u_m128::multi_dot(&packedTris.e1[0], q);

	// reject based on det being close to zero
	// NYI

	// inv_det
	u_m128 f = u_m128::div_ps(ref_oneM128, a);

	// distance from v1 to ray origin
	// T
	u_m128 s[3];
	u_m128::multi_sub(s, &m_origin[0], &packedTris.v0[0]);

	// Calculate u parameter and test bound
	u_m128 u;
	u_m128 dot_sq;
	dot_sq = u_m128::multi_dot(s, q);
	u = u_m128::mul_ps(f, dot_sq);

	// the intersection lies outside triangle
	// NYI

	// Prepare to test v parameter
	// Q
	u_m128 r[3];
	u_m128::multi_cross(r, s, &packedTris.e1[0]);

	// calculate V parameter and test bound
	// v
	u_m128 dot_dir_r;
	u_m128 v;
	dot_dir_r = u_m128::multi_dot(&m_direction[0], r);
	v = u_m128::mul_ps(f, dot_dir_r);

	// intersection outside of triangles?
	// NYI

	// t
	u_m128 t, dot_edge_r;
	dot_edge_r = u_m128::multi_dot(&packedTris.e2[0], r);
	t = u_m128::mul_ps(f, dot_edge_r);

	// if t > epsilon, hit

	// Failure conditions
	// determinant close to zero?
	u_m128 failed = u_m128::mm_and(u_m128::ref_cmpge_ps(a, ref_negativeEpsilonM128), u_m128::ref_cmple_ps(a, ref_positiveEpsilonM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(u, ref_zeroM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(v, ref_zeroM128));


	failed = u_m128::mm_or(failed, u_m128::ref_cmpge_ps(u_m128::add_ps(u, v), ref_oneM128));
	failed = u_m128::mm_or(failed, u_m128::ref_cmple_ps(t, ref_zeroM128));

	failed = u_m128::mm_or(failed, packedTris.inactiveMask);

	//bool bHit  = false;
	int mask = failed.movemask_ps();
    if (mask != 15) // first 4 bits set
	{
		// there is at least one winner, so no need to set this
		int winner = 0;

		float * pT = (float*) &t;

		// use simplified way of getting results
		uint32_t * pResults = (uint32_t*) &failed;

		for (int n=0; n<4; n++)
		{
			// hit!
	//		if (!pResults[n] && n<packedTris.num_tris)
			if (!pResults[n])
			{
				if (pT[n] < nearest_dist)
				{
					nearest_dist = pT[n];
					//r_winner = n;
					winner = n+1;
					//bHit = true;
				}
			}
		} // for n

		return winner;
	} // mask

	return 0;
}


#endif // REFERENCE IMPLEMENTATION

} // namespace



