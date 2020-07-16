#include "llighttests_simd.h"

namespace LM {



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

const __m128 oneM128 = _mm_set1_ps(1.0f);
const __m128 minusOneM128 = _mm_set1_ps(-1.0f);
const __m128 positiveEpsilonM128 = _mm_set1_ps(0.000001f);
const __m128 negativeEpsilonM128 = _mm_set1_ps(-0.000001f);
const __m128 zeroM128 = _mm_set1_ps(0.0f);


void PackedRay::Create(const Ray &ray)
{
	//m_OrigRay = ray;

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
//		// find normals manually
//		// there is at least one winner, so no need to set this
//		int winner = 0;

//		float * pT = (float*) &t;
//		for (int n=0; n<4; n++)
//		{
//			// calculate normal
//			Tri test_tri;
//			packedTris.ExtractTriangle(n, test_tri);
//			Vector3 norm;
//			test_tri.FindNormal_EdgeForm(norm);

//			float dot = m_OrigRay.d.dot(norm);
//			m_OrigRay.d.dot(norm);

//		} // for n

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

	/*
	__m128 tResults = _mm256_blendv_ps(t, minusOneM256, failed);

	int mask = _mm256_movemask_ps(tResults);
	if (mask != 0xFF)
	{
		// There is at least one intersection
		result.idx = -1;

		float* ptr = (float*)&tResults;
		for (int i = 0; i < 8; ++i)
		{
			if (ptr[i] >= 0.0f && ptr[i] < result.t)
			{
				result.t = ptr[i];
				result.idx = i;
			}
		}

		return result.idx != -1;
	}
	*/

	//    return false;
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
//	return pray.Intersect(ptris, r_nearest_t, r_winner);
}


bool LightTests_SIMD::TestIntersect4(const Tri *tris[4], const Ray &ray, float &r_nearest_t, int &r_winner) const
{
	/*
	PackedTriangles ptris;

	// load up .. this could cause all kinds of cache misses..
	for (int m=0; m<3; m++)
	{
		// NOTE : set ps takes input IN REVERSE!!!!
		ptris.e1[m].mm128 = _mm_set_ps(tris[3]->pos[0].coord[m],
		tris[2]->pos[0].coord[m],
		tris[1]->pos[0].coord[m],
		tris[0]->pos[0].coord[m]);

		ptris.e2[m].mm128 = _mm_set_ps(tris[3]->pos[1].coord[m],
		tris[2]->pos[1].coord[m],
		tris[1]->pos[1].coord[m],
		tris[0]->pos[1].coord[m]);

		ptris.v0[m].mm128 = _mm_set_ps(tris[3]->pos[2].coord[m],
		tris[2]->pos[2].coord[m],
		tris[1]->pos[2].coord[m],
		tris[0]->pos[2].coord[m]);
	}
//	ptris.inactiveMask = _mm_set1_ps(0.0f);

	// simd ray
	PackedRay pray;
	pray.m_origin[0] = _mm_set1_ps(ray.o.x);
	pray.m_origin[1] = _mm_set1_ps(ray.o.y);
	pray.m_origin[2] = _mm_set1_ps(ray.o.z);

	pray.m_direction[0] = _mm_set1_ps(ray.d.x);
	pray.m_direction[1] = _mm_set1_ps(ray.d.y);
	pray.m_direction[2] = _mm_set1_ps(ray.d.z);

	return pray.intersect(ptris, r_nearest_t, r_winner);
	*/
	return true;
}

} // namespace


/*
#define LLIGHT_SIMD_CASE(INDEX) {float check_dist = pT[INDEX];\
if (check_dist < nearest_dist)\
{\
nearest_dist = check_dist;\
winner = INDEX+1;\
}\
}

		mask = (mask ^ 15) & 15;

		// use the mask to switch
		// reverse order?
		switch (mask)
		{
		default:
			{
				// should never happen
				winner = 0;
			}
			break;
		case 1:
			{
				LLIGHT_SIMD_CASE(0)
			}
			break;
		case 2:
			{
				LLIGHT_SIMD_CASE(1)
			}
			break;
		case 3:
			{
				LLIGHT_SIMD_CASE(0)
				LLIGHT_SIMD_CASE(1)
			}
			break;
		case 4:
			{
				LLIGHT_SIMD_CASE(2)
			}
			break;
		case 5:
			{
				LLIGHT_SIMD_CASE(2)
				LLIGHT_SIMD_CASE(0)
			}
			break;
		case 6:
			{
				LLIGHT_SIMD_CASE(2)
				LLIGHT_SIMD_CASE(1)
			}
			break;
		case 7:
			{
				LLIGHT_SIMD_CASE(2)
				LLIGHT_SIMD_CASE(1)
				LLIGHT_SIMD_CASE(0)
			}
			break;
		case 8:
			{
				LLIGHT_SIMD_CASE(3)
			}
			break;
		case 9:
			{
				LLIGHT_SIMD_CASE(3)
				LLIGHT_SIMD_CASE(0)
			}
			break;
		case 10:
			{
				LLIGHT_SIMD_CASE(3)
				LLIGHT_SIMD_CASE(1)
			}
			break;
		case 11:
			{
				LLIGHT_SIMD_CASE(3)
				LLIGHT_SIMD_CASE(1)
				LLIGHT_SIMD_CASE(0)
			}
			break;
		case 12:
			{
				LLIGHT_SIMD_CASE(3)
				LLIGHT_SIMD_CASE(2)
			}
			break;
		case 13:
			{
				LLIGHT_SIMD_CASE(3)
				LLIGHT_SIMD_CASE(2)
				LLIGHT_SIMD_CASE(0)
			}
			break;
		case 14:
			{
				LLIGHT_SIMD_CASE(3)
				LLIGHT_SIMD_CASE(2)
				LLIGHT_SIMD_CASE(1)
			}
			break;
		case 15:
			{
				LLIGHT_SIMD_CASE(3)
				LLIGHT_SIMD_CASE(2)
				LLIGHT_SIMD_CASE(1)
				LLIGHT_SIMD_CASE(0)
			}
			break;
		}

#undef LLIGHT_SIMD_CASE
*/
