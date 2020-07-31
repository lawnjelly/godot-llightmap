#pragma once

#include "core/math/vector3.h"
#include "lvector.h"

namespace LM
{

struct MiniList
{
	uint32_t first;
	uint32_t num;
};


class Vec3i
{
public:
	Vec3i() {}
	Vec3i(int xx, int yy, int zz) {x = xx; y = yy; z = zz;}
	int32_t x, y, z;

	Vec3i &operator-=(const Vec3i &v) {x -= v.x; y -= v.y; z -= v.z; return *this;}
	int SquareLength() const {return (x * x) + (y * y)  + (z * z);}
	float Length() const {return sqrtf(SquareLength());}

	void Set(int xx, int yy, int zz) {x = xx; y = yy; z = zz;}
	void Set_Round(const Vector3 &p) {x = (int) (p.x + 0.5f); y = (int) (p.y + 0.5f); z = (int) (p.z + 0.5f);}
	void To(Vector3 &p) const {p.x = x; p.y = y; p.z = z;}
	String ToString() const {return itos(x) + ", " + itos(y) + ", " + itos(z);}
};

class Vec2_i16
{
public:
	Vec2_i16() {}
	Vec2_i16(int xx, int yy) {x = xx; y = yy;}
	int16_t x;
	int16_t y;
	void Set(int xx, int yy) {x = xx; y = yy;}
	bool IsZero() const {return (x == 0) && (y == 0);}
	bool IsNonZero() const {return !IsZero();}
};

class Tri
{
public:
	Vector3 pos[3];

	void FlipWinding() {Vector3 temp = pos[0]; pos[0] = pos[2]; pos[2] = temp;}
	void FindBarycentric(const Vector3 &pt, float &u, float &v, float &w) const;
	void InterpolateBarycentric(Vector3 &pt, const Vector3 &bary) const {InterpolateBarycentric(pt, bary.x, bary.y, bary.z);}
	void InterpolateBarycentric(Vector3 &pt, float u, float v, float w) const;
	void FindNormal_EdgeForm(Vector3 &norm) const
	{
		norm = -pos[0].cross(pos[1]);
		norm.normalize();
	}
	void FindNormal(Vector3 &norm) const
	{
		Vector3 e0 = pos[1] - pos[0];
		Vector3 e1 = pos[2] - pos[1];
		norm = -e0.cross(e1);
		norm.normalize();
	}
	void ConvertToEdgeForm()
	{
		Tri t = *this;
		// b - a
		pos[0] = t.pos[1] - t.pos[0];
		// c - a
		pos[1] = t.pos[2] - t.pos[0];
		// a
		pos[2] = t.pos[0];
	}
};


class Ray
{
public:
	Vector3 o; // origin
	Vector3 d; // direction

	bool TestIntersect(const Tri &tri, float &t) const;
	bool TestIntersect_EdgeForm(const Tri &tri_edge, float &t) const;

	void FindIntersect(const Tri &tri, float t, float &u, float &v, float &w) const;
	bool IntersectAAPlane(int axis, Vector3 &pt, float epsilon = 0.0001f) const;

};

// a ray in flight .. in the ray bank to be processed
struct FHit
{
	void SetNoHit() {tx = UINT16_MAX;}
	bool IsNoHit() const {return tx == UINT16_MAX;}
	// texel hit points
	uint16_t tx, ty;
};

struct FColor
{
	float r, g, b;
	void Set(float v) {r = v; g = v; b = v;}
	void Set(float rr, float gg, float bb) {r = rr; g = gg; b = bb;}
	void Set(const Color &col) {r = col.r; g = col.g; b = col.b;}
	float Max() const {return MAX(r, MAX(g, b));}
	FColor operator*(float v) const {FColor s; s.r = r * v; s.g = g * v; s.b = b * v; return s;}
	FColor operator/(float v) const {FColor s; s.r = r / v; s.g = g / v; s.b = b / v; return s;}
	FColor &operator+=(const FColor &v) {r += v.r; g += v.g; b += v.b; return *this;}
	FColor operator*(const FColor &o) const {FColor s; s.r = r * o.r; s.g = g * o.g; s.b = b * o.b; return s;}
};

struct FRay
{
	Ray ray;
	int num_rays_left;

	// the color of the ray at the moment
	FColor color;

	// the color of the ray that is reflected after hitting a surface
	FColor bounce_color;

	// hit texel
	FHit hit;
};


class Line2
{
public:
	Vector2 a;
	Vector2 b;
	float Distance(const Vector2 &pt) const
	{
		Vector2 dir = b - a;
		Vector2 norm = Vector2(dir.y, -dir.x);
		norm.normalize();

		Vector2 test_vec = pt - a;
		return test_vec.dot(norm);
	}
};

class UVTri
{
public:
	Vector2 uv[3];

	void FlipWinding() {Vector2 temp = uv[0]; uv[0] = uv[2]; uv[2] = temp;}
	void FindUVBarycentric(Vector2 &res, const Vector3 &bary) const {FindUVBarycentric(res, bary.x, bary.y, bary.z);}
	void FindUVBarycentric(Vector2 &res, float u, float v, float w) const;
	bool ContainsPoint(const Vector2 &pt, float epsilon = 0.0f) const;
	bool ContainsTexel(int tx, int ty, int width, int height) const;

	void FindBarycentricCoords(const Vector2 &pt, float &u, float &v, float &w) const;
	void FindBarycentricCoords(const Vector2 &pt, Vector3 &bary) const {FindBarycentricCoords(pt, bary.x, bary.y, bary.z);}
	bool IsWindingCW() const {return CalculateTwiceArea() < 0.0f;}
	float CalculateTwiceArea() const
	{
		const Vector2 &a = uv[0]; const Vector2  &b = uv[1]; const Vector2 &c = uv[2];
		return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
	}
};



bool BarycentricInside(float u, float v, float w)
{
	if ((u < 0.0f) || (u > 1.0f) ||
	(v < 0.0f) || (v > 1.0f) ||
	(w < 0.0f) || (w > 1.0f))
		return false;

	return true;
}

bool BarycentricInside(const Vector3 &bary)
{
	return BarycentricInside(bary.x, bary.y, bary.z);
}

bool IsMeshInstanceSuitable(const MeshInstance &mi)
{
	// must be set to bake in lightmap
	if (!mi.get_flag(GeometryInstance::FLAG_USE_BAKED_LIGHT))
		return false;

	Ref<Mesh> rmesh = mi.get_mesh();
	Array arrays = rmesh->surface_get_arrays(0);
	if (!arrays.size())
		return false;

	PoolVector<Vector3> verts = arrays[VS::ARRAY_VERTEX];
	if (!verts.size())
		return false;

	PoolVector<Vector3> norms = arrays[VS::ARRAY_NORMAL];
	if (!norms.size())
		return false;

//	PoolVector<int> indices = arrays[VS::ARRAY_INDEX];
//	if (!indices.size())
//		return false;

	PoolVector<Vector2> uv1s = arrays[VS::ARRAY_TEX_UV];
	PoolVector<Vector2> uv2s = arrays[VS::ARRAY_TEX_UV2];

	bool uv1 = uv1s.size() != 0;
	bool uv2 = uv2s.size() != 0;

	if (!uv1 && !uv2)
		return false;

	return true;
}

float Triangle_CalculateTwiceAreaSquared(const Vector3 &a, const Vector3 &b, const Vector3 &c)
{
	// compute the area squared. Which is the result of the cross product of two edges. If it's near zero, bingo
	Vector3 edge1 = b-a;
	Vector3 edge2 = c - a;

	Vector3 vec = edge1.cross(edge2);
	return vec.length_squared();
}


void CheckForInvalidIndices(LVector<int> &r_invalid_tris, const PoolVector<int> &inds, const PoolVector<Vector3> &verts)
{
	// check for invalid tris
//	bool res = true;

	int nTris = inds.size() / 3;
	int indCount = 0;
	for (int t=0; t<nTris; t++)
	{
		int i0 = inds[indCount++];
		int i1 = inds[indCount++];
		int i2 = inds[indCount++];

		bool ok = true;
		if (i0 == i1) ok = false;
		if (i1 == i2) ok = false;
		if (i0 == i2) ok = false;

		// check positions
		if (ok)
		{
			// vertex positions
			const Vector3 &p0 = verts[i0];
			const Vector3 &p1 = verts[i1];
			const Vector3 &p2 = verts[i2];

			float area = Triangle_CalculateTwiceAreaSquared(p0, p1, p2);
			if (area < 0.00001f)
			{
				print_line("\t\tdetected zero area triangle, ignoring");
				ok = false;
			}
		}

		if (ok)
		{
			//copy.push_back(i0); copy.push_back(i1); copy.push_back(i2);
		}
		else
		{
			r_invalid_tris.push_back(t);
			//res = false;
		}
	}

	// if any were invalid, copy the new list
//	if (!res)
//	{
//		indices.resize(copy.size());
//		for (int n=0; n<copy.size(); n++)
//		{
//			indices.set(n, copy[n]);
//		}
//	}

//	return res;
}

bool EnsureIndicesValid(PoolVector<int> &indices, const PoolVector<Vector3> &verts)
{
	// no indices? create some
	if (!indices.size())
	{
		// indices are blank!! let's create some, assuming the mesh is using triangles
		indices.resize(verts.size());
		PoolVector<int>::Write write = indices.write();
		int * pi = write.ptr();

		for (int n=0; n<verts.size(); n++)
		{
			*pi = n;
			pi++;
		}
	}

	LVector<int> invalid_tris;
	CheckForInvalidIndices(invalid_tris, indices, verts);

	if (!invalid_tris.size())
		return true;

	// we have found invalid tris.

	// check for duplicated inds in a triangle.
	LVector<int> copy;
	copy.reserve(indices.size());

	int nTris = indices.size() / 3;
	int indCount = 0;

	int invalid_count = 0;
	int next_invalid_id = invalid_tris[0];

	for (int t=0; t<nTris; t++)
	{
		// is this tri invalid?
		if (t == next_invalid_id)
		{
			invalid_count++;
			if (invalid_count < invalid_tris.size())
			{
				next_invalid_id = invalid_tris[invalid_count];
			}
			continue;
		}


		int i0 = indices[indCount++];
		int i1 = indices[indCount++];
		int i2 = indices[indCount++];
		copy.push_back(i0); copy.push_back(i1); copy.push_back(i2);
	}

	// copy the new list
	indices.resize(copy.size());
	{
		PoolVector<int>::Write write = indices.write();
		int * pi = write.ptr();

		for (int n=0; n<copy.size(); n++)
		{
			pi[n] = copy[n];
		}
	}

	return false;
}


// somewhat more complicated test to try and get all the tris that hit any part of the texel.
// really this should be a box / tri test.
inline bool UVTri::ContainsTexel(int tx, int ty, int width, int height) const
{
	float unit_x = 1.0f / width;
	float unit_y = 1.0f / height;
	float fx = tx * unit_x;
	float fy = ty * unit_y;
	float half_x = unit_x * 0.5f;
	float half_y = unit_y * 0.5f;

	// centre first as most likely
	if (ContainsPoint(Vector2(fx + half_x, fy + half_y)))
		return true;

	// 4 corners
	if (ContainsPoint(Vector2(fx, fy)))
		return true;
	if (ContainsPoint(Vector2(fx + unit_x, fy + unit_y)))
		return true;
	if (ContainsPoint(Vector2(fx + unit_x, fy)))
		return true;
	if (ContainsPoint(Vector2(fx, fy + unit_y)))
		return true;

	return false;
}

inline bool UVTri::ContainsPoint(const Vector2 &pt, float epsilon) const
{
	const Vector2 &a = uv[0];
	const Vector2 &b = uv[1];
	const Vector2 &c = uv[2];

	// if the point is in front of any of the edges, it is outside the tri,
	// if not, it must be inside

	// N.B. the triangle must be ordered clockwise for this to work
	Line2 edge;
	edge.a = b;
	edge.b = a;

	if (edge.Distance(pt) < epsilon)
		return false;

	edge.a = c;
	edge.b = b;

	if (edge.Distance(pt) < epsilon)
		return false;

	edge.a = a;
	edge.b = c;

	if (edge.Distance(pt) < epsilon)
		return false;

	return true;
}

inline void UVTri::FindBarycentricCoords(const Vector2 &pt, float &u, float &v, float &w) const
{
	const Vector2 &a = uv[0];
	const Vector2 &b = uv[1];
	const Vector2 &c = uv[2];
	Vector2 v0 = b - a, v1 = c - a, v2 = pt - a;
	float d00 = v0.dot(v0);
	float d01 = v0.dot(v1);
	float d11 = v1.dot(v1);
	float d20 = v2.dot(v0);
	float d21 = v2.dot(v1);
	float invDenom = 1.0f / (d00 * d11 - d01 * d01);
	v = (d11 * d20 - d01 * d21) * invDenom;
	w = (d00 * d21 - d01 * d20) * invDenom;
	u = 1.0f - v - w;
}


inline void UVTri::FindUVBarycentric(Vector2 &res, float u, float v, float w) const
{
	res = (uv[0] * u) + (uv[1] * v) + (uv[2] * w);
}
//////////////////////////////////////////////////////////

inline void Tri::FindBarycentric(const Vector3 &pt, float &u, float &v, float &w) const
{
	const Vector3 &a = pos[0];
	const Vector3 &b = pos[1];
	const Vector3 &c = pos[2];

	// we can alternatively precache
	Vector3 v0 = b - a, v1 = c - a, v2 = pt - a;
	float d00 = v0.dot(v0);
	float d01 = v0.dot(v1);
	float d11 = v1.dot(v1);
	float d20 = v2.dot(v0);
	float d21 = v2.dot(v1);
	float invDenom = 1.0f / (d00 * d11 - d01 * d01);
	v = (d11 * d20 - d01 * d21) * invDenom;
	w = (d00 * d21 - d01 * d20) * invDenom;
	u = 1.0f - v - w;
}

inline void Tri::InterpolateBarycentric(Vector3 &pt, float u, float v, float w) const
{
	const Vector3 &a = pos[0];
	const Vector3 &b = pos[1];
	const Vector3 &c = pos[2];
	pt.x = (a.x * u) + (b.x * v) + (c.x * w);
	pt.y = (a.y * u) + (b.y * v) + (c.y * w);
	pt.z = (a.z * u) + (b.z * v) + (c.z * w);
}


// returns barycentric, assumes we have already tested and got a collision
inline void Ray::FindIntersect(const Tri &tri, float t, float &u, float &v, float &w) const
{
	Vector3 pt = o + (d * t);
	tri.FindBarycentric(pt, u, v, w);
}

// returns false if doesn't cut x.
// the relevant axis should be filled on entry (0 is x, 1 is y, 2 is z) with the axis aligned plane constant.
bool Ray::IntersectAAPlane(int axis, Vector3 &pt, float epsilon) const
{
	// x(t) = o.x + (d.x * t)
	// y(t) = o.y + (d.y * t)
	// z(t) = o.z + (d.z * t)

	// solve for x
	// find t
	// d.x * t = x(t) - o.x
	// t = (x(t) - o.x) / d.x

	// just aliases
	//	const Vector3 &o = ray.o;
	//	const Vector3 &d = ray.d;

	switch (axis)
	{
	case 0:
		{
			if (fabsf(d.x) < epsilon)
				return false;

			float t = pt.x - o.x;
			t /= d.x;

			pt.y = o.y + (d.y * t);
			pt.z = o.z + (d.z * t);
		}
		break;
	case 1:
		{
			if (fabsf(d.y) < epsilon)
				return false;

			float t = pt.y - o.y;
			t /= d.y;

			pt.x = o.x + (d.x * t);
			pt.z = o.z + (d.z * t);
		}
		break;
	case 2:
		{
			if (fabsf(d.z) < epsilon)
				return false;

			float t = pt.z - o.z;
			t /= d.z;

			// now we have t, we can calculate y and z
			pt.x = o.x + (d.x * t);
			pt.y = o.y + (d.y * t);
		}
		break;
	default:
		break;
	}

	return true;
}


// same as below but with edges pre calced
inline bool Ray::TestIntersect_EdgeForm(const Tri &tri_edge, float &t) const
{
	//Edge1, Edge2
	const Vector3 &e1 = tri_edge.pos[0];
	const Vector3 &e2 = tri_edge.pos[1];
	const Vector3 &a = tri_edge.pos[2];

	const float cfEpsilon = 0.000001f;
	Vector3 P, Q, T;
	float det, inv_det, u, v;

	//Begin calculating determinant - also used to calculate u parameter
	P = d.cross(e2);

	//if determinant is near zero, ray lies in plane of triangle
	det = e1.dot(P);

	//NOT CULLING
	if(det > -cfEpsilon && det < cfEpsilon) return false;//IR_NONE;
	inv_det = 1.f / det;

	//calculate distance from V1 to ray origin
	T = o - a;

	//Calculate u parameter and test bound
	u = T.dot(P) * inv_det;

	//The intersection lies outside of the triangle
	if(u < 0.f || u > 1.f) return false; // IR_NONE;

	//Prepare to test v parameter
	Q = T.cross(e1);

	//Calculate V parameter and test bound
	v = d.dot(Q) * inv_det;

	//The intersection lies outside of the triangle
	if(v < 0.f || u + v  > 1.f) return false; //IR_NONE;

	t = e2.dot(Q) * inv_det;

	if(t > cfEpsilon)
	{ //ray intersection
		//	*out = t;
		return true; //IR_HIT;
	}

	// No hit, no win
	return false; //IR_BEHIND;
}

// moller and trumbore test (actually calculated uv and t but we are not returning uv)
inline bool Ray::TestIntersect(const Tri &tri, float &t) const
{
	const Vector3 &a = tri.pos[0];
	const Vector3 &b = tri.pos[1];
	const Vector3 &c = tri.pos[2];

	const float cfEpsilon = 0.000001f;
	//#define EPSILON 0.000001
	Vector3 e1, e2;  //Edge1, Edge2
	Vector3 P, Q, T;
	float det, inv_det, u, v;
	//	float t;

	//Find vectors for two edges sharing V1
	e1 = b - a;
	e2 = c - a;

	//Begin calculating determinant - also used to calculate u parameter
	P = d.cross(e2);

	//if determinant is near zero, ray lies in plane of triangle
	det = e1.dot(P);

	//NOT CULLING
	if(det > -cfEpsilon && det < cfEpsilon) return false;//IR_NONE;
	inv_det = 1.f / det;

	//calculate distance from V1 to ray origin
	T = o - a;

	//Calculate u parameter and test bound
	u = T.dot(P) * inv_det;

	//The intersection lies outside of the triangle
	if(u < 0.f || u > 1.f) return false; // IR_NONE;

	//Prepare to test v parameter
	Q = T.cross(e1);

	//Calculate V parameter and test bound
	v = d.dot(Q) * inv_det;

	//The intersection lies outside of the triangle
	if(v < 0.f || u + v  > 1.f) return false; //IR_NONE;

	t = e2.dot(Q) * inv_det;

	if(t > cfEpsilon)
	{ //ray intersection
		//	*out = t;
		return true; //IR_HIT;
	}

	// No hit, no win
	return false; //IR_BEHIND;
}


}
