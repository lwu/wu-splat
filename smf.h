#ifndef SMFDEMO_INCLUDED // -*- C++ -*-
#define SMFDEMO_INCLUDED

/************************************************************************

  Example code for reading SMF files.

  $Id: smf.h,v 1.1 2001/02/26 03:36:55 garland Exp $

 ************************************************************************/

#include <vector>
#include <gfx/vec4.h>
#include "interleaved.h"
#include "abnormal.h"

typedef TVec3<int> Face;

struct Vert
{
	Vert() {}
	Vert(const Vec3f& v) : pos(v), radius(0.0f), cone_dot(1.0f) {
	}

	Vec3f	pos;
	float	radius;
	Vec3f	normal;
#if defined(__ABNORMAL)
	Norm	abnormal;
#endif // __ABNORMAL
	float	cone_dot;
};

typedef InterleavedArray<Vec3f, Vert, offsetof(Vert,pos)>    VertexArray;
typedef InterleavedArray<float, Vert, offsetof(Vert,radius)> RadiusArray;
typedef InterleavedArray<Vec3f, Vert, offsetof(Vert,normal)> NormalArray;

typedef std::vector<Vert>::iterator VertIt;

template <typename T>
T get_interleaved(VertIt begin, VertIt end)
{
	return T(begin, end);
}

class TriMesh
{
public:
	VertexArray get_vertex_array() {
		return get_interleaved<VertexArray>(vertex.begin(), vertex.end());
	}

	RadiusArray get_radius_array() {
		return get_interleaved<RadiusArray>(vertex.begin(), vertex.end());
	}

	NormalArray get_normal_array() {
		return get_interleaved<NormalArray>(vertex.begin(), vertex.end());
	}

	std::vector<Vert> vertex;
    std::vector<Face> face;	

	void clear() {
		vertex.clear();
		face.clear();		
	}
};

extern void read_smf_object(TriMesh& tm, const char *filename);
extern void read_smf_object_fast(TriMesh& tm, const char *filename);

// SMFDEMO_INCLUDED
#endif
