/************************************************************************

  Example code for reading SMF files

  $Id: smf.cxx,v 1.2 2001/02/26 03:36:55 garland Exp $

 ************************************************************************/

#pragma warning (disable: 4786) 
#pragma warning (disable: 4788) 

#include "smf.h"
#include <gfx/script.h>
#include <gfx/geom3d.h>

static TriMesh *current_mesh = NULL;

/*
static int script_v(const char *op, char *argline, Scripting *)
{
    Vec3f v;

    if(string_shift_numbers(&argline, v, 3)<3) return SCRIPT_ERR_SYNTAX;
    current_mesh->vertex.push_back(v);
    return SCRIPT_OK;
}
*/

/*
static int script_f(const char *op, char *argline, Scripting *)
{
    Face f;

    if(string_shift_numbers(&argline, f, 3)<3) return SCRIPT_ERR_SYNTAX;

    // Input indices are numbered from 1.  Subtract 1 to get 0-based
    // indices ala C++ arrays.
    f[0] -= 1; f[1] -= 1; f[2] -= 1;
    current_mesh->face.push_back(f);

	//VertexArray va = current_mesh->get_vertex_array();

    //const Vec3f& v1 = va[f[0]];
    //const Vec3f& v2 = va[f[1]];
    //const Vec3f& v3 = va[f[2]];
    //Vec3f n = triangle_normal( v1, v2, v3 );
    //current_mesh->normal.push_back(n);
    
    return SCRIPT_OK;
}
*/

 /*
static script_defn handlers[] =
{
    {"v", script_v},
    {"f", script_f},
};

REGISTER_COMMANDS(handlers);

void read_smf_object(TriMesh& tm, const char *filename)
{
	current_mesh = &tm;

    script_do_file(filename);
    
    current_mesh = NULL;    
}
 */

#include <iostream>
#include <fstream>

using namespace std;

// This is typically 3-4x faster since it only
//  reads vertices/faces and avoids script parsing overhead

void read_smf_object_fast(TriMesh& tm, const char *filename)
{
	ifstream is(filename);

	if (!is) {
		cerr << "[Couldn't open " << filename << " for reading!]" << endl;
		return;
	}

	Vec3f v;
	Face f;
	while (is) {
		char c;
		is >> c;

		switch (c) {
		case 'v':
			is >> v[0] >> v[1] >> v[2];
			tm.vertex.push_back(v);
			break;
		case 'f':
			is >> f[0] >> f[1] >> f[2];
			f[0] -= 1; f[1] -= 1; f[2] -= 1;
			tm.face.push_back(f);
			break;
		default:
			cerr << "(Parsing error?)";
			break;
		}
	}

	if (!tm.face.empty()) tm.face.pop_back();
}
