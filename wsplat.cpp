#include <gfx/gl.h>
#include <GL/glut.h>
#include <iostream>
#include <float.h>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include "jacobi.h"
#include "wsplat.h"
#include "smf.h"
#include "abnormal.h"
#include "util.h"

using namespace std;

WSplat::WSplat() 
{	
	// GUI options
	draw_splats = true;
	draw_points = false;
	draw_spheres = false;
	draw_triangles = false;
	draw_normals = false;
	round_points = true;

	backface_culling = true;
	always_subdivide = false;
	change_colors = false;
	psychedelic = false;
	abnormal = false;
	project_size = true;
	allow_quads = true;
	quads_only = false;

	draw_edges_only = false;
	silhouette_edge_thresh = 0.3f;	
	max_point_size = 63;

	// tree-generation options
	oriented_partition = false;
	subdivide_large_tris = true;
	subdivide_a_lot = true;

	splat_tree = 0;
	num_nodes = 0;
	num_ignored = 0;
	num_zero_norms = 0;

	max_pixel_error = 4.0f;
	num_points_drawn = 0;	
}

void WSplat::destroy() 
{
	delete splat_tree;
	splat_tree = 0;

	mesh.clear();

	num_nodes = 0;
	num_ignored = 0;
	num_zero_norms = 0;
}

void WSplat::get_bounding_sphere(Vec3f& centre, float& radius)
{
	centre = centroid;
	radius = splat_tree ? splat_tree->node.radius : 0.0f;
}

ColorChanger::ColorChanger()
{
	current_index = 0;

	const int N = 3;

	for (float r=0.0f; r <= 1.0f; r += 1.0f/(float)N) {
		for (float g=0.0f; g <= 1.0f; g += 1.0f/(float)N) {
			for (float b=0.0f; b <= 1.0f; b += 1.0f/(float)N) {
				colors.push_back( Vec3f(r,g,b) );
			}
		}
	}

	colors[0] = Vec3f(1.0f, 1.0f, 1.0f);
}

void ColorChanger::change()
{
	current_index++;
	if (current_index > colors.size()) {
		current_index = 0;
	}

	set_color();
}

inline void ColorChanger::set_color() const 
{
	glColor3fv(colors[current_index]);
}

void ColorChanger::reset_index() 
{
	current_index = 0;
}

float dist_sqr(const Vec3f& v0, const Vec3f& v1)
{	
	return norm2(v1 - v0);
}



int depth = 0;
int max_depth = 0;

// Load a .SMF file, building up the WSplat tree
void WSplat::load(const char* fn)
{
	Norm::generate_table();

	if (!fn) {
		cout << "Invalid filename!" << endl;
		return;
	}

	// cleanup old splat data
	static bool cleanup = false;
	if (cleanup) {
		cout << "\nDeracinating splat tree... ";
	}
	cleanup = true;
	destroy();

	// disk I/O
	Timer t;
	cout << "\nReading '" << fn << "' from disk... ";	
	read_smf_object_fast(mesh, fn);
	double dt = t.get_elapsed();
	cout.precision(2);
	cout << "(" << dt << "s)";
	
	cout << "\nRead " << mesh.vertex.size() << " " << verts(mesh.vertex.size()) << ", "
		<< mesh.face.size() << " faces." << endl;

	VertexArray va = mesh.get_vertex_array();	
	centroid = accumulate( va.begin(), va.end(), Vec3f() ) / va.size();	

	// preprocessing
	Timer t_process;
	if (subdivide_large_tris) {
		tessellate_faces();
	}
	compute_bounding_spheres();
	store_normals();
	grow_tree();

	int bytes = num_nodes * sizeof(Vert);
	int num_kb = bytes >> 10;

	cout << "\n\tCreated " << num_nodes << " nodes (" << num_kb << "KB total)" << endl;
	if (num_ignored > 0) {
		cout << "\tIgnored " << num_ignored << " 0-norm " << verts(num_ignored) << endl;
	}
	if (num_zero_norms > 0) {
		cout << "\tFound " << num_zero_norms << " 0-norm " << verts(num_zero_norms)
			<< " in the final tree." << endl;
	}
	cout << "\tTree depth = " << max_depth << "." << endl;
	cout << "\tPre-processing time " << t_process.get_elapsed() << "s" << endl;

	cout.precision(5);
}


// Functors for computing bounding box info

template <typename T, int DIM>
struct MinMax
{
	MinMax() {
		fill( &min[0], &min[DIM], 1e8 );
		fill( &max[0], &max[DIM], -1e8 );
	}

	T min[DIM], max[DIM];
};

template <typename Container, typename T, int DIM>
struct ComputeMinMax
{
	void operator()(const Container& c) {
		for (int i=0; i < DIM; i++) {
			T data = c[i];
			mm.min[i] = min(mm.min[i], data);
			mm.max[i] = max(mm.max[i], data);
		}
	}

	operator MinMax<T,DIM>() {
		return mm;
	}
private:
	MinMax<T,DIM> mm;
};

template <typename T, int DIM>
int LongestAxis(const MinMax<T,DIM>& mm)
{
	T maxi(0);
	int max_axis = 0;

	for (int i=0; i < DIM; i++) {
		T delta = mm.max[i] - mm.min[i];
		if (delta > maxi) {
			maxi = delta;
			max_axis = i;
		}
	}

	return max_axis;
}

struct ShorterThan {
	ShorterThan(float t, VertexArray v) : verts(v) {
		thresh = t;		
	}

	bool operator()(const Face& f) {
		return (max_len(f) <= thresh);
	}

	float max_len(const Face& f) {
		float d01 = dist_sqr( verts[f[0]], verts[f[1]] );
		float d12 = dist_sqr( verts[f[1]], verts[f[2]] );
		float d20 = dist_sqr( verts[f[2]], verts[f[0]] );

		return max( d01, max(d12, d20) );
	}
private:
	float thresh;
	VertexArray verts;
};

// Recursively subdivides longest edge of triangle until threshold met
void tessellate_triangle(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2,
						vector<Vec3f>& verts, vector<Face>& faces, float thresh)
{
	float d01 = dist_sqr( v0, v1 ); // Yes, this could be much faster..
	float d12 = dist_sqr( v1, v2 );
	float d20 = dist_sqr( v2, v0 );
	
	float max_len = max(d01, max(d12, d20));

	if (max_len <= thresh) {
		// add new triangle
		int new_index = verts.size();
		verts.push_back(v0);
		verts.push_back(v1);
		verts.push_back(v2);
		faces.push_back( Face(new_index, new_index+1, new_index+2) );

		//cout << "New tri: " << v0 << "," << v1 << "," << v2 << endl;
	} else {
		int longest[2];

		if (d01 > d12) {
			if (d01 > d20) {
				// 01 longest
				longest[0] = 0; longest[1] = 1;
			} else {
				// 20 longest
				longest[0] = 2; longest[1] = 0;
			}
		} else {
			if (d12 > d20) {
				// 12 longest
				longest[0] = 1; longest[1] = 2;
			} else {
				// 20 longest
				longest[0] = 2; longest[1] = 0;
			}
		}

		int shortest = 3 - (longest[0] + longest[1]);
		const Vec3f* vecs[3] = {&v0, &v1, &v2};

		// subdivide							
		Vec3f midpt = 0.5f * (*vecs[longest[0]] + *vecs[longest[1]]);
		tessellate_triangle(*vecs[longest[0]], midpt, *vecs[shortest], verts, faces, thresh);
		tessellate_triangle(midpt, *vecs[longest[1]], *vecs[shortest], verts, faces, thresh);
	}
}

// Examines model and tessellates triangles that are too large
//  in comparison to longest AABB axis
void WSplat::tessellate_faces()
{	
	cout << "Examining face tessellation statistics... ";

	float max_edge_len = 0.0f;
	float min_edge_len = +1e6;
	double avg_edge_len = 0.0;
	VertexArray verts = mesh.get_vertex_array();

	MinMax<float,3> mm = for_each( verts.begin(), verts.end(),
		ComputeMinMax<Vec3f,float,3>() );
	
	int max_index = LongestAxis(mm);

	float max_axis_dim = mm.max[max_index] - mm.min[max_index];

	cout << "Length of max. bbox axis = " << max_axis_dim << endl;

	int num_faces = mesh.face.size();	

	int cnt=0;

	float MAGIC_SCALAR = subdivide_a_lot ? 
		(1.0f / 10000.0f) : (1.0f / 1000.0f);
	const float magic_thresh = max_axis_dim*MAGIC_SCALAR;
	
	ShorterThan lengthp(magic_thresh, verts);

	vector<Face>::iterator it = 
		partition( mesh.face.begin(), mesh.face.end(), lengthp);

	cnt = distance(it, mesh.face.end() );
	cout << (float)cnt/(float)num_faces*100.0f << "% of faces are larger than "
		<< 100.0f*MAGIC_SCALAR << "% size of longest AABB axis" << endl;

	vector<Face> faces(it, mesh.face.end());

	mesh.face.erase(it, mesh.face.end());
	
	int tri_cnt = 0, face_cnt = 0;

	for (int i=0; i < faces.size(); i++) {
		VertexArray vs = mesh.get_vertex_array();
		const Face& f = faces[i];
		vector<Vec3f> verts;
		vector<Face>  faces;
		tessellate_triangle(vs[f[0]], vs[f[1]], vs[f[2]],
			verts, faces, magic_thresh);

		tri_cnt += verts.size();
		face_cnt += faces.size();

		int new_index0 = mesh.vertex.size();
		copy(verts.begin(), verts.end(),
			back_inserter(mesh.vertex));

		for (int i=0; i < faces.size(); i++) {
			faces[i][0] += new_index0;
			faces[i][1] += new_index0;
			faces[i][2] += new_index0;
			mesh.face.push_back(faces[i]);
		}
	}

	cout << "Added " << tri_cnt << " new vertices, " << face_cnt << " new faces." << endl;
	//cout << mesh.face.size() << " faces, " << mesh.vertex.size() << " vertices." << endl;
}


struct my_sqrt : unary_function<double,double> {
	double operator()(double x) { return sqrt(x); }
};

// Compute minimal per-vertex bounding spheres
void WSplat::compute_bounding_spheres() 
{
	cout << "Computing per-vertex bounding spheres... ";

	int num_vertices = mesh.vertex.size();
	int num_faces = mesh.face.size();

	VertexArray verts = mesh.get_vertex_array();
	RadiusArray rads  = mesh.get_radius_array();

	// Set minimum bounding sphere radius of a vertex
	//  as the maximum length of the edge leaving that vertex
	for (int i=0; i < num_faces; i++) {
		const Face& f = mesh.face[i];

		float d01 = dist_sqr( verts[f[0]], verts[f[1]] );
		float d12 = dist_sqr( verts[f[1]], verts[f[2]] );
		float d20 = dist_sqr( verts[f[2]], verts[f[0]] );

		rads[f[0]] = max(rads[f[0]], 
			max(d01, d20));
		rads[f[1]] = max(rads[f[1]], 
			max(d01, d12));
		rads[f[2]] = max(rads[f[2]], 
			max(d12, d20));
	}	

	// sqrt() since we previously only computed distance^2
	transform(rads.begin(), rads.end(),
		rads.begin(), my_sqrt() );	
}

bool is_zero(const Vec3f& v)
{
	return (v[0] == 0.0f && v[1] == 0.0f && v[2] == 0.0f);
}

// Store normals for each vertex using face data
void WSplat::store_normals()
{
	vector<int> normal_count( mesh.vertex.size() );

	NormalArray na = mesh.get_normal_array();
	VertexArray va = mesh.get_vertex_array();

	int face_cnt = 0;

	// For every face, add normal to each vertex
	vector<Face>::iterator it;
	for ( it = mesh.face.begin(); it != mesh.face.end(); ++it) {
		const Face& f = *it;

		Vec3f n = triangle_normal( va[f[0]], va[f[1]], va[f[2]] );
		na[f[0]] += n;
		na[f[1]] += n;
		na[f[2]] += n;
		
		++normal_count[f[0]];
		++normal_count[f[1]];
		++normal_count[f[2]];
	}

	// Find average normal by dividing by # of normals added
	NormalArray::iterator ni = na.begin();
	int i=0;
	
	int num_degenerates = 0;
	int num_zero_norms = 0;

	for ( ; ni != na.end(); ++ni) {
		if (normal_count[i] == 0) {
			++num_degenerates;
		} else {
			*ni /= normal_count[i];		
			*ni /= norm(*ni);
#if defined(__ABNORMAL)
			mesh.vertex[i].abnormal = *ni;
#endif
		}
		
		if (is_zero(*ni)) {
			num_zero_norms++;
		}

		++i;
	}

	if (num_degenerates > 0 || num_zero_norms > 0) {
		cout << "\n";
	}

	if (num_degenerates > 0) {
		cout << "\tWarning: found " << num_degenerates 
			<< " faceless " << verts(num_degenerates) << endl;
	}

	if (num_zero_norms > 0) {
		cout << "\tWarning: found " << num_zero_norms << " "
			<< verts(num_zero_norms) << " with 0 norm." << endl;
	}
}


struct PlaneContainsPoint {
	PlaneContainsPoint(const Vec3f& origin, const Vec3f& dir) {
		plane_origin = origin;
		plane_dir = dir;
	}

	bool operator()(const Vert& v) {
		Vec3f dv = v.pos - plane_origin;
		return (dv*plane_dir >= 0.0f);
	}
private:
	Vec3f plane_origin;
	Vec3f plane_dir;
};

struct AboveCoord {
	AboveCoord(int coord, float value) {
		coord_index = coord;
		thresh_value = value;
	};

	bool operator()(const Vert& v) {
		return (v.pos[coord_index] < thresh_value);
	}
private:
	int   coord_index;
	float thresh_value;
};

void WSplat::grow_tree()
{
	cout << "Growing tree...";
	splat_tree = build_tree( mesh.vertex.begin(), mesh.vertex.end() );

	// draw normals as a fraction of largest bounding sphere radius
	normal_scale = splat_tree->node.radius / 10.0f;
}


struct Counter {
	Counter() {
		depth++;
		max_depth= max(max_depth, depth);
	}

	~Counter() {
		depth--;
	}
};


float Camera::calc_projected_dist(const Vec3f& sphere_pos)
{
	Vec3f dp = sphere_pos - this->pos;
	float dist = dir * dp;

	return dist;
}

void WSplat::set_camera(const Vec3f& pos, const Vec3f& dir, float screen_width)
{
	cam.pos = pos;
	cam.dir = dir / norm(dir);
	cam.screen_width = screen_width;

	// unit vectors
	Vec3f up(0.0f, 1.0f, 0.0f);
	cam.left = cross(up, cam.dir);
	cam.up = cross(cam.dir, cam.left);

	cam.left /= norm(cam.left);
	cam.up   /= norm(cam.up);
}

void WSplat::set_error(float tau) 
{
	max_pixel_error = tau;
}



struct sum_covariance {
	sum_covariance(const Vec3f& m) {
		mean = m;
	}
	void operator()(const Vec3f& v) {
		Vec3f df = v - mean;
		Vec3 d(df);
		CV += outer_product(d,d);
	}
	
	operator Mat3() {
		return CV;
	}
private:
	Mat3 CV;
	Vec3f mean;
};


VertTree* WSplat::build_tree(VertIt begin, VertIt end)
{
	Counter c;
	num_nodes++;

	size_t dist = distance(begin, end);

	if (dist == 0) {
		return 0;
	} else if (dist == 1) {
		if (is_zero(begin->normal)) {
			num_ignored++;
			return 0;
		}

		VertTree* tree = new VertTree( *begin );
#if defined(__ABNORMAL)
		tree->node.abnormal = tree->node.normal;
#endif
		return tree;
	}

	VertexArray verts(begin, end);

	Vec3f avg;
	vector<Vert>::iterator pos;

	if (oriented_partition) {
		avg = accumulate( verts.begin(), verts.end(), Vec3f() ) / verts.size();

		Mat3 CV =  // Compute the covariance matrix
			for_each( verts.begin(), verts.end(), sum_covariance(avg) );
		
		// And extract the eigenvalues/eigenvectors
		Vec3 evals, evecs[3];
		if( !eigen(CV, evals, evecs) ) {
			// couldn't extract eigenvectors
			oriented_partition = false;
		} else {
			// evecs[0] is the direction to split in			
			Vec3f dir(evecs[0]);

			// Partition
			pos = partition( begin, end, 
				PlaneContainsPoint(avg, dir));
		}
	} 
	
	// Find longest axis
	if (!oriented_partition) {	
		MinMax<float,3> mm = for_each( verts.begin(), verts.end(),
			ComputeMinMax<Vec3f,float,3>() );
		
		int max_index = LongestAxis(mm);
		
		for (int i=0; i < 3; i++) {		
			avg[i] = (mm.max[i] + mm.min[i]) * 0.5f;	
		}	

		// Partition
		pos = partition( begin, end, 
			AboveCoord(max_index, avg[max_index]) );
	}
	

	if (pos == begin || pos == end) {
		// Pleonastic vertices. Just return one vertex.
		// (In the case of retessellation, we generate
		//  duplicate vertices. We want to use those
		//  instead of faceless vertices we left behind..)
		VertIt it = begin;

		for ( ; it != end; ++it) {
			if (!is_zero(pos->normal)) {
				return new VertTree(*it);
			}
		}
		
		// Better return nothing than a ab-normal vertex
		return 0;
	}


	VertTree* left  = build_tree(begin, pos);
	VertTree* right = build_tree(pos, end);
	
	VertTree* tree = 0;
	if ( !(left || right) ) {
		return tree;
	} else {
		tree = new VertTree( avg );
		tree->left_subtree = left;
		tree->right_subtree = right;
	}	

	// compute new normal
	Vec3f normal;
	int num_normal = 0;

	float dist_left;
	if (left) {
		dist_left = dist_sqr( left->node.pos, avg );
		normal += left->node.normal;
		num_normal++;		
	} else {
		dist_left = 0.0f;		
	}

	float dist_right;
	if (right) {
		dist_right = dist_sqr( right->node.pos, avg );
		normal += right->node.normal;
		num_normal++;
	} else {
		dist_right = 0.0f;	
	}	

	assert(left || right);
	assert(num_normal > 0);

	if (num_normal == 2) {
		if (!is_zero(normal)) {
			normal *= 0.5f;
			normal /= norm(normal);
		} else {
			// this only happens 1/1e6 of the time
			normal = left->node.normal;
		}
	}

	if (is_zero(normal)) {
		num_zero_norms++;
	}

	tree->node.normal = normal;
#if defined(__ABNORMAL)
	tree->node.abnormal = normal;
#endif

	assert (!_isnan(normal[0]));

	// store normal cone dot product 
	float dot_left;
	if (left) {
		// avoid NaN for acos(1.0+epsilon)
		float dot = normal*left->node.normal * 0.999f;
		dot_left = cos (acos(dot) + acos(left->node.cone_dot));
	} else {
		dot_left = 0.0f;
	}

	float dot_right;
	if (right) {
		float dot = normal*right->node.normal * 0.999f;
		dot_right = cos (acos(dot) + acos(right->node.cone_dot));
	} else {
		dot_right = 0.0f;
	}
	
	tree->node.cone_dot = min(dot_left, dot_right);

	assert (!_isnan(tree->node.cone_dot));

	VertTree* larger = dist_right > dist_left ? right : left;
	float max_dist = sqrt( max(dist_right, dist_left) );
	
	max_dist += larger->node.radius;

	tree->node.radius = max_dist;

	return tree;
}



void WSplat::draw()
{
	int num_faces = mesh.face.size();
	int num_elements = num_faces*3;

	glPointSize(6.0f);
	if (round_points) {
		glEnable(GL_POINT_SMOOTH);
		glAlphaFunc(GL_GREATER, 0.5f);
		glEnable(GL_ALPHA_TEST);
	} else {
		glDisable(GL_POINT_SMOOTH);
		glDisable(GL_ALPHA_TEST);
	}

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	/*if (draw_triangles) {
		glInterleavedArrays(GL_V3F, sizeof(Vert), &mesh.vertex[0]);
		glColor3f(0.0f, 1.0f, 0.0f);
		glDrawElements(GL_TRIANGLES, num_elements, 
			GL_UNSIGNED_INT, &mesh.face[0]);		
	}*/

	if (!psychedelic) {
		glColor3f(1.0f, 1.0f, 1.0f);
		rainbow.reset_index();
	}

	if (draw_points) {		
		glInterleavedArrays(GL_V3F, sizeof(Vert), &mesh.vertex[0]);
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT, sizeof(Vert), &mesh.vertex[0].normal); 
		glDrawArrays(GL_POINTS, 0, mesh.vertex.size());
		glDisableClientState(GL_NORMAL_ARRAY);
	}

	VertexArray va = mesh.get_vertex_array();
	RadiusArray ra  = mesh.get_radius_array();

	if (draw_spheres) {
		glColor3f(1.0f, 1.0f, 0.0f);
		int num_vertices = mesh.vertex.size();
		for (int i=0; i < num_vertices; i++) {
			glPushMatrix();
			const Vec3f v( va[i] );
			glTranslatef( v[0], v[1], v[2] );
			glutSolidSphere( ra[i], 8,8 );
			glPopMatrix();
		}
	}

	num_points_drawn = 0;

	if (draw_splats) {
		if (splat_tree) {
			NormConeDir cdir = DIR_UNDECIDED;
			draw_tree(splat_tree, cdir);
		}
	}
}

void WSplat::draw_quad(const Vec3f& pos, const Vec3f& normal, float radius)
{
	Vec3f left = cam.left * radius;
	Vec3f up   = cam.up * radius;

	glBegin(GL_QUADS);
		glNormal3fv(normal);
		glVertex3fv(pos - left + up);
		glVertex3fv(pos + left + up);
		glVertex3fv(pos + left - up);
		glVertex3fv(pos - left - up);
	glEnd();
}

extern int max_user_depth;
Vec3f rotate_direction(const Vec3f& moving, float angle, const Vec3f& axis);

void WSplat::draw_tree(VertTree* tree, NormConeDir cdir) 
{
	Counter c;

	bool has_child = false;
	bool subdivide = true;

	if (backface_culling && cdir == DIR_UNDECIDED) {
		float dot = tree->node.normal * cam.dir;
		float ddot = -tree->node.cone_dot + 1.0f;
		
		if (draw_edges_only) {
			// reject normal cones that don't fall within 
			//  silhouette edge threshold, or point away from viewer
			if ((dot >= 0.0f && dot - ddot >= 0.0f/*edge_dot_thresh*/) ||
				(dot <= 0.0f && dot + ddot <= -silhouette_edge_thresh)) {
				return;
			} 
		} else {
			// normal backface culling	
			if (dot >= 0.0f && dot - ddot >= 0.0f) {		
				// all normals in cone point away; don't draw
				return;
			}
			
			if (dot <= 0.0f && dot + ddot <= 0.0f) {
				// all normals in cone are visible; don't check anymore
				cdir = DIR_TOWARD;
			}
		}
	}

	if (depth <= max_user_depth && change_colors) {
		rainbow.change();
	}

	float dist;	
	float screen_size;
	
	if (project_size) {
		dist = cam.calc_projected_dist(tree->node.pos);
		screen_size = tree->node.radius / dist * cam.screen_width * 2.0f;
	} else {
		screen_size = 4.0f;
	}

	// draw children
	if ((depth <= max_user_depth && screen_size >= max_pixel_error) || always_subdivide) {
		if (tree->left_subtree) {
			draw_tree(tree->left_subtree, cdir);
			has_child = true;
		}

		if (tree->right_subtree) {
			has_child = true;
			draw_tree(tree->right_subtree, cdir);
		}
	} 
	
	// ... or draw current node
	if (!has_child) {				
		rainbow.set_color();

		num_points_drawn++;

		if (!quads_only && (!allow_quads || screen_size < max_point_size)) {			
			// draw normal point
			glPointSize(screen_size);
			glBegin(GL_POINTS);
#if defined(__ABNORMAL)
				if (abnormal) {
					Vec3f norm = tree->node.abnormal;
					glNormal3fv(norm);
				} else 
#endif
				{
					glNormal3fv(tree->node.normal);
				}
				glVertex3fv(tree->node.pos);
			glEnd();
		} else {
			// draw quad
			draw_quad(tree->node.pos, tree->node.normal, tree->node.radius);
		}

		if (draw_normals) {	
			glBegin(GL_LINES);
				glVertex3fv(tree->node.pos);
				glVertex3fv(tree->node.pos + tree->node.normal*normal_scale);
			glEnd();				
		}
								
		if (draw_cones) {
			Vec3f norm(tree->node.normal);
			static Vec3f up(0.0f, 1.0f, 0.0f);
			Vec3f left = cross(norm, up);
			float half_theta = 0.5f * acos(tree->node.cone_dot)*180.0f/M_PI;
			Vec3f cone_dir0 = rotate_direction(norm, half_theta, left) * normal_scale;
			Vec3f cone_dir1 = rotate_direction(norm, -half_theta, left) * normal_scale;
			
			glBegin(GL_LINE_LOOP);
			glVertex3fv(tree->node.pos);
			glVertex3fv(tree->node.pos + cone_dir0);
			glVertex3fv(tree->node.pos + cone_dir1);			
			glEnd();
		}
	}
}