
#ifndef WSPLAT_H
#define WSPLAT_H

#include <vector>
#include <gfx/geom3d.h>
#include "smf.h"

struct VertTree {
	VertTree() {
		left_subtree = right_subtree = 0;
	}	

	VertTree(const Vert& v) : node(v) {
		left_subtree = right_subtree = 0;
	}

	~VertTree() {
		delete left_subtree;
		delete right_subtree;
	}

	Vert node;

	VertTree* left_subtree;
	VertTree* right_subtree;
};

struct Camera {
	float calc_projected_dist(const Vec3f& pos);

	Vec3f pos;
	Vec3f dir;

	Vec3f up, left;

	float screen_width;
};

class ColorChanger {
public:
	ColorChanger();

	void change();
	void reset_index();
	void set_color() const;

private:
	std::vector<Vec3f> colors;
	int current_index;
};

class WSplat {
public:
	WSplat();
	
	void load(const char* fn);
	void draw();	

	void get_bounding_sphere(Vec3f& centroid, float& radius);

	void set_camera(const Vec3f& pos, const Vec3f& dir, float screen_width);
	void set_error(float tau);

	// options exposed to GUI
	bool draw_triangles;
	bool draw_points;
	bool draw_spheres;
	bool draw_splats;
	bool draw_normals;
	bool round_points;
	bool change_colors;
	bool psychedelic;
	bool always_subdivide;
	bool abnormal;
	bool draw_cones;
	bool backface_culling;
	bool project_size;
	bool allow_quads;
	bool quads_only;
	bool draw_edges_only;
	int  max_point_size;

	float silhouette_edge_thresh;
	Camera cam;

	// tree generation options
	bool oriented_partition;
	bool subdivide_large_tris;
	bool subdivide_a_lot;

	int num_points_drawn;
private:
	TriMesh mesh;	

	// tree construction
	void tessellate_faces();
	void compute_bounding_spheres();
	void store_normals();
	void grow_tree();
	VertTree* build_tree(VertIt begin, VertIt end);

	// tree rendering
	enum NormConeDir { DIR_TOWARD, DIR_UNTOWARD, DIR_UNDECIDED };
	void draw_tree(VertTree* tree, NormConeDir cdir);
	void draw_quad(const Vec3f& pos, const Vec3f& normal, float radius);

	// other functions
	void destroy();

	// data
	Vec3f centroid;	

	// splat tree of vertices
	VertTree* splat_tree;	
	int num_nodes;
	int num_ignored;
	int num_zero_norms;

	float max_pixel_error;

	// rendering params
	ColorChanger rainbow;	
	float normal_scale;
};

#endif //WSPLAT_H