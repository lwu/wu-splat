// An implementation of the QSplat multiresolution point rendering system in libgfx/OpenGL

#include <iostream>
#include <string>
#include <gfx/gui.h>
#include <gfx/raster.h>
#include <gfx/mat4.h>
#include <FL/FL_ask.H>
#include <FL/FL_file_chooser.H>
#include <FL/FL_Value_Slider.H>
#include <FL/FL_Multiline_Output.H>
#include <FL/FL_Check_Button.H>
#include <FL/FL_Counter.H>
#include <FL/FL_Color_Chooser.H>

#include "wsplat.h"
#include "util.h"

using namespace std;

class GUI : public MxGUI
{
public:
	GUI();
	~GUI();

	void init(int argc, char* argv[]);
	void file_open_dialog();
	void get_user_options();

	virtual void setup_for_drawing();
	virtual void draw_contents();
	virtual bool key_press(int key);

	virtual bool mouse_drag(int *where, int *last, int which);
	virtual bool mouse_down(int *where, int which);
	virtual bool mouse_up(int *where, int which);
	virtual void update_animation();

	virtual void add_upper_controls(int& yfill, const int pad);
	virtual void add_lower_controls(int& yfill, const int pad);

	void apply_camera();
	void update_pos();
	void reset_cam();

	void benchmark();

public:
	Fl_Value_Slider* tau_slider;
	Fl_Check_Button* btn_thresh;
	Fl_Counter*		 subdiv_level;

	Fl_Check_Button* btn_lighting;
	Fl_Check_Button* btn_normals;	
	Fl_Check_Button* btn_backface;
	Fl_Check_Button* btn_update_cam;

	int mouse_pressed;
	int mouse_x0, mouse_y0;
	int mouse_x1, mouse_y1;	

	WSplat splatter;

	Vec3f view_pos;
	Vec3f cynosure;
	Vec3f cynosure_offset;
	float view_radius;
	float max_view_radius;
	float theta, phi;
	float dtheta;

	Vec4f light_color;
	Vec4f background_color;
	bool enable_lighting;
	bool track_cam;
	bool upside_up;
};

extern int max_depth;

int max_user_depth = 0;
GUI gui;
const float MAX_FPS = 120.0f;

GUI::GUI()
{
	mouse_pressed = 0;
	view_radius = 1.0f;
	max_view_radius = 50.0f;	

	theta = phi = 0.0f;
	dtheta = 0.0f;	

	enable_lighting = true;

	update_pos();

	// widgets
	tau_slider = 0;
	btn_thresh = 0;
	subdiv_level = 0;

	btn_lighting = 0;
	btn_normals = 0;
	btn_backface = 0;	
	btn_update_cam = 0;

	light_color = Vec4f(1.0f, 1.0f, 1.0f, 1.0f);
	background_color = Vec4f(37.0f/255.0f, 34.0f/255.0f, 172.0f/255.0f, 0.0f);	

	track_cam = true;
	upside_up = true;

	view_pos = Vec3f(0.0f, 0.0f, 1.0f);
}

GUI::~GUI()
{
	delete tau_slider;
	delete btn_thresh;
	delete subdiv_level;

	delete btn_lighting;
	delete btn_normals;	
	delete btn_backface;
	delete btn_update_cam;
}

void GUI::setup_for_drawing()
{		
	glClearColor(background_color[0], background_color[1], background_color[2], 0.0f);
	glEnable(GL_DEPTH_TEST);

    // Lighting
	if (enable_lighting) {
		glEnable(GL_LIGHTING);
		Vec4f ambient_light(1.0, 1.0, 1.0, 1.0);
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, (float *)ambient_light);
		
		static bool init_light = false;
		if (!init_light) {
			const Vec4f light0_pos(0.0f, 0.5f, 1.0f, 0.0f);
			glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
			glEnable(GL_LIGHT0);
			init_light = true;
		}
		
		const Vec4f rgb4 = light_color;
		const Vec4f r_amb  = 0.1*rgb4;
		const Vec4f r_diff = 1.0*rgb4;
		const Vec4f r_spec = 0.3*rgb4;
		
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, r_amb);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, r_diff);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, r_spec);
		glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 100);
	} else {
		glDisable(GL_LIGHTING);
	}	

	static bool init = false;
	if (!init) {
		GLint psize_range[2];
		glGetIntegerv(GL_POINT_SIZE_RANGE, psize_range);
		cout << "GL_POINT_SIZE_RANGE is (" << psize_range[0]
			<< ", " << psize_range[1] << ")" << endl;
		splatter.max_point_size = psize_range[1];
		init = true;
	}
}

void GUI::reset_cam()
{
	theta = phi = 0.0f;
	cynosure_offset = Vec3f();
}

void GUI::apply_camera()
{
	float aspect = (float)canvas->w() / (float)canvas->h();

	float radius;
	gui.splatter.get_bounding_sphere(cynosure, radius);

	const float VIEW_ANGLE = 60.0f;

	const float zfar  = radius * 8.0f;
	const float znear = radius / 8.0f;

	glMatrixMode(GL_PROJECTION);
	gluPerspective(VIEW_ANGLE, aspect, znear, zfar);

	glMatrixMode(GL_MODELVIEW);
	
	Vec3f real_cynosure = cynosure + cynosure_offset;
	
	Vec3f eye = view_pos + real_cynosure;
	Vec3f at = real_cynosure;
	float up[] = { 0.0, 1.0, 0.0 };

	if (!upside_up) {
		up[1] = -1.0f;
	}

	gluLookAt(eye[0], eye[1], eye[2],
		  at[0], at[1], at[2],
		  up[0], up[1], up[2]);

	if (track_cam) {
		splatter.set_camera(eye, at - eye, canvas->w());
	}
	splatter.set_error( tau_slider->value() );
}

void GUI::draw_contents()
{
	max_user_depth = (int)subdiv_level->value();

	double t_0 = get_cpu_time();

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);  glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);   glLoadIdentity();

	apply_camera();
	
	splatter.draw();
	
	double t_1 = get_cpu_time();		
	double dt = t_1 - t_0;
	float fps = dt > 0 ? 1.0f / (float)dt : MAX_FPS;
	
	status("fps = %05.1f, %d points, max. user depth = %d/%d", fps, 
		splatter.num_points_drawn, max_user_depth, max_depth);
}

Vec3f rotate_direction(const Vec3f& moving, float angle, const Vec3f& axis)
{
	Vec3 daxis(axis[0], axis[1], axis[2]);
	Vec3 dmoving(moving[0], moving[1], moving[2]);
    Mat4 M = rotation_matrix_deg(angle, daxis);
    Vec4 v = Vec4(dmoving, 0.0);
    Vec4 vnew = M*v;

    return Vec3f(vnew[0], vnew[1], vnew[2]);
}

void GUI::update_pos() {
	
	Vec3f pos(0.0f, 0.0f, view_radius);
	Vec3f up(0.0f, 1.0f, 0.0f);

	view_pos = rotate_direction(pos, theta, up);	

	Vec3f left = cross(up, view_pos);
	left /= norm(left);
	view_pos = rotate_direction(view_pos, phi, left);	
}

bool GUI::mouse_down(int *where, int which)
{
	mouse_x0 = where[0];
	mouse_y0 = where[1];
	mouse_x1 = mouse_x0;
	mouse_y1 = mouse_y0;

	mouse_pressed = which;

	return true;
}

bool GUI::mouse_up(int *where, int which)
{
	mouse_pressed = 0;
	return true;
}

bool GUI::mouse_drag(int *where, int *last, int which)
{	
	mouse_x1 = where[0]; mouse_y1 = where[1];
	return true;
}

void bound_range(float* x, float x_low, float x_high)
{
	if (*x < x_low) {
		*x = x_low;
	} else if (*x > x_high) {
		*x = x_high;
	}
}

void GUI::update_animation()
{
	const float sensitivity = 0.02f;
	const float view_scale = gui.view_radius / 100.0f;
	const float pan_sensitivity = sensitivity / 4.0f;

	if (mouse_pressed > 0) {
		int x0 = mouse_x0, y0 = mouse_y0;
		int x1 = mouse_x1, y1 = mouse_y1;
		int dx = x1 - x0, dy = y1 - y0;
		
		if (mouse_pressed != 2) {
			dtheta = -(dx * sensitivity);			
			theta += dtheta;

			if (mouse_pressed == 1) {
				// LMB also rotates phi			
				phi -= dy * sensitivity;		
				//bound_range(&phi, -89.0f, 89.0f);

				if (phi > 180.0f) {
					phi -= 360.0f;
				} else if (phi < -180.0f) {
					phi += 360.0f;
				}

				upside_up = (fabs(phi) <= 90.0f);
			} else {
				// RMB zooms				
				view_radius += dy * sensitivity * view_scale;
				view_radius = max(0.001f, view_radius);
				view_radius = min(max_view_radius, view_radius);			
			}
		} else {
			// pan
			Vec3f dp = (splatter.cam.left * -dx +
				splatter.cam.up * -dy) * pan_sensitivity * view_scale;

			cynosure_offset += dp;
		}

		update_pos();
	}
}

bool GUI::key_press(int key)
{
	const float view_radius_delta = 0.02f;

    switch( key )
	{
		case FL_Left: 
			theta += 5.0f;
			update_pos();
			break;
        case FL_Right:
			theta -= 5.0f;
			update_pos();
			break;

        case FL_Up:     
			view_radius -= view_radius_delta;
			update_pos();
			break;
        case FL_Down:   
			view_radius += view_radius_delta; 
			update_pos();
			break;

		case ',':
			max_user_depth = max(0, max_user_depth-1);
			subdiv_level->value(max_user_depth);
			break;
		case '.':
			max_user_depth = min(max_depth, max_user_depth+1);
			subdiv_level->value(max_user_depth);
			break;
		default:
			return false;
	}

	canvas->redraw();
	return true;
}

void GUI::benchmark()
{
	const int N = 20;

	theta = phi = 0.0f;
	float dtheta = 360.0f / (float)N;	

	Vec3f centre;
	float radius;
	splatter.get_bounding_sphere(centre, radius);

	Let<bool> L0(track_cam, true);
	Let<bool> L1(splatter.draw_normals, false);	
	Let<float> L2(view_radius, radius);
	Let<Vec3f> L3(cynosure, centre);	
	Let<Vec3f> L4(cynosure_offset, Vec3f());
	
	Timer t;

	for (int i=0; i < N; i++) {		
		draw_contents();		
		theta += dtheta;		
	}

	double dt = t.get_elapsed();
	float fps = (dt==0) ? MAX_FPS : (float)N / dt;

	fl_message("Timedemo: %2.1f average fps", fps);
}

void redraw()
{	
	gui.setup_for_drawing();
	gui.canvas->redraw();
}

static void cb_redraw(Fl_Menu_ *, void *)
{
	redraw();
}

static void cb_thresh(Fl_Menu_ *, void *)
{
	gui.splatter.always_subdivide = !gui.splatter.always_subdivide;
	redraw();
}

static void cb_toggle_lighting(Fl_Menu_ *, void *)
{
	gui.enable_lighting = !gui.enable_lighting;
	redraw();
}

static void cb_toggle_normals(Fl_Menu_ *, void *)
{
	gui.splatter.draw_normals = !gui.splatter.draw_normals;
	redraw();
}

static void cb_toggle_backface(Fl_Menu_ *, void *)
{
	gui.splatter.backface_culling = !gui.splatter.backface_culling;
	redraw();
}

static void cb_toggle_update_cam(Fl_Menu_ *, void *)
{
	gui.track_cam = !gui.track_cam;
}


void GUI::add_upper_controls(int& yfill, const int pad) 
{
	int size_x = 280;
	int size_y = 16;

	// Add various sliders
	yfill += pad;
	int x = 20;
	tau_slider = new Fl_Value_Slider(20,yfill, size_x, size_y, "Error (pixels)");
	tau_slider->type(FL_HOR_NICE_SLIDER);
	tau_slider->minimum(1.0f);
	tau_slider->maximum(256.0f);	
	tau_slider->value(16.0f);
	x += size_x;

	int thresh_width = 120;
	x += 0;
	btn_thresh = new Fl_Check_Button(x, yfill, thresh_width, 20, "Enabled");
	btn_thresh->set();
	btn_thresh->callback( (Fl_Callback*)cb_thresh );
	x += thresh_width;

	int level_width = 120;
	x += 20;
	subdiv_level = new Fl_Counter(x, yfill, level_width, 20, "Max depth");	
	subdiv_level->range(0,32);
	subdiv_level->step(1.0, 4.0);
	subdiv_level->callback( (Fl_Callback*)cb_redraw );

	int tau_yfill = tau_slider->h() + tau_slider->labelsize();
	yfill += tau_yfill;
}


void GUI::add_lower_controls(int& yfill, const int pad) 
{
	yfill += pad;

	// check buttons
	int x = 0;
	int light_width = 80;	
	x += 20;
	btn_lighting = new Fl_Check_Button(x, yfill, light_width, 20, "Lighting");
	btn_lighting->value( gui.enable_lighting );
	btn_lighting->callback( (Fl_Callback*)cb_toggle_lighting );
	x += light_width;

	int normals_width = 80;
	x += 20;
	btn_normals = new Fl_Check_Button(x, yfill, light_width, 20, "Normals");
	btn_normals->value( gui.splatter.draw_normals );
	btn_normals->callback( (Fl_Callback*)cb_toggle_normals );
	x += normals_width;

	int cam_width = 120;
	x += 20;
	btn_update_cam = new Fl_Check_Button(x, yfill, cam_width, 20, "Track Camera");
	btn_update_cam->value( gui.track_cam );
	btn_update_cam->callback( (Fl_Callback*)cb_toggle_update_cam );
	x += cam_width;

	int backface_width = 120;
	x += 20;
	btn_backface = new Fl_Check_Button(x, yfill, light_width, 20, "Backface Culling");
	btn_backface->value( gui.splatter.backface_culling );
	btn_backface->callback( (Fl_Callback*)cb_toggle_backface );
	x += backface_width;

	yfill += btn_lighting->h();	
}

const char* choose_file(const char* msg, const char* pattern, const char* default_fname)
{
	//return default_fname;
	return fl_file_chooser(msg, pattern, default_fname);
}

void sleep(int n) {
	float t0 = get_cpu_time();

	while (get_cpu_time() - t0 < n) {}
}

class Popup {
public:
	Popup(const char* window_label) :
		wnd(width, height, window_label),
		output(1,1, width-1, height-1, "")
	{
			wnd.set_modal();
			wnd.add(output);
	}

	void text(const char* txt) 
	{
		output.value(txt);
		wnd.show();
		Fl::wait();
	}
private:
	Fl_Window wnd;
	Fl_Multiline_Output output;
	static int width;
	static int height;
};

int Popup::width = 384;
int Popup::height = 48;


static void cb_partition(Fl_Menu_ *, void *)
{
	gui.splatter.oriented_partition = !gui.splatter.oriented_partition;
}

void GUI::get_user_options()
{
	int rv = fl_choice("Partitioning method:", 
		"Oriented", "Axis-Aligned", NULL);

	splatter.oriented_partition = (rv == 0);

	rv = fl_ask("Recursively tessellate as needed?");

	splatter.subdivide_large_tris = (rv == 1);

	if (splatter.subdivide_large_tris) {
		rv = fl_choice("Subdivision level:",
			"A little", "A lot", NULL);
		splatter.subdivide_a_lot = (rv == 1);
	}
}

void GUI::file_open_dialog() {
	const char *fname = choose_file("Model", "*.smf", "bunny.smf");

	if (!fname) return;

	get_user_options();

	Popup p("Status");
	string status = string("Loading ") + fname;
	p.text(status.c_str());

	gui.splatter.load(fname);

	// initially at max. depth
	max_user_depth = max_depth;
	subdiv_level->range(0, max_user_depth);
	subdiv_level->value(max_depth);

	// get largest bounding sphere
	const float view_dist_mult = 1.3f;
	gui.splatter.get_bounding_sphere(cynosure, view_radius);
	
	max_view_radius = view_radius * 10.0f;
	view_radius *= view_dist_mult;

	gui.update_pos();
}

int main(int argc, char *argv[])
{		
	gui.init(argc, argv);
	gui.file_open_dialog();

	return gui.run();
}


// GUI data and functions

static void cb_load(Fl_Menu_ *, void *)
{
	gui.file_open_dialog();
	gui.canvas->redraw(); // force redraw
}

static void cb_run_timedemo(Fl_Menu_*, void *)
{
	gui.benchmark();
}

static void cb_reset_cam(Fl_Menu_*, void *)
{
	gui.reset_cam();
	gui.update_pos();
	redraw();
}

static void cb_exit(Fl_Menu_ *, void *) // copied from gui.cxx:120
	{ MxGUI::current->cleanup_for_exit(); exit(0); }

static void cb_snapshot(Fl_Menu_ *m, int format) // copied from gui.cxx:141
{
    MxGUI::current->canvas->redraw();		// don't want to snap menu
    MxGUI::current->snapshot_to_file(format);	// snapshot what's drawn
}

static void cb_display_size(Fl_Menu_ *m, int xw) // copied from gui.cxx:147
{
    if( MxGUI::current->toplevel->resizable() )
    {
	int yw = (3*xw)/4;
	MxGUI::current->resize_canvas(xw, yw);
    }
}

static void cb_animate(Fl_Menu_ *m, void *) // copied from gui.cxx:123
    { MxGUI::current->animate(m->mvalue()->value()!=0); }

float get_input(const char* prompt, float orig)
{
	MxGUI *gui = MxGUI::current;

    // Convert default to a string
    static char def[64]; sprintf(def, "%.1f", orig);

    const char *result = fl_input(prompt, def);
	float rv = 1.0f;
    if( result ) {
		rv = atof(result);
    }

	return rv;
}

static void cb_fps(Fl_Menu_ *m, void *) // copied from gui.cxx:126
{
    MxGUI *gui = MxGUI::current;
	float fps = get_input("Number of frames per second to draw", gui->default_fps);
	if( gui->target_fps>0 ) gui->target_fps=fps;
}

static void cb_silhouette(Fl_Menu_ *m, void *)
{
	float thresh = get_input("Silhouette edge threshold [0.0, 1.0f]", gui.splatter.silhouette_edge_thresh);

	thresh = max(0.0f, thresh);
	thresh = min(1.0f, thresh);

	gui.splatter.silhouette_edge_thresh = thresh;
}


static void cb_about(Fl_Menu_ *m, void *)
{
	fl_message("\tWuSplat\n\tby Leslie Wu");
}

void pick_color(const char* title, Vec4f& v)
{
	uchar r = (uchar)(v[0]*255.0f);
	uchar g = (uchar)(v[1]*255.0f);
	uchar b = (uchar)(v[2]*255.0f);	
	if (fl_color_chooser(title, r,g,b)) {
		v[0] = (float)r/255.0f;
		v[1] = (float)g/255.0f;
		v[2] = (float)b/255.0f;		
	}
}

static void cb_choose_color(Fl_Menu_ *m, void *)
{
	pick_color("Light color", gui.light_color);
	gui.setup_for_drawing();
}

static void cb_choose_background(Fl_Menu_ *m, void *)
{
	pick_color("Background color", gui.background_color);
	gui.setup_for_drawing();
}


void GUI::init(int argc, char* argv[])
{	
	// Setup menu
	Fl_Menu_Item menu_layout[] =
	{
		MXGUI_BEGIN_MENU("&File")
		{"L&oad", FL_CTRL + 'o', (Fl_Callback *)cb_load, NULL, FL_MENU_DIVIDER},
        MXGUI_BEGIN_MENU("Sna&pshot to")
#if defined(HAVE_LIBPNG)
        {"&PNG", FL_CTRL+'p', (Fl_Callback *)cb_snapshot, (void *)IMG_PNG},
#else
        {"&PNG", FL_CTRL+'p', NULL, NULL, FL_MENU_INACTIVE },
#endif
#if defined(HAVE_LIBTIFF)
        {"&TIFF", FL_CTRL+'P', (Fl_Callback *)cb_snapshot, (void *)IMG_TIFF},
#else
        {"&TIFF", FL_CTRL+'P', NULL, NULL, FL_MENU_INACTIVE},
#endif
        {"PP&M", 0, (Fl_Callback *)cb_snapshot, (void *)IMG_PNM},
#if defined(HAVE_LIBJPEG)
        {"&JPEG", 0, (Fl_Callback *)cb_snapshot, (void *)IMG_JPEG},
#else
        {"&JPEG", 0, NULL, NULL, FL_MENU_INACTIVE},
#endif
        MXGUI_END_MENU
		{"E&xit", FL_CTRL + 'q',  (Fl_Callback *)cb_exit, NULL},
		MXGUI_END_MENU
		
		MXGUI_BEGIN_MENU("&View")
		{"&Animate", FL_CTRL+'a', (Fl_Callback *)cb_animate, NULL, FL_MENU_TOGGLE},
		{"Animation speed ...", FL_CTRL+'r', (Fl_Callback *)cb_fps, NULL},		
		
		MXGUI_BEGIN_MENU("Display &size")
		{"&320x240", 0, (Fl_Callback *)cb_display_size, (void*)320},
		{"&640x480", 0, (Fl_Callback *)cb_display_size, (void*)640},
		{"&800x600", 0, (Fl_Callback *)cb_display_size, (void*)800},
		{"&1024x768", 0, (Fl_Callback *)cb_display_size, (void*)1024},
		MXGUI_END_MENU
		MXGUI_END_MENU

		MXGUI_BEGIN_MENU("&Rendering")
		{"Choose light color", FL_CTRL+'c', (Fl_Callback *)cb_choose_color, NULL},
		{"Choose background color", FL_CTRL+'b', (Fl_Callback *)cb_choose_background, NULL},
		{"Silhouette threshold ...", FL_CTRL+'h', (Fl_Callback *)cb_silhouette, NULL},		
		MXGUI_END_MENU

		MXGUI_BEGIN_MENU("&Options")
		{"Run timedemo", FL_CTRL+'t', (Fl_Callback *)cb_run_timedemo, NULL},
		{"Reset viewer", FL_CTRL+'r', (Fl_Callback *)cb_reset_cam, NULL},
		MXGUI_END_MENU
		
		MXGUI_FINISH_MENUBAR
	};

	gui.initialize(argc, argv, &menu_layout[0]);

	toplevel->label("WuSplat");

	// Add menu toggles
	add_toggle_menu("&Rendering/Silhouette edges", FL_CTRL+'y',
			splatter.draw_edges_only);
    add_toggle_menu("&Rendering/Draw Splats", FL_CTRL+'d',
			splatter.draw_splats);
    //add_toggle_menu("&Rendering/Draw Triangles", FL_CTRL+'t',
	//		splatter.draw_triangles);
    add_toggle_menu("&Rendering/Draw Points", FL_CTRL+'p',
			splatter.draw_points);
    add_toggle_menu("&Rendering/Draw Spheres", FL_CTRL+'s',
			splatter.draw_spheres);

	add_toggle_menu("&Options/Project Radius", FL_CTRL+'p',
			splatter.project_size);
	add_toggle_menu("&Options/Allow Quads", FL_CTRL+'u',
			splatter.allow_quads);
	add_toggle_menu("&Options/Quads only", FL_CTRL+'l',
			splatter.quads_only);
	add_toggle_menu("&Options/Abnormal", FL_CTRL+'b',
			splatter.abnormal);
	
	add_toggle_menu("&OpenGL/Round Points", FL_CTRL+'r',
			splatter.round_points);
	add_toggle_menu("&Visualization/Draw Cones", FL_CTRL+'e',
			splatter.draw_cones);
	add_toggle_menu("&Visualization/Modulate Colors", FL_CTRL+'m',
			splatter.change_colors);	
	add_toggle_menu("&Fun/Psychedelic", FL_CTRL+'y',
			splatter.psychedelic);

	menu_bar->add("&Help/About", 0, (Fl_Callback*)cb_about);
}
