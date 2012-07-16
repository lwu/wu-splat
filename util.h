
#ifndef UTIL_H
#define UTIL_H

double get_cpu_time();

// Timer helper class
class Timer {
public:
	Timer() {
		t0 = get_cpu_time();
	}

	double get_elapsed() const {
		return get_cpu_time() - t0;
	}
private:
	double t0;	
};


// Grammar geekiness
const char* verts(int n);


// Temporarily set one object to be another
template<class T>
struct Let {
	Let(T& t, const T& let_t) : t_ref(t) {
		old_t = t;
		t = let_t;
	}

	~Let() {
		t_ref = old_t;
	}
private:
	T old_t;	
	T& t_ref;
};


#endif // UTIL_H