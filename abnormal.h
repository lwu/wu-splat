
#ifndef ABNORMAL_H
#define ABNORMAL_H

#include <gfx/vec3.h>

#define BIT_DIM 14

// uncomment to enable abnormal quantization
#define __ABNORMAL

void test_abnormals();

class Norm {
public:
	Norm();
	Norm(const Vec3f& v);

	void operator+=(const Vec3f& v);
	void operator=(const Vec3f& v);

	operator Vec3f();

	static void generate_table();

	typedef unsigned index_type;
	struct data {
		index_type index_x : BIT_DIM;
		index_type index_y : BIT_DIM;
		unsigned x : 1;
		unsigned y : 1;
		unsigned z : 1;
	};
	
	data self;

private:


	static Vec3f vec_table[BIT_DIM][BIT_DIM];
};

#endif // ABNORMAL_H