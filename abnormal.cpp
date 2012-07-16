
#include <assert.h>
#include <iostream>
#include "abnormal.h"

using namespace std;

float error(const Vec3f& v0, const Vec3f& v1)
{
	float cos_theta = v0*v1;
	return acos(cos_theta) * (180.0f / M_PI);
}

void explicit_test(const Vec3f& v)
{
	Norm in = v;
	Vec3f out = in;

	cout << v << " -> " << out << ",\n\terror = "
		<< error(v, out) << " degrees." << endl;
}

void test_abnormals()
{
	Norm::generate_table();

	explicit_test( Vec3f(1.0f, 0.0f, 0.0f) );
	explicit_test( Vec3f(0.0f, 1.0f, 0.0f) );
	explicit_test( Vec3f(0.0f, 0.0f, 1.0f) );

	const int N = 8;
	double error_sum = 0;

	for (int i=0; i < N; i++) {
		float x = (rand() / (double)RAND_MAX) * 0.5f - 0.25f;
		float y = (rand() / (double)RAND_MAX) * 0.5f - 0.25f;
		float z = (rand() / (double)RAND_MAX) * 0.5f - 0.25f;
		Vec3f v( x,y,z );
		v /= norm(v);

		Norm in = v;
		Vec3f out = in;

		explicit_test(v);

		error_sum += error(v,out);
	}

	error_sum /= (double)N;

	cout << "(" << BIT_DIM*2 << " bits per normal.)" << endl;
	cout << "Average error over " << N << " random normals = "
		<< error_sum << " degrees." << endl;
}


const float bit_scale = (float)(BIT_DIM-1);

Vec3f Norm::vec_table[BIT_DIM][BIT_DIM];

void Norm::generate_table()
{
	const float bit_inv = 1.0f / bit_scale;

	std::cout << "Generating norm table..." << std::endl;
	for (int x=0; x < BIT_DIM; x++) {
		for (int y=0; y < BIT_DIM; y++) {
			if (x == 13 && y == 0) {
				bool boob = true;
			}

			if (x == 6 && y == 2) {
				bool boob = true;
			}

			float fx = (float)x * bit_inv;
			float fy = (float)y * bit_inv;
			
			float z_sqr = 1.0f - fx*fx - fy*fy;
			
			float fz = 0.0f;
			if (z_sqr >= 0.0f) {
				fz = sqrt(z_sqr);
			}
				
			vec_table[x][y] = Vec3f(fx, fy, fz);

			if (z_sqr < 0.0f) {
				vec_table[x][y] /= norm(vec_table[x][y]);
			}

		}
	}
}

Norm::Norm()
{
	self.index_x = self.index_y = 0;
}

Norm::Norm(const Vec3f& v) 
{
	operator=(v);
}

void Norm::operator=(const Vec3f& v)
{
	assert(norm2(v) <= 1.001f);
	float n2 = norm2(v);	
	
	float x, y;

	if (v[0] >= 0.0f) {
		self.x = 1;
		x = v[0];
	} else {
		self.x = 0;
		x = -v[0];
	}

	if (v[1] >= 0.0f) {
		self.y = 1;
		y = v[1];
	} else {
		self.y = 0;
		y = -v[1];
	}
	
	self.z = (v[2] >= 0.0f);

	self.index_x = (index_type)(x*bit_scale+0.5f);
	self.index_y = (index_type)(y*bit_scale+0.5f);	

	assert(self.index_x >= 0 && self.index_x < BIT_DIM);
	assert(self.index_y >= 0 && self.index_y < BIT_DIM);
}

Norm::operator Vec3f()
{
	Vec3f v = vec_table[self.index_x][self.index_y];
	if (!self.x) {
		v[0] = -v[0];
	}

	if (!self.y) {
		v[1] = -v[1];
	}

	if (!self.z) {
		v[2] = -v[2];
	}

	return v;
}

void Norm::operator+=(const Vec3f& v) 
{
	Vec3f new_v = this->operator Vec3f();
	new_v += v;

	this->operator=(new_v);
}