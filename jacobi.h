/************************************************************************

  Jacobi iteration routines for computing eigenvalues/eigenvectors.

  NOTE: This code was adapted from VTK source code (vtkMath.cxx) which
        seems to have been adapted directly from Numerical Recipes in C.

  $Id: jacobi.cxx,v 1.2 2001/09/24 21:57:17 garland Exp $

 ************************************************************************/

#include <gfx/gfx.h>
#include <gfx/mat3.h>
#include <gfx/mat4.h>

#define ROT(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau)

#define MAX_ROTATIONS 60

// Description:
// Jacobi iteration for the solution of eigenvectors/eigenvalues of a NxN
// real symmetric matrix. Square NxN matrix a; output eigenvalues in w;
// and output eigenvectors in v. Resulting eigenvalues/vectors are sorted
// in decreasing order; eigenvectors are normalized.
//
template<class Vec, class Mat>
bool __jacobi(Mat &a, Vec &w, Mat &v)
{
    const int N = Vec::dim();
    
    int i, j, k, iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c;
    double tmp;
    Vec b, z;

    // initialize
    for (ip=0; ip<N; ip++) 
    {
	for (iq=0; iq<N; iq++) v[ip][iq] = 0.0;
	v[ip][ip] = 1.0;
    }
    for (ip=0; ip<N; ip++) 
    {
	b[ip] = w[ip] = a[ip][ip];
	z[ip] = 0.0;
    }

    // begin rotation sequence
    for (i=0; i<MAX_ROTATIONS; i++) 
    {
	sm = 0.0;
	for (ip=0; ip<2; ip++) 
	{
	    for (iq=ip+1; iq<N; iq++) sm += fabs(a[ip][iq]);
	}
	if (sm == 0.0) break;

	if (i < 4) tresh = 0.2*sm/(9);
	else tresh = 0.0;

	for (ip=0; ip<2; ip++) 
	{
	    for (iq=ip+1; iq<N; iq++) 
	    {
		g = 100.0*fabs(a[ip][iq]);
		if (i > 4 && (fabs(w[ip])+g) == fabs(w[ip])
		    && (fabs(w[iq])+g) == fabs(w[iq]))
		{
		    a[ip][iq] = 0.0;
		}
		else if (fabs(a[ip][iq]) > tresh) 
		{
		    h = w[iq] - w[ip];
		    if ( (fabs(h)+g) == fabs(h)) t = (a[ip][iq]) / h;
		    else 
		    {
			theta = 0.5*h / (a[ip][iq]);
			t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
			if (theta < 0.0) t = -t;
		    }
		    c = 1.0 / sqrt(1+t*t);
		    s = t*c;
		    tau = s/(1.0+c);
		    h = t*a[ip][iq];
		    z[ip] -= h;
		    z[iq] += h;
		    w[ip] -= h;
		    w[iq] += h;
		    a[ip][iq]=0.0;
		    for (j=0;j<ip-1;j++) 
		    {
			ROT(a,j,ip,j,iq);
		    }
		    for (j=ip+1;j<iq-1;j++) 
		    {
			ROT(a,ip,j,j,iq);
		    }
		    for (j=iq+1; j<N; j++) 
		    {
			ROT(a,ip,j,iq,j);
		    }
		    for (j=0; j<N; j++) 
		    {
			ROT(v,j,ip,j,iq);
		    }
		}
	    }
	}

	for (ip=0; ip<N; ip++) 
	{
	    b[ip] += z[ip];
	    w[ip] = b[ip];
	    z[ip] = 0.0;
	}
    }

    if ( i >= MAX_ROTATIONS )
	return false;

    // sort eigenfunctions
    for (j=0; j<N; j++) 
    {
	k = j;
	tmp = w[k];
	for (i=j; i<N; i++)
	{
	    if (w[i] >= tmp) 
	    {
		k = i;
		tmp = w[k];
	    }
	}
	if (k != j) 
	{
	    w[k] = w[j];
	    w[j] = tmp;
	    for (i=0; i<N; i++) 
	    {
		tmp = v[i][j];
		v[i][j] = v[i][k];
		v[i][k] = tmp;
	    }
	}
    }

    // VTK addition to original Numerical Recipes code:
    //    insure eigenvector consistency (i.e., Jacobi can compute
    //    vectors that are negative of one another (.707,.707,0) and
    //    (-.707,-.707,0). This can wreak havoc in
    //    hyperstreamline/other stuff. We will select the most
    //    positive eigenvector.
    int numPos;
    for (j=0; j<N; j++)
    {
	for (numPos=0, i=0; i<N; i++) if ( v[i][j] >= 0.0 ) numPos++;
	if ( numPos < 2 ) for(i=0; i<N; i++) v[i][j] *= -1.0;
    }

    return true;
}

#undef ROT
#undef MAX_ROTATIONS

bool eigen(const Mat3& m, Vec3& eig_vals, Vec3 eig_vecs[3])
{
    Mat3 a, v;  Vec3 w;
    int i,j;

    for(i=0;i<3;i++) for(j=0;j<3;j++) a[i][j] = m(i,j);

    bool result = __jacobi(a, w, v);
    if( result )
    {
	for(i=0;i<3;i++) eig_vals[i] = w[i];

	for(i=0;i<3;i++) for(j=0;j<3;j++) eig_vecs[i][j] = v[j][i];
    }

    return result;
}

bool eigen(const Mat4& m, Vec4& eig_vals, Vec4 eig_vecs[4])
{
    Mat4 a, v;  Vec4 w;
    int i,j;

    for(i=0;i<4;i++) for(j=0;j<4;j++) a[i][j] = m(i,j);

    bool result = __jacobi(a, w, v);
    if( result )
    {
	for(i=0;i<4;i++) eig_vals[i] = w[i];

	for(i=0;i<4;i++) for(j=0;j<4;j++) eig_vecs[i][j] = v[j][i];
    }

    return result;
}
