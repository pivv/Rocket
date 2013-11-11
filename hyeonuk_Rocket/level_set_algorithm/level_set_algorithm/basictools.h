#include <iostream>
#include <cstdlib>
#include <time.h>
#include <cmath>


#define CHECK_TIME_START __int64 freq, start, end; if (QueryPerformanceFrequency((_LARGE_INTEGER*)&freq)) {QueryPerformanceCounter((_LARGE_INTEGER*)&start);
// a는 float type milli second이고 b가 FALSE일때는 에러입니다
#define CHECK_TIME_END(a,b) QueryPerformanceCounter((_LARGE_INTEGER*)&end); a=(float)((double)(end - start)/freq*1000); b=TRUE; } else b=FALSE;

#define MAXSIZE 1
#define MINERROR 1e-7

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define PI 3.141592653589793238


inline double sign(double r)
{
	if(r>0) return 1;
	else if(r==0) return 0;
	else return -1;
}

inline double minmod(double a, double b)
{
	return (sign(a)+sign(b))/2 * min(abs(a), abs(b));
}


double distance_circle_point(double x, double y, double r0, double x0, double y0)
{
    double r = sqrt(pow(x-x0,2) + pow(y-y0,2));
    return r-r0;
}

double vec_inner_product(int dim, double *v1, double *v2, double *w1, double *w2)
{
	double r = 0;
	for(int i=0; i<dim; i++) r += (v2[i]-v1[i])*(w2[i]-w1[i]);
	return r;
}

double *vec_normalization(int dim, double *v)
{
	double r = 0;
	for(int i=0; i<dim; i++) r += v[i]*v[i];
	for(int i=0; i<dim; i++) v[i] = v[i]/r;
	return v;
}

inline double vec_rot(double *v1, double *v2, double *w1, double *w2) {return (v2[0]-v1[0])*(w2[1]-w1[1]) - (v2[1]-v1[1])*(w2[0]-w1[0]);} // dim = 2

double *vec_curl(double *v1, double *v2, double *w1, double *w2) // dim = 3
{
	double *r = new double[3];
	r[0] = (v2[1]-v1[1])*(w2[2]-w1[2]) - (v2[2]-v1[2])*(w2[1]-w1[1]);
	r[1] = (v2[2]-v1[2])*(w2[0]-w1[0]) - (v2[0]-v1[0])*(w2[2]-w1[2]);
	r[2] = (v2[0]-v1[0])*(w2[1]-w1[1]) - (v2[1]-v1[1])*(w2[0]-w1[0]);
	return r;
}

double distance_face_point(int dim, double *v, double *w1, double *w2, double *w3)
{
	if (dim==2)
	{		
		double d = abs((w2[1]-w1[1])*(v[0]-w1[0]) - (w2[0]-w1[0])*(v[1]-w1[1]))/sqrt(pow(w2[0]-w1[0],2) + pow(w2[1]-w1[1],2));
		return d;
	}
	if (dim==3)
	{

		if(w3!=NULL)
		{
			double a = (w2[1]-w1[1])*(w3[2]-w1[2]) - (w2[2]-w1[2])*(w3[1]-w1[1]);
			double b = (w2[2]-w1[2])*(w3[0]-w1[0]) - (w2[0]-w1[0])*(w3[2]-w1[2]);
			double c = (w2[0]-w1[0])*(w3[1]-w1[1]) - (w2[1]-w1[1])*(w3[0]-w1[0]);
			double d = abs(a*(v[0]-w1[0]) + b*(v[1]-w1[1]) + c*(v[2]-w1[2]))/sqrt(pow(a,2)+pow(b,2)+pow(c,2));
			return d;
		}
		else
		{
			double l2 = pow(w2[0]-w1[0],2) + pow(w2[1]-w1[1],2) + pow(w2[2]-w1[2],2);
			double inner = (w1[0]-v[0])*(w2[0]-w1[0]) + (w1[1]-v[1])*(w2[1]-w1[1]) + (w1[2]-v[2])*(w2[2]-w1[2]);
			double a = w1[0]-v[0] - inner/l2 * (w2[0]-w1[0]);
			double b = w1[1]-v[1] - inner/l2 * (w2[1]-w1[1]);
			double c = w1[2]-v[2] - inner/l2 * (w2[2]-w1[2]);
			double d = sqrt(pow(a,2)+pow(b,2)+pow(c,2));
			return d;
		}
	}
	return 0;
}