#include <stdio.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include "mex.h"
#include <ctime>
#include "matrix.h"
#include <algorithm>
#include <time.h>
#include <vector>
#include <conio.h>

#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
//#define G(i,j,r,c) ((i)*(c)+(j))
#define G(i,j,r,c) ((i)+(j)*(r))
//#define H(i,j,k, r,c,d) ((i)*(c)*(d)+(j)*(d)+(k))
//#define H(i,j,k, r,c,d) ((i)*(d)+(j)*(r*d)+(k))
#define H(i,j,k, r,c,d) ((i)+(j)*(r)+(k)*(r)*(c))
//#define GG(i,j)		((i)*(C)+(j))
#define GG(i,j)		((i)+(j)*(R))
//#define HH(i,j,k,d)	((i)*(C)*(d)+(j)*(d)+(k))
//#define HH(i,j,k,d)	((i)*(d)+(j)*(R)*(d)+(k))
#define HH(i,j,k,d)		((i) + (j)*(R) + (k)*R*C)
#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180))
#define sq(x) ((x)*(x))



class fij {
private:
	int data[2];
public:
	fij();
	fij(int*);
	fij(int, int);
	~fij();
	int norm2();
	int& operator [](int);
	friend fij operator +(fij&, fij&);
	friend fij operator -(fij&, fij&);
	friend int operator *(fij&, fij&);
	friend fij operator *(fij&, int&);
};



fij::fij() {
	data[0] = 0;
	data[1] = 0;
}
fij::fij(int* d) {
	data[0] = d[0];
	data[1] = d[1];
}
fij::fij(int a, int b) {
	data[0] = a;
	data[1] = b;
}
fij::~fij() {}
int fij::norm2() {
	return data[0] * data[0] + data[1] * data[1];
}
int& fij::operator [](int idx) {
	if (idx > 1 || idx < 0) {
		printf("index exceeds fij length.\n");
		return data[0];
	}
	return data[idx];
}
//fij operator +(fij& b) {
//	fij ret;
//	ret[0] = data[0] + b[0];
//	ret[1] = data[1] + b[1];
//	return ret;
//}
fij operator +(fij& f1, fij& f2) {
	fij ret;
	ret[0] = f1[0] + f2[0];
	ret[1] = f1[1] + f2[1];
	return ret;
}
fij operator -(fij& f1, fij& f2) {
	fij ret;
	ret[0] = f1[0] - f2[0];
	ret[1] = f1[1] - f2[1];
	return ret;
}
int operator *(fij& f1, fij& f2) {
	return f1[0] * f2[0] + f1[1] * f2[1];
}
fij operator *(fij& f, int & a) {
	fij ret(f[0] * a, f[1] * a);
	return ret;
}