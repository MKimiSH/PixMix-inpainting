#include <stdio.h>
#include <math.h>
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
// Function Declarations
void fillOneLevel(int* initf, float* I, const bool* M, float* D, int level, int useline, int);

class fij;
void fillIwithF(); 
double calcSpatialCost(int i, int j, fij ff);
double calcAppearanceCost(int i, int j, fij ff);
inline fij  getFij(int i, int j);
inline void setFij(int i, int j, fij ff);
inline bool isMask(int i, int j);
inline bool isValid(int i, int j);
inline bool isValid(fij ff);
inline fij	randFij(int i, int j, double dist);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// global variables
int		*f;
float	*im;
bool	*mask;
int		R=0, C=0;
//#define GG(i,j)		((i)*(C)+(j))
#define GG(i,j)		((i)+(j)*(R))
//#define HH(i,j,k,d)	((i)*(C)*(d)+(j)*(d)+(k))
//#define HH(i,j,k,d)	((i)*(d)+(j)*(R)*(d)+(k))
#define HH(i,j,k,d)		((i) + (j)*(R) + (k)*R*C)

class fij {
private:
	int data[2];
public:
	fij() {
		data[0] = 0;
		data[1] = 0;
	}
	fij(int* d) {
		data[0] = d[0];
		data[1] = d[1];
	}
	fij(int a, int b) {
		data[0] = a;
		data[1] = b;
	}
	~fij(){}
	int norm2() {
		return data[0] * data[0] + data[1] * data[1];
	}
	int& operator [](int idx) {
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
	friend fij operator +(fij& f1, fij& f2) {
		fij ret;
		ret[0] = f1[0] + f2[0];
		ret[1] = f1[1] + f2[1];
		return ret;
	}
	friend fij operator -(fij& f1, fij& f2) {
		fij ret;
		ret[0] = f1[0] - f2[0];
		ret[1] = f1[1] - f2[1];
		return ret;
	}
};



// main procedure
void fillOneLevel(int* initf, float* I, const bool* M, float* D, int level, int useline, int iternum) {
	printf("Entering fillOneLevel\n");
	double alphaSp = .0005, alphaAp = 1 - alphaSp;
	double wrs = .5;
	int rrs = MIN(R, C);
	int i, j, k, r = R, c = C;

	// fill the image with initial guess of f
	//fillIwithF();

	// propagation and random search
	int numIter = iternum;
	for(int it = 0; it<numIter; ++it) {
		fillIwithF();

		//printf("begin propagation and search.\n");
		if (it % 2 == 0) {
			fij fijOld, fijLeft, fijUp, fijRS;
			double costOld, costLeft, costUp;
			double spcOld, spcLeft, spcUp;
			double apcOld, apcLeft, apcUp;
			double drs, spcRS, apcRS, costRS;
			for (i = 0; i < R; ++i) {
				for (j = 0; j < C; ++j) {
					if (isMask(i,j)) {
						// propagate forward
						fijOld = getFij(i, j);
						costLeft = 1e10;
						costUp = 1e10;
						spcOld = calcSpatialCost(i, j, fijOld);
						apcOld = calcAppearanceCost(i, j, fijOld);
						costOld = alphaAp*apcOld + alphaSp*spcOld;
						fijLeft = getFij(i, j - 1) + fij(0,1);
						fijUp = getFij(i - 1, j) + fij(1, 0);
						if (isMask(i, j - 1) && isValid(fijLeft))
						{
							spcLeft = calcSpatialCost(i, j, fijLeft);
							apcLeft = calcAppearanceCost(i, j, fijLeft);
							costLeft = spcLeft*alphaSp + apcLeft*alphaAp;
						}
						if (isMask(i - 1, j) && isValid(fijUp))
						{
							spcUp = calcSpatialCost(i, j, fijUp);
							apcUp = calcAppearanceCost(i, j, fijUp);
							costUp = spcUp*alphaSp + apcUp*alphaAp;
						}
						if (costUp < costOld && costUp < costLeft) 
						{
							costOld = costUp;
							setFij(i, j, fijUp);
						}
						else if (costLeft < costOld) {
							costOld = costLeft;
							setFij(i, j, fijLeft);
						}
						//random search
						drs = wrs*rrs;
						srand((unsigned)time(NULL));
						
						while (drs>D[GG(i, j)] + 3) {
							fijRS = randFij(i, j, drs);
							costRS = 1e10;
							int numTol = 3;
							while (!isValid(fijRS) && numTol > 0) {
								fijRS = randFij(i, j, drs);
								numTol--;
							}
							if (numTol <= 0) {
								continue;
							}
							spcRS = calcSpatialCost(i, j, fijRS);
							apcRS = calcAppearanceCost(i, j, fijRS);
							costRS = spcRS*alphaSp + apcRS*alphaAp;
							if (costRS < costOld) {
								setFij(i, j, fijRS);
								costOld = costRS;
							}
							drs *= wrs;
						}

					}
				}
			}

		}
		else {
			fij fijOld, fijRight, fijDown, fijRS;
			double costOld, costRight, costDown;
			double spcOld, spcRight, spcDown;
			double apcOld, apcRight, apcDown;
			double drs, spcRS, apcRS, costRS;
			for (j = C - 1; j >= 0; --j) {
				for (i = R - 1; i >= 0; --i) {
					if (isMask(i, j)) {
						// propagate forward
						fijOld = getFij(i, j);
						costRight = 1e10;
						costDown = 1e10;
						spcOld = calcSpatialCost(i, j, fijOld);
						apcOld = calcAppearanceCost(i, j, fijOld);
						costOld = alphaAp*apcOld + alphaSp*spcOld;
						fijRight = getFij(i, j + 1) - fij(0, 1);
						fijDown = getFij(i + 1, j) - fij(1, 0);
						if (isMask(i, j + 1) && isValid(fijRight))
						{
							spcRight = calcSpatialCost(i, j, fijRight);
							apcRight = calcAppearanceCost(i, j, fijRight);
							costRight = spcRight*alphaSp + apcRight*alphaAp;
						}
						if (isMask(i + 1, j) && isValid(fijDown))
						{
							spcDown = calcSpatialCost(i, j, fijDown);
							apcDown = calcAppearanceCost(i, j, fijDown);
							costDown = spcDown*alphaSp + apcDown*alphaAp;
						}
						if (costDown < costOld && costDown < costRight)
						{
							costOld = costDown;
							setFij(i, j, fijDown);
						}
						else if (costRight < costOld) {
							costOld = costRight;
							setFij(i, j, fijRight);
						}
						//random search
						drs = wrs*rrs;
						srand((unsigned)time(NULL));

						while (drs>D[GG(i, j)] + 3) {
							fijRS = randFij(i, j, drs);
							costRS = 1e10;
							int numTol = 3;
							while (!isValid(fijRS) && numTol > 0) {
								fijRS = randFij(i, j, drs);
								numTol--;
							}
							if (numTol <= 0) {
								continue;
							}
							spcRS = calcSpatialCost(i, j, fijRS);
							apcRS = calcAppearanceCost(i, j, fijRS);
							costRS = spcRS*alphaSp + apcRS*alphaAp;
							if (costRS < costOld) {
								setFij(i, j, fijRS);
								costOld = costRS;
							}
							drs *= wrs;
						}

					}
				}
			}
		}
		//printf("end propagation and search.\n");
	}



}

void fillIwithF() {
	//printf("entering fillIwithF\n");
	int i, j, k, r = R, c = C;
	//for (i = 0; i < R; ++i) {
	//	for (j = 0; j < C; ++j) {
	//		printf("%d ", f[HH(i,j,0,2)]);
	//	}
	//	printf("\n");
	//}

	for (i = 0; i < R; ++i) {
		for (j = 0; j < C; ++j) {
			//if (mask[G(i, j, r, c)])
			if(isMask(i,j))
			{
				//int fij[2] = { f[H(i,j,0,r,c,2)], f[H(i,j,1,r,c,2)] };
				fij ff = getFij(i, j);
				//printf("ff, %d, %d\n", ff[0], ff[1]);
				for (k = 0; k < 3; ++k) {
					im[H(i, j, k, r, c, 3)] = im[H(ff[0], ff[1], k, r, c, 3)];
				}
			}
		}
	}
}

// calculate cost_spatial using 8-neighbor
double calcSpatialCost(int i, int j, fij ff){

	double spc = 0;
	double w = 0.125;
	double tao = 200;
	int r_begin = MAX(i - 1, 0), r_end = MIN(i + 1, R-1);
	int c_begin = MAX(j - 1, 0), c_end = MIN(j + 1, C-1);
	int r, c;
	for (r = r_begin; r <= r_end; ++r) {
		for (c = c_begin; c <= c_end; ++c) {
			//% f(p) = fij, f(p + v) = f(r, c, :), v = [r - i], c - j], f(p) + v = fij + v
			fij v = fij(r, c) - fij(i, j), fv = getFij(r,c);
			fij diff = fv - (ff + v);
			spc = spc + MIN(diff.norm2(), tao);
		}
	}
	spc = spc*w;
	return spc;
}

//calculate cost_appear using 5x5 patch
double calcAppearanceCost(int i, int j, fij ff) {

	double apc = 0;
	double w = 1.0 / 25;
	int r_begin = MAX(i - 2, 0), r_end = MIN(i + 2, R - 1);
	int c_begin = MAX(j - 2, 0), c_end = MIN(j + 2, C - 1);
	int r, c;
	float i1[3] = {}, i2[3] = {};
	for (r = r_begin; r <= r_end; ++r) {
		for (c = c_begin; c <= c_end; ++c) {
			fij v(r-i, c-j), fv = ff+v;
			for (int k = 0; k < 3; ++k) {
				i1[k] = im[HH(r, c, k, 3)];
			}
			if (isValid(fv)) {
				for (int k = 0; k < 3; ++k) {
					i2[k] = im[HH(fv[0], fv[1], k, 3)];
				}
			}
			apc = apc + (i1[0] - i2[0])*(i1[0] - i2[0]) +
				(i1[1] - i2[1])*(i1[1] - i2[1]) + (i1[2] - i2[2])*(i1[2] - i2[2]);
		}
	}
	apc = apc*w;
	return apc;
}

inline fij  getFij(int i, int j) {
	fij ret;
	ret[0] = f[HH(i, j, 0, 2)];
	ret[1] = f[HH(i, j, 1, 2)];
	return ret;
}
inline void setFij(int i, int j, fij ff) {
	f[HH(i, j, 0, 2)] = ff[0];
	f[HH(i, j, 1, 2)] = ff[1];
}
inline bool isMask(int i, int j) {
	return (i >= 0 && i < R) && (j >= 0 && j < C) && mask[GG(i, j)] != 0;
}
inline bool isValid(int i, int j) {
	return (i >= 0 && i < R) && (j >= 0 && j < C) && mask[GG(i, j)] == 0;
}
inline bool isValid(fij ff) {
	return isValid(ff[0], ff[1]);
}
//产生随机搜索的映射
inline fij	randFij(int i, int j, double dist) {
	fij ret;
	ret[0] = (int)(((rand() / (double)RAND_MAX) - 0.5) * 2 * dist) + i;
	ret[1] = (int)(((rand() / (double)RAND_MAX) - 0.5) * 2 * dist) + j;
	return ret;
}

//=====================================================================
//mexFunction Entry point
//=====================================================================

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#define FILLEDI		plhs[0]
#define RETF		plhs[1]

#define INITF		prhs[0]
#define	INI			prhs[1]
#define INM			prhs[2]
#define	IND			prhs[3]
#define LEVEL		prhs[4]
#define USELINE		prhs[5]
#define ITERNUM		prhs[6]

	// input
	int*	initf = (int*)mxGetData(INITF);
	float*	I =	(float*)mxGetData(INI);
	bool*	M = (bool*)mxGetData(INM);
	float*	D = (float*)mxGetData(IND);
	
	int		level = (int)mxGetScalar(LEVEL);
	int		useline = (int)mxGetScalar(USELINE);
	int		iternum = (int)mxGetScalar(ITERNUM);

	// output
	R = mxGetM(INM);
	C = mxGetN(INM);
	printf("R=%d, C=%d\n", R, C);

	//FILLEDI = mxCreateNumericMatrix(R*C * 3, 1, mxSINGLE_CLASS, mxREAL);
	const int dim3[3] = { R, C, 3 };
	//int* ddd = dim3;
	FILLEDI = mxCreateNumericArray(3, dim3, mxSINGLE_CLASS, mxREAL);
	float*	filledI = (float*)mxGetPr(FILLEDI);

	const int dim2[3] = { R,C,2 };
	RETF = mxCreateNumericArray(3, dim2, mxINT32_CLASS, mxREAL);
	int*	retf = (int*)mxGetPr(RETF);

	// global variables

	f = initf;
	im = I;
	mask = M;

	//for (int i = 0; i < R; ++i) {
	//	for (int j = 0; j < C; ++j) {
	//		printf("%d ", initf[i * 2 + j*R * 2]);
	//	}
	//	printf("\n");
	//}
	//printf("\n");
	//for (int i = 0; i < R; ++i) {
	//	for (int j = 0; j < C; ++j) {
	//		printf("%d ", mask[i + j*R]);
	//	}
	//	printf("\n");
	//}
	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			//printf("%d ", initf[i * 2 + j*R * 2]);
			f[i * 2 + j*R * 2] -= 1;
			f[i * 2 + j*R * 2 + 1] -= 1;
		}
		//printf("\n");
	}

	fillOneLevel(initf, I, M, D, level, useline, iternum);

	double sumspc = 0, sumapc = 0;
	for (int j = 0; j < C; ++j) {
		for (int i = 0; i < R; ++i) {
			sumspc += calcSpatialCost(i, j, getFij(i, j));
			sumapc += calcAppearanceCost(i, j, getFij(i, j));
		}
	}
	printf("total spc = %lf, total apc = %lf\n", sumspc, sumapc);
	//filledI = im; 软拷贝不行
	//retf = f;	同上
	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			//printf("%d ", initf[i * 2 + j*R * 2]);
			f[i * 2 + j*R * 2] += 1;
			f[i * 2 + j*R * 2 + 1] += 1;
		}
		//printf("\n");
	}

	memcpy(filledI, im, sizeof(float)*R*C * 3);
	memcpy(retf, f, sizeof(int)*R*C * 2);


	return;
}