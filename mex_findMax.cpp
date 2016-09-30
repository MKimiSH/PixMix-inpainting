#include <stdio.h>
#include <math.h>
#include "mex.h"
#include <ctime>
#include "matrix.h"
#include <algorithm>
//#include "opencv2/core/core.hpp"
//#include "opencv2/highgui/highgui.hpp"
//#include "pthread.h"
#include <time.h>
#include <vector>
#include <conio.h>

// Function Declarations
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// global variables
// NONE

void findMax(float* A, float* B, float* out, int n, int m) {
#define MAX(a,b)	((a>b)?a:b)
	//int m = mxGetM(A);
	//int n = mxGetN(A);

	int i, j;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			float a = A[i*n + j];
			float b = B[i*n + j];
			out[i*n + j] = MAX(a,b);
		}
	}
}


//=====================================================================
//mexFunction Entry point
//=====================================================================

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#define OUT			plhs[0]

#define INA			prhs[0]
#define	INB			prhs[1]

	float* A = (float*)mxGetData(INA);
	float* B = (float*)mxGetData(INB);

	int m = mxGetM(INA);
	int n = mxGetN(INA);
	int mb = mxGetM(INB);
	int nb = mxGetN(INB);
	if (m != mb || n != nb) {
		mexErrMsgTxt("matrix dimension disagree.");
	}

	OUT = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
	float* out = (float*)mxGetPr(OUT);
	
	printf("m=%d, n=%d\n", m, n);
	findMax(A, B, out, n, m);
	return;
}