#include"inpaint.h"


//void getLines(const mxArray*);
void getParams(const mxArray*);

void fillIWithF();
float calcTotalCost(int, int, fij);
//float calcLineCost(int, int, fij);
float calcSpatialCost(int i, int j, fij ff);
float calcAppearanceCost(int i, int j, fij ff);
float calcAdditionalAppCost(int, int, fij);
float calcSumTotalCost();
inline fij  getFij(int i, int j);
inline void setFij(int i, int j, fij ff);
inline bool isMask(int i, int j);
//inline bool isNUMask(int i, int j);
inline bool isValid(int i, int j);
inline bool isValid(fij ff);
inline bool checkConvergence(float, float);
inline fij	randFij(int i, int j, float dist);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


// global variables
typedef struct linestr {
	float* point1;
	float* point2;
	float theta;
	float rho;
	linestr() {
		point1 = new float[2];
		point2 = new float[2];
		theta = 0;
		rho = 0;
	}
	linestr(float* p1, float* p2, float r, float t) :point1(p1), point2(p2), rho(r), theta(t) {}

	/*~linestr() { delete [] point1; delete [] point2; }*/
} line; // unused

//line *lines;
//int nlines;
// Function Declarations
void fillOneLevel(int*, float*, const bool*, float*, float*, line*, int);

int		*f;
float	*im;
float	*refim;
bool	*mask;
//bool	*numask;
mxArray *linesptr;
int		R = 0, C = 0;
float alphaSp = .005, alphaAp = 0.5, alphaAddt = 1 - alphaSp - alphaAp;//alphaStr = 1 - alphaSp - alphaAp;
float cs_imp = 1, cs_rad = 1; // unused
float kappa = MAX(R, C); // unused
const float TOL = 10; // for convergence


void getParams(const mxArray* pparams) {
	mxArray *tmp;
	tmp = mxGetField(pparams, 0, "alphaSp");
	if (tmp) {
		alphaSp = (float)mxGetScalar(tmp);
	}
	tmp = mxGetField(pparams, 0, "alphaAp");
	if (tmp) {
		alphaAp = (float)mxGetScalar(tmp);
	}
	//alphaStr = 1 - alphaSp - alphaAp;
	alphaAddt = 1 - alphaSp - alphaAp;
	tmp = mxGetField(pparams, 0, "cs_imp");
	if (tmp) {
		cs_imp = (float)mxGetScalar(tmp);
	}
	tmp = mxGetField(pparams, 0, "cs_rad");
	if (tmp) {
		cs_rad = (float)mxGetScalar(tmp);
	}
	tmp = mxGetField(pparams, 0, "kappa");
	if (tmp) {
		kappa = (float)mxGetScalar(tmp);
		if (kappa == 0) {}
		kappa = MAX(R, C);
	}
	else {
		kappa = MAX(R, C);
	}
	printf("params: %lf, %lf, %lf, %lf, %lf\n", alphaSp, alphaAp, cs_imp, cs_rad, kappa);

}


inline bool checkConvergence(float a, float b) {
	//return sq(a - b) < TOL;
	return false;
}

// main procedure
void fillOneLevel(int* initf, float* I, const bool* M, float* D, float* refI, line* Lines, int iternum) {
	printf("Entering fillSecondImage\n");
	float wrs = .5;
	int rrs = MAX(R, C);
	int i, j, k, r = R, c = C;

	// fill the image with initial guess of f
	//fillIwithF();

	// propagation and random search
	int numIter = iternum;
	float prevCost = 1e20, currCost = 0;
	for (int it = 0; it<numIter; ++it) {
		fillIWithF();
		calcSumTotalCost();
		if (checkConvergence(prevCost, currCost)) {
			printf("converged before iteration %d\n", it);
			return;
		}
		prevCost = currCost;
		currCost = 0;
		//printf("begin propagation and search.\n");
		if (it % 2 == 0) {
			fij fijOld, fijLeft, fijUp, fijRS;
			float costOld, costLeft, costUp;
			//float spcOld, spcLeft, spcUp;
			//float apcOld, apcLeft, apcUp;
			float drs, spcRS, apcRS, costRS;
			for (i = 0; i < R; ++i) {
				for (j = 0; j < C; ++j) {
					if (isMask(i, j)) {
						// propagate forward
						fijOld = getFij(i, j);
						costLeft = 1e10;
						costUp = 1e10;
						costOld = calcTotalCost(i, j, fijOld);
						fijLeft = getFij(i, j - 1) + fij(0, 1);
						fijUp = getFij(i - 1, j) + fij(1, 0);
						if (isMask(i, j - 1) && isValid(fijLeft))
						{
							costLeft = calcTotalCost(i, j, fijLeft);
						}
						if (isMask(i - 1, j) && isValid(fijUp))
						{
							costUp = calcTotalCost(i, j, fijUp);
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
						fijOld = getFij(i, j);
						while (drs>D[GG(i, j)] + 3) {
							fijRS = randFij(fijOld[0], fijOld[1], drs);
							// should search with center f(i,j) rather than(i,j)
							costRS = 1e10;
							int numTol = 3;
							while (!isValid(fijRS) && numTol > 0) {
								fijRS = randFij(fijOld[0], fijOld[1], drs);
								numTol--;
							}
							if (numTol <= 0) {
								continue;
							}
							costRS = calcTotalCost(i, j, fijRS);
							if (costRS < costOld) {
								setFij(i, j, fijRS);
								costOld = costRS;
								break;
							}
							drs *= wrs;
						}
						currCost += costOld;
					}
				}
			}

		}
		else {
			fij fijOld, fijRight, fijDown, fijRS;
			float costOld, costRight, costDown;
			//float spcOld, spcRight, spcDown;
			//float apcOld, apcRight, apcDown;
			float drs, spcRS, apcRS, costRS;
			for (j = C - 1; j >= 0; --j) {
				for (i = R - 1; i >= 0; --i) {
					if (isMask(i, j)) {
						// propagate forward
						fijOld = getFij(i, j);
						costRight = 1e10;
						costDown = 1e10;
						costOld = calcTotalCost(i, j, fijOld);
						fijRight = getFij(i, j + 1) - fij(0, 1);
						fijDown = getFij(i + 1, j) - fij(1, 0);
						if (isMask(i, j + 1) && isValid(fijRight))
						{
							costRight = calcTotalCost(i, j, fijRight);
						}
						if (isMask(i + 1, j) && isValid(fijDown))
						{
							costDown = calcTotalCost(i, j, fijDown);
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
						fijOld = getFij(i, j);

						while (drs>D[GG(i, j)] + 3) {
							fijRS = randFij(fijOld[0], fijOld[1], drs);
							costRS = 1e10;
							int numTol = 3;
							while (!isValid(fijRS) && numTol > 0) {
								fijRS = randFij(fijOld[0], fijOld[1], drs);
								numTol--;
							}
							if (numTol <= 0) {
								continue;
							}
							costRS = calcTotalCost(i, j, fijRS);
							if (costRS < costOld) {
								setFij(i, j, fijRS);
								costOld = costRS;
								break;
							}
							drs *= wrs;
						}
						currCost += costOld;
					}
				}
			}
		}
		//printf("end propagation and search.\n");
	}



}

void fillIWithF() {
	//printf("entering fillIwithF\n");
	int i, j, k, r = R, c = C;

	for (i = 0; i < R; ++i) {
		for (j = 0; j < C; ++j) {
			//if (mask[G(i, j, r, c)])
			if (isMask(i, j))
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

float calcTotalCost(int i, int j, fij ff) {
	// if (numask != NULL && isNUMask(i, j)) {
		// return 1e5;
	// }
	float spc = calcSpatialCost(i, j, ff);
	float apc = calcAppearanceCost(i, j, ff);
	float addtc = calcAdditionalAppCost(i, j, ff);
	return spc*alphaSp + apc*alphaAp + addtc*alphaAddt;
}

// calculate cost_spatial using 8-neighbor
// seemingly, the neighbor pixel which is out of the mask should not be considered
float calcSpatialCost(int i, int j, fij ff) {
	float spc = 0;
	float w = 0.125;
	int maxdist = MAX(R, C) / 2; // this should be like this way...
	float tao = maxdist*maxdist;
	int r_begin = MAX(i - 1, 0), r_end = MIN(i + 1, R - 1);
	int c_begin = MAX(j - 1, 0), c_end = MIN(j + 1, C - 1);
	int r, c, cnt = 0;
	for (r = r_begin; r <= r_end; ++r) {
		for (c = c_begin; c <= c_end; ++c) {
			//% f(p) = fij, f(p + v) = f(r, c, :), v = [r - i], c - j], f(p) + v = fij + v
			fij v = fij(r, c) - fij(i, j), fv = getFij(r, c);
			fij diff = fv - (ff + v);
			spc = spc + MIN(diff.norm2(), tao);
			cnt++;
		}
	}
	spc = spc*w / cnt;
	return spc;
}

//calculate cost_appear using 5x5 patch
float calcAppearanceCost(int i, int j, fij ff) {

	float apc = 0;
	float w = 1.0 / 25;
	int r_begin = MAX(i - 2, 0), r_end = MIN(i + 2, R - 1);
	int c_begin = MAX(j - 2, 0), c_end = MIN(j + 2, C - 1);
	int r, c, cnt = 0;
	float i1[3] = {}, i2[3] = {};
	for (r = r_begin; r <= r_end; ++r) {
		for (c = c_begin; c <= c_end; ++c) {
			fij v(r - i, c - j), fv = ff + v;
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
			cnt++;
		}
	}
	apc = apc*w;
	return apc * 255 / cnt; // 8-bit grayscale, [0,1] is too small...compared to spc..
}

float calcAdditionalAppCost(int i, int j, fij ff) {
	float addtc = 0;
	float w = 1.0 / 25;
	int r_begin = MAX(i - 2, 0), r_end = MIN(i + 2, R - 1);
	int c_begin = MAX(j - 2, 0), c_end = MIN(j + 2, C - 1);
	int r, c, cnt = 0;
	float i1[3] = {}, i2[3] = {};
	for (r = r_begin; r <= r_end; ++r) {
		for (c = c_begin; c <= c_end; ++c) {
			fij v(r - i, c - j), fv = ff + v;
			for (int k = 0; k < 3; ++k) {
				i1[k] = im[HH(r, c, k, 3)];
				i2[k] = refim[HH(r, c, k, 3)];
			}
			addtc += sq(i1[0] - i2[0]) + sq(i1[1] - i2[1]) + sq(i1[2] - i2[2]);
			cnt++;
		}
	}
	addtc = addtc * w;
	return addtc * 255 / cnt; // 8-bit grayscale, [0,1] is too small...compared to spc..
}

float calcSumTotalCost() {
	float sumspc = 0, sumapc = 0, sumstr = 0, sumaddtc = 0;
	for (int j = 0; j < C; ++j) {
		for (int i = 0; i < R; ++i) {
			if (isMask(i, j)) {
				sumspc += calcSpatialCost(i, j, getFij(i, j));
				sumapc += calcAppearanceCost(i, j, getFij(i, j));
				//sumstr += calcLineCost(i, j, getFij(i, j));
				sumaddtc += calcAdditionalAppCost(i, j, getFij(i, j));
			}
		}
	}
	//printf("total spc = %lf, total apc = %lf, total strc = %lf\n", sumspc, sumapc, sumstr);
	float ret = sumspc*alphaSp + sumapc*alphaAp + sumaddtc*alphaAddt;//sumstr*alphaStr;
	printf("total cost = %lf, total spc = %lf, total apc = %lf, total addtc = %lf\n", ret, sumspc, sumapc, sumaddtc);
	return ret;
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
inline fij	randFij(int i, int j, float dist) {
	fij ret;
	ret[0] = (int)(((rand() / (float)RAND_MAX) - 0.5) * 2 * dist) + i;
	ret[1] = (int)(((rand() / (float)RAND_MAX) - 0.5) * 2 * dist) + j;
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
#define INREFI		prhs[4]
#define LINES		prhs[5]
#define ITERNUM		prhs[6]
#define INPARAMS	prhs[7]
//#define INNUM		prhs[3]
//#define LEVEL		prhs[5]

	// input
	printf("???\n");
	int*	initf = (int*)mxGetData(INITF);
	float*	I = (float*)mxGetData(INI);
	bool*	M = (bool*)mxGetData(INM);
	float*	D = (float*)mxGetData(IND);
	float*	rI = (float*)mxGetData(INREFI);

	R = mxGetM(INM);
	C = mxGetN(INM);
	printf("R=%d, C=%d\n", R, C);

	getParams(INPARAMS);
	int		iternum = (int)mxGetScalar(ITERNUM);


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
	refim = rI;
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

	fillOneLevel(initf, I, M, D, rI, (line*)NULL, iternum);

	calcSumTotalCost();

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