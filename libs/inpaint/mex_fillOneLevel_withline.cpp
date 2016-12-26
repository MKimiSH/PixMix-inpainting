#include"inpaint.h"


void getLines(const mxArray*);
void getParams(const mxArray*);
void getBegEnd();

void fillIWithF();
float calcTotalCost(int, int, fij);
float calcLineCost(int, int, fij);
float calcSpatialCost(int i, int j, fij ff);
float calcAppearanceCost(int i, int j, fij ff);
float calcSumTotalCost();
inline fij  getFij(int i, int j);
inline void setFij(int i, int j, fij ff);
inline bool isMask(int i, int j);
inline bool isNUMask(int i, int j);
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
} line;

line *lines;
int nlines;
// Function Declarations
void fillOneLevel(int* initf, imdata* I, const bool* M, float* D, int level, line*, int);

int		*f;
imdata	*im;
bool	*mask;
bool	*numask;
mxArray *linesptr;
int		R = 0, C = 0;
int		Rbeg, Cbeg, Rend, Cend, nnzMask;
int		RSRounds = 3;
float alphaSp = .005, alphaAp = 0.5, alphaStr = 1 - alphaSp - alphaAp;
float cs_imp = 1, cs_rad = 1;
float kappa = MAX(R, C); // I have no idea how to set this thing...
const float TOL = 1e-2; // for convergence

void getLines(const mxArray* plines) {
	nlines = mxGetNumberOfElements(plines);
	lines = new line[nlines];
	mxArray *tmp;
	double* ft;
	for (int i = 0; i < nlines; ++i) {
		tmp = mxGetFieldByNumber(plines, i, 0); //point1
		ft = (double*)mxGetData(tmp);
		//printf("ft=%x\n", (int)ft);
		//printf("*ft=%lf\n", *ft);
		lines[i].point1[0] = ft[0]-1;
		lines[i].point1[1] = ft[1]-1;

		tmp = mxGetFieldByNumber(plines, i, 1); //point2
		ft = (double*)mxGetData(tmp);
		lines[i].point2[0] = ft[0]-1;
		lines[i].point2[1] = ft[1]-1; // move origin to (0,0)

		tmp = mxGetFieldByNumber(plines, i, 2); //theta
		lines[i].theta = (float)mxGetScalar(tmp);

		tmp = mxGetFieldByNumber(plines, i, 3); //rho
		lines[i].rho = (float)mxGetScalar(tmp);
		printf("line%d: (%f, %f)--(%f, %f), theta = %f, rho = %f\n", i, lines[i].point1[0], lines[i].point1[1],
			lines[i].point2[0], lines[i].point2[1], lines[i].theta, lines[i].rho);
	}
}

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
	alphaStr = 1 - alphaSp - alphaAp;
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
		if (kappa == 0){}
			kappa = MAX(R, C);
	}
	else {
		kappa = MAX(R, C);
	}
	printf("params: %lf, %lf, %lf, %lf, %lf\n", alphaSp, alphaAp, cs_imp, cs_rad, kappa);

}

// 当然可以再快一些，先这样吧
void getBegEnd() {
	//Rbeg = Cbeg = 0; Rend = R - 1; Cend = C - 1; return;
	Rbeg = Cbeg = MAX(R, C);
	Rend = Cend = nnzMask = 0;
	for (int j = 0; j < C; ++j) {
		for (int i = 0; i < R; ++i) {
			if (isMask(i, j)) {
				nnzMask++;
				if (Rbeg > i)Rbeg = i;
				if (Cbeg > j)Cbeg = j;
				if (Rend < i)Rend = i;
				if (Cend < j)Cend = j;
			}
		}
	}
	printf("(%d, %d)->(%d, %d), nnzMask = %d\n", Rbeg, Cbeg, Rend, Cend, nnzMask);
}

// relative convergence criterion
inline bool checkConvergence(float a, float b) {
	
	double x = MAX(a, b), y = MIN(a,b);
	return (x - y) / x < TOL;
}

// main procedure
void fillOneLevel(int* initf, imdata* I, const bool* M, float* D, int level, line* Lines, int iternum) {
	printf("Entering fillOneLevel\n");
	float wrs = .5;
	int rrs = MAX(R, C);
	int i, j, k, r = R, c = C;

	// fill the image with initial guess of f
	//fillIwithF();

	// propagation and random search
	int numIter = iternum;
	double prevCost = 1e20, currCost = 0;
	for (int it = 0; it<numIter; ++it) {
		fillIWithF();
		currCost = calcSumTotalCost();
		if (checkConvergence(prevCost, currCost)) {
			printf("converged before iteration %d\n", it);
			return;
		}
		float prevAvgCost = currCost / nnzMask;
		int calcedPx = 0;
		prevAvgCost *= 1.2;
		prevCost = currCost;
		currCost = 0;
		//printf("begin propagation and search.\n");
		if (it % 2 == 0) {
			fij fijOld, fijLeft, fijUp, fijRS;
			float costOld, costLeft, costUp;
			//float spcOld, spcLeft, spcUp;
			//float apcOld, apcLeft, apcUp;
			float drs, spcRS, apcRS, costRS;
			for (i = Rbeg; i <= Rend; ++i) {
				for (j = Cbeg; j <= Cend; ++j) {
					if (isMask(i, j)) {
						// propagate forward
						fijOld = getFij(i, j);
						costOld = calcTotalCost(i, j, fijOld);

						costLeft = 1e10;
						costUp = 1e10;
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
						
						if (costOld < prevAvgCost) { // loss小就不搜索！！！
							currCost += costOld;
							continue;
						}
						calcedPx++;
						//random search
						drs = wrs*rrs;
						srand((unsigned)time(NULL));
						fijOld = getFij(i, j);
						bool found = false;
						while (drs>D[GG(i, j)] + 2 && !found) {
							int thisRound = RSRounds;
							while (!found && thisRound--) {
								fijRS = randFij(fijOld[0], fijOld[1], drs);
								// should search with center f(i,j) rather than(i,j)
								costRS = 1e10;
								int numTol = 3;
								while (!isValid(fijRS) && numTol--) {
									fijRS = randFij(fijOld[0], fijOld[1], drs);
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
							}
							drs *= wrs;
						}
						currCost += costOld;
					}
				}
			}
			printf("%d pixels updated\n", calcedPx);
		}
		else {
			fij fijOld, fijRight, fijDown, fijRS;
			float costOld, costRight, costDown;
			//float spcOld, spcRight, spcDown;
			//float apcOld, apcRight, apcDown;
			float drs, spcRS, apcRS, costRS;
			for (j = Cend; j >= Cbeg; --j) {
				for (i = Rend; i >= Rbeg; --i) {
					if (isMask(i, j)) {
						// propagate forward
						fijOld = getFij(i, j);
						costOld = calcTotalCost(i, j, fijOld);

						costRight = 1e10;
						costDown = 1e10;
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
						if (costOld < prevAvgCost) {
							currCost += costOld;
							continue;
						}
						calcedPx++;
						//random search
						drs = wrs*rrs;
						srand((unsigned)time(NULL));
						fijOld = getFij(i, j);
						bool found = false;
						while (drs>D[GG(i, j)] + 2 && !found) {
							int thisRound = RSRounds;
							while (!found && thisRound--) {
								fijRS = randFij(fijOld[0], fijOld[1], drs);
								// should search with center f(i,j) rather than(i,j)
								costRS = 1e10;
								int numTol = 3;
								while (!isValid(fijRS) && numTol--) {
									fijRS = randFij(fijOld[0], fijOld[1], drs);
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
							}
							drs *= wrs;
						}
						currCost += costOld;
					}
				}
			}
			printf("%d pixels updated\n", calcedPx);
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
	if (numask != NULL && isNUMask(i, j)) {
		return 1e5;
	}
	float spc = calcSpatialCost(i, j, ff);
	float apc = calcAppearanceCost(i, j, ff);
	float strc = calcLineCost(i, j, ff);
	return spc*alphaSp + apc*alphaAp + strc*alphaStr;
}

// calculate cost_spatial using 8-neighbor
// seemingly, the neighbor pixel which is out of the mask should not be considered
float calcSpatialCost(int i, int j, fij ff) {

	int spc = 0;
	float w = 0.125;
	int maxdist = MAX(R, C) / 2; // this should be like this way...
	float tao = maxdist*maxdist;
	int r_begin = MAX(i - 1, 0), r_end = MIN(i + 1, R - 1);
	int c_begin = MAX(j - 1, 0), c_end = MIN(j + 1, C - 1);
	int r, c, cnt = 0;
	for (c = c_begin; c <= c_end; ++c) {
		for (r = r_begin; r <= r_end; ++r) {
			//% f(p) = fij, f(p + v) = f(r, c, :), v = [r - i], c - j], f(p) + v = fij + v
			fij v = fij(r, c) - fij(i, j), fv = getFij(r, c);
			fij diff = fv - (ff + v);
			spc = spc + MIN(diff.norm2(), tao);
			cnt++;
		}
	}
	return (double)spc * w / cnt;
}


//calculate cost_appear using 5x5 patch
float calcAppearanceCost(int i, int j, fij ff) {

	int apc = 0;
	float w = 1.0 / 25;
	int r_begin = MAX(i - 2, 0), r_end = MIN(i + 2, R - 1);
	int c_begin = MAX(j - 2, 0), c_end = MIN(j + 2, C - 1);
	int r, c, cnt = 0;
	//imdata i1[3] = {}, i2[3] = {};
	for (c = c_begin; c <= c_end; ++c) {
		for (r = r_begin; r <= r_end; ++r) {
			fij v(r - i, c - j), fv = ff + v;
			if (!isValid(fv)) {
				apc = apc + 190000;
				continue;
			}
			for (int k = 0; k < 3; ++k) {
				//i2[k] = im[HH(fv[0], fv[1], k, 3)];
				apc += sq((int)(im[HH(r, c, k, 3)] - im[HH(fv[0], fv[1], k, 3)]));
			}

			//apc = apc + (i1[0] - i2[0])*(i1[0] - i2[0]) +
				//(i1[1] - i2[1])*(i1[1] - i2[1]) + (i1[2] - i2[2])*(i1[2] - i2[2]);
			cnt++;
		}
	}
	return (double)apc*w / cnt / 255; // 8-bit grayscale, [0,1] is too small...compared to spc..
}

float d_cs(fij ff, line ln) {
	int x1 = ln.point1[0], y1 = ln.point1[1];
	int x2 = ln.point2[0], y2 = ln.point2[1];
	//printf("(%d, %d)---(%d, %d)\n", x1, y1, x2, y2);
	fij f1(y1, x1), f2(y2, x2);
	int lenline = (f1 - f2).norm2();
	int ffproj = (f1 - f2) * (ff - f2);
	//printf("lenline = %d, ffproj = %d\n", lenline, ffproj);
	if (ffproj > lenline || ffproj < 0)
		return kappa;

	float x = ff[1], y = ff[0];
	float A = cosd(ln.theta), B = sind(ln.theta), C = -ln.rho;
	float dist = A*x + B*y + C;
	return dist;
}

float calcLineCost(int r, int c, fij ff) {
	//float costs[nlines] = {};
	float omega, dp, dfp, curr;
	float ret = 0, maxomega = 0;
	for (int i = 0; i < nlines; ++i) {
		dp = d_cs(fij(r,c), lines[i]);
		dfp = d_cs(ff, lines[i]);
		omega = cs_imp * exp(-8 * sq(dp / cs_rad));
		curr = sq(dp - dfp) * omega;
		if (maxomega < omega)
		{
			ret = curr;
			maxomega = omega;
		}
		//costs[i] = sq(dp - dfp) * omega;
	}
	return ret;
}

float calcSumTotalCost() {
	float sumspc = 0, sumapc = 0, sumstr = 0;
	for (int j = 0; j < C; ++j) {
		for (int i = 0; i < R; ++i) {
			if (isMask(i, j)) {
				sumspc += calcSpatialCost(i, j, getFij(i, j));
				sumapc += calcAppearanceCost(i, j, getFij(i, j));
				sumstr += calcLineCost(i, j, getFij(i, j));
			}
		}
	}
	//printf("total spc = %lf, total apc = %lf, total strc = %lf\n", sumspc, sumapc, sumstr);
	float ret = sumspc*alphaSp + sumapc*alphaAp + sumstr*alphaStr;
	printf("total cost = %lf, total spc = %lf, total apc = %lf, total strc = %lf\n", ret, sumspc, sumapc, sumstr);
	printf("avg px cost = %lf\n", ret / nnzMask);
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
inline bool isNUMask(int i, int j) {
	return numask[GG(i, j)];
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
#define INNUM		prhs[3]
#define	IND			prhs[4]
#define LEVEL		prhs[5]
#define LINES		prhs[6]
#define ITERNUM		prhs[7]
#define INPARAMS	prhs[8]

	// input
	printf("???\n");
	int*	initf = (int*)mxGetData(INITF);
	imdata*	I = (imdata*)mxGetData(INI);
	printf("??????\n");
	bool*	M = (bool*)mxGetData(INM);
	//numask = (bool*)mxGetData(INNUM);
	int t = mxGetNumberOfElements(INNUM);
	if (t == 0) 
		numask = NULL;
	else numask = (bool*)mxGetData(INNUM);

	float*	D = (float*)mxGetData(IND);

	R = mxGetM(INM);
	C = mxGetN(INM);
	printf("R=%d, C=%d\n", R, C);

	int		level = (int)mxGetScalar(LEVEL);
	getLines(LINES);
	getParams(INPARAMS);
	int		iternum = (int)mxGetScalar(ITERNUM);


	//FILLEDI = mxCreateNumericMatrix(R*C * 3, 1, mxSINGLE_CLASS, mxREAL);
	//const int dim3[3] = { R, C, 3 };
	const int dim3[3] = { 3, R, C };
	//int* ddd = dim3;
	FILLEDI = mxCreateNumericArray(3, dim3, mxUINT8_CLASS, mxREAL);
	imdata*	filledI = (imdata*)mxGetPr(FILLEDI);

	//const int dim2[3] = { R,C,2 };
	const int dim2[3] = { 2,R,C };
	RETF = mxCreateNumericArray(3, dim2, mxINT32_CLASS, mxREAL);
	int*	retf = (int*)mxGetPr(RETF);

	// global variables

	f = initf;
	im = I;
	mask = M;
	getBegEnd();

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

	fillOneLevel(initf, I, M, D, level, lines, iternum);

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

	memcpy(filledI, im, sizeof(imdata)*R*C * 3);
	memcpy(retf, f, sizeof(int)*R*C * 2);
	if (lines)
		delete [] lines;

	return;
}