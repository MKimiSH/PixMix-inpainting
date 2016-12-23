// NOT WORKING, USED VECTORIZATION INSTEAD

#include"inpaint.h"

inline int randi(int);
inline fij  getFij(int i, int j);
inline void setFij(int i, int j, fij ff);
inline void setFij(int, int, int, int);
inline bool isMask(int i, int j);
inline bool isValid(int i, int j);
inline bool isValid(fij ff);

int* orif;
int* f;
bool* mask;
int R, C, numelf, R1, C1;

inline int randi(int d){
    return (int)((rand() / (float)RAND_MAX) * d);
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
inline void setFij(int i, int j, int ii, int jj){
 	f[HH(i, j, 0, 2)] = ii;
	f[HH(i, j, 1, 2)] = jj;   
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

void initMap(){
    printf("initMap!!\n");
    int i, j, cnt=0;
    if(numelf == 0){
        for(j=0; j<C; ++j){
            for(i=0; i<R; ++i){
				if (isMask(i, j)) {
					do {
						/*fij curf(randi(R), randi(C));*/
						setFij(i, j, randi(R), randi(C));
					} while (!isValid(getFij(i, j)));
				}
                if(!isMask(i,j)){
                    setFij(i,j, i,j);
                }
            }
        }
        return;
    }
	printf("orif exists! %d = %d x %d \n", numelf, R1, C1);
    // orif ~= []
	//for (int i = 0; i < R1; ++i) {
	//	for (int j = 0; j < C1; ++j) {
	//		//printf("%d ", initf[i * 2 + j*R * 2]);
	//		//orif[i * 2 + j*R1 * 2] -= 1;
	//		//orif[i * 2 + j*R1 * 2 + 1] -= 1;
	//		orif[H(i, j, 0, R1, C1, 2)]--;
	//		orif[H(i, j, 1, R1, C1, 2)]--;
	//	}
	//	//printf("\n");
	//}
	printf("orif minus one!\n");
    for (j=0; j<C1; ++j){
        for(i=0; i<R1; ++i){
            fij tlfij(orif[H(i, j, 0, R1, C1, 2)], orif[H(i, j, 1, R1, C1, 2)]);
            setFij(i*2, j*2, tlfij);
			if (isValid(i * 2, j * 2 + 1)) {
				setFij(i * 2, j * 2 + 1, tlfij + fij(0, 1));
			}
            if(isValid(i*2+1, j*2)){
                setFij(i*2+1, j*2, tlfij+fij(1,0));
                if(isValid(i*2, j*2+1)) 
                    setFij(i*2+1, j*2+1, tlfij+fij(1,1));
            }
        }
    }
    for(j=0; j<C; ++j){
        for(i=0; i<R; ++i){
            while(isMask(i,j) && !isValid(getFij(i,j))){
                fij curf(randi(R), randi(C));
                setFij(i,j, curf);
                cnt = cnt+1;
            }
            if(!isMask(i,j)){
                setFij(i,j, i,j);
            }
        }
    }
    printf("%d violations, end init\n", cnt);
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
#define RETF    plhs[0]

#define ORIF    prhs[0]
#define	M       prhs[1]
	//freopen("aa.txt", "w", stdout);
	printf("aaaaa!!!\n");
    // input
    orif = (int*)mxGetData(ORIF);
    mask = (bool*)mxGetData(M);
  
    numelf = mxGetNumberOfElements(ORIF);
    if(numelf){
        R1 = mxGetM(ORIF);
        C1 = mxGetN(ORIF);
    }
	R = mxGetM(M);
	C = mxGetN(M);
	printf("R=%d, C=%d\n", R, C);

	const int dim2[3] = { R,C,2 };
	RETF = mxCreateNumericArray(3, dim2, mxINT32_CLASS, mxREAL);
	int* retf = (int*)mxGetPr(RETF);
	f = new int[R * C * 2];

	initMap();

	//filledI = im; 软拷贝不行
	//retf = f;	同上
	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			//printf("%d ", initf[i * 2 + j*R * 2]);
			//f[i * 2 + j*R * 2] += 1;
			//f[i * 2 + j*R * 2 + 1] += 1;
			f[HH(i, j, 0, 2)] += 1;
			f[HH(i, j, 1, 2)] += 1;
		}
		//printf("\n");
	}
	memcpy(retf, f, sizeof(int)*R*C * 2);
	delete[] f;

	return;
}