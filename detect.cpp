#include <mex.h>
#include <math.h>



bool pointInPolygon(double x,double y,double polyX[], double polyY[],int polySides) {
    int   j = polySides-1;
    bool  oddNodes = false;

    for (int i = 0;i < polySides; i++)
    {
        if((polyY[i] < y && polyY[j] >= y
            ||   polyY[j] < y && polyY[i] >= y)
           && (polyX[i] <= x || polyX[j] <= x))
        {
            oddNodes ^= (polyX[i] + (y - polyY[i]) / (polyY[j] - polyY[i]) * (polyX[j] - polyX[i]) < x);
        }
        j = i;
    }
    return oddNodes;
}

int dis(double a[], double b[], double v[], int l)
{
    int d = 0;
    for (int i = 0; i < 3; i++) {
        if(abs(a[i] - b[i*l]) > v[i])
        {
            d = 1;
            break;
        }
    }
    return d;
}

/**
   prhs[0] : landmarks(,:2)
   prhs[1] : landmarks(,:1)
   prhs[2] : m3
   prhs[3] : m7
   prhs[4] : U
   prhs[5] : v
   prhs[6] : [minx, maxx, miny, maxy]
   prhs[7] : gamma
**/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *landmarks1 = mxGetPr(prhs[0]);
    double *landmarks2 = mxGetPr(prhs[1]);
    double *m3 = mxGetPr(prhs[2]);
    double *m7 = mxGetPr(prhs[3]);
    double *U = mxGetPr(prhs[4]);
    double *v = mxGetPr(prhs[5]);
    double *board = mxGetPr(prhs[6]);
    double gamma = mxGetScalar(prhs[7]);

    //#define OUT plhs[0];
    int l = mxGetM(prhs[0]);
    
    size_t K = mxGetNumberOfDimensions(prhs[2]);
    const mwSize *N = mxGetDimensions(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(N[0], N[1], mxREAL);
    double *out = mxGetPr(plhs[0]);
    
    /* 
    for (int i = 0; i< m; i++) {
        for (int j = 0; j < n; j++) {
            mexPrintf("%lf ",landmarks[i+j*m]);
        }
    }
    */
    
    for (int i = (int)board[0]; i < board[1]; i++) {
        for (int j = (int)board[2]; j < board[3]; j++) {
            if(pointInPolygon(i,j, landmarks1, landmarks2, l))
            {
                double a[3];
                a[0] =  (m3[i+N[0]*(j+N[1]*0)]*2 + m7[i+N[0]*(j+N[1]*0)])/3;
                a[1] =  (m3[i+N[0]*(j+N[1]*1)]*2 + m7[i+N[0]*(j+N[1]*1)])/3;
                a[2] =  (m3[i+N[0]*(j+N[1]*2)]*2 + m7[i+N[0]*(j+N[1]*2)])/3;
                int sumU = 0;
                for (int k = 0; k < l; k++) {
                    sumU += dis(a, &(U[k]), v, l);
                }
                // mexPrintf("%d\n", sumU);

                
                if(sumU > l * gamma)
                {
                    out[i+j*N[0]] = 1;
                }
            }
        }
    }
}