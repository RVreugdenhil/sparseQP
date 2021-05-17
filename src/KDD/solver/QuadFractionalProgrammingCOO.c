#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <blas.h>    // dgemv
#include <lapack.h>  // dstegr

#define  abs1(a)         ((a) < 0.0 ? -(a) : (a))
#define  sign1(a)        ((a)==0) ? 0 : (((a)>0.0)?1:(-1))
#define  max1(a,b)       ((a) > (b) ? (a) : (b))
#define  min1(a,b)       ((a) < (b) ? (a) : (b))

double dot1(const double*a,const double*b,const ptrdiff_t n)
{
    double ret = 0;
    ptrdiff_t i;
    for(i=0;i<n;i++)
        ret +=a[i]*b[i];
    return ret;
}
void AxSym(const double*A, const double*x, double*result, const ptrdiff_t n)
{
    /*
     * input:
     * A: n x n
     * x: n x 1
     * output:
     * result = A*x: n x 1
     * NOTE: A is a symmetric matrix, only the the lower triangular part of A is to be referenced
     *       x and the result can not be the same space
     */
    char *UPLO="L";ptrdiff_t N;double ALPHA; ptrdiff_t LDA; ptrdiff_t INCX; double BETA;ptrdiff_t INCY;
    N = n; ALPHA = 1; LDA = n; INCX=1; BETA=0; INCY=1;
    dsymv(UPLO,&N,&ALPHA,A,&LDA,x,&INCX,&BETA,result,&INCY);
}

double computeObj(double s, double u0, double u1, double u2,double d0,double d1,double d2)
{
    return (u2*s*s + u1*s + u0)  /  (d2*s*s + d1*s + d0);
}
double quadfrac2(double u0, double u1, double u2,double d0,double d1,double d2)
{
    double x1;double x2;
    double c2 = u2*d1 - d2*u1;
    double c1 = 2*u2*d0 - 2*d2*u0;
    double c0 = u1*d0 - d1*u0;
    double dd = sqrt(c1*c1 - 4*c2*c0);
    
    if(c2==0)
    {
        x1 = -c0/c1;
        x2 = x1;
    }
    else
    {
        x1 = (-c1+dd)/(2*c2);
        x2 = (-c1-dd)/(2*c2);
    }
    if(computeObj(x1,u0,u1,u2,d0,d1,d2)<computeObj(x2,u0,u1,u2,d0,d1,d2))
    {
        return x1;
    }
    else
    {
        return x2;
    }
    
}

ptrdiff_t pos_max(const double*a,const ptrdiff_t n)
{
    double val = abs1(a[0]);
    ptrdiff_t pos=0;
    ptrdiff_t i ;
    for (i=1;i<n;i++)
    {
        if(abs1(a[i]) > val)
        {
            val = abs1(a[i]);
            pos = i;
        }
    }
    return pos;
}
void pvec(const double*a,const ptrdiff_t n)
{
    // print the vector a
    ptrdiff_t i;
    printf("[ ");
    for (i=0;i<n;i++)
    {
        printf("%f ",a[i]);
    }
    printf("]\n");
}

void solve(double *x,const double *in_x, const double*A, const double*b, const double c, const double *Q, const double *r,const double s,const ptrdiff_t n,const ptrdiff_t max_iter)
{
    
    double *Ax = (double*)malloc(n*sizeof(double));
    double *Qx = (double*)malloc(n*sizeof(double));
    double *up_g =  (double*)malloc(n*sizeof(double));
    double *down_g =  (double*)malloc(n*sizeof(double));
    double *grad =  (double*)malloc(n*sizeof(double));
    double xAx = 0;double bx = 0; double xQx = 0; double rx = 0;
	ptrdiff_t iter, i,j;
 
    
    memcpy(x,in_x,n*sizeof(double));
    AxSym(A, x, Ax, n);
    xAx = dot1(x,Ax,n);
    bx = dot1(x,b,n);
    AxSym(Q,x,Qx,n);
    xQx = dot1(x,Qx,n);
    rx = dot1(r,x,n);
    
    for ( iter=1;iter<=max_iter;iter++)
    {
        double u0,u1,u2,d0,d1,d2,alpha;
        double up,down,fobj;

        up = 0.5*xAx + bx + c;
        down = 0.5*xQx + rx + s;
        fobj =  up / down;
        
        for (j=0;j<n;j++)
        {
            up_g[j] = Ax[j]+b[j];
            down_g[j] = Qx[j]+r[j];
            grad[j] = up_g[j]/down - fobj * down_g[j]/down ;
        }
        i = pos_max(grad,n);
        u0 = 0.5*xAx + bx + c;
        u1 = Ax[i] + b[i];
        u2 = 0.5*A[i*n+i];
        d0 = 0.5*xQx + rx + s;
        d1 = Qx[i] + r[i];
        d2 = 0.5*Q[i*n+i];
        
        alpha =  quadfrac2(u0,u1,u2,d0,d1,d2);
        x[i] = x[i] + alpha;
        
        xAx = xAx + 2*alpha*Ax[i] + alpha*alpha*A[i*n+i];
        xQx = xQx + 2*alpha*Qx[i] + alpha*alpha*Q[i*n+i];
        for (j=0;j<n;j++)
        {
            Ax[j] = Ax[j] + alpha*A[i*n+j];
            Qx[j] = Qx[j] + alpha*Q[i*n+j];
        }
        bx = bx + alpha*b[i];
        rx = rx + alpha*r[i];
        
    }
    
    for (ptrdiff_t i=0;i<n;i++)if(x[i]==0)x[i]=1e-100;
    free(Ax);
    free(Qx);
    free(up_g);
    free(down_g);
    free(grad);
    
}


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* x =            mxGetPr(prhs[0]);
    double* A =            mxGetPr(prhs[1]);
    double* b =            mxGetPr(prhs[2]);
    double c =            mxGetScalar(prhs[3]);
    double* D =            mxGetPr(prhs[4]);
    double* e =            mxGetPr(prhs[5]);
    double f =            mxGetScalar(prhs[6]);
    ptrdiff_t n =         mxGetScalar(prhs[7]);
    ptrdiff_t max_iter =         mxGetScalar(prhs[8]);
    double *out_x;
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    out_x = mxGetPr(plhs[0]);
    solve(out_x,x,A,b,c,D,e,f,n,max_iter);
}

