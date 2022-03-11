#include "_old/utils.h"
#include <math.h>

int dlange_mf(char Norm, int *M, int *N, double *a, int *LDA, double*work)
{
    //mexPrintf("\n -> entering dlange mf \n");
    int m = *M;
    int n = *N;
    int j= 0;
    int i= 0;
    int sum  = 0;
    int res = 0;
    for( j = 0;  j <  n; j++)
    {
        res = 0; 
        sum  = 0;
        for(i = 0; i < m; i++)
        {
            sum = sum + abs(a[i+j]);
        }
        if (sum > res)
        {
            res = sum;
        }
    }
    //mexPrintf("\n -> exit dlange mf \n");
    return res;
}

double dnrm2_mf(int *n, double *x, int *incx)
{
  //mexPrintf("\n -> entering dnrm2_mf mf \n");
  long int ix, nn, iincx;
  double norm, scale, absxi, ssq, temp;

/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DNRM2 := sqrt( x'*x )   

    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to SLASSQ.   
       Sven Hammarling, Nag Ltd.   */

  /* Dereference inputs */
  nn = *n;
  iincx = *incx;

  if( nn > 0 && iincx > 0 )
  {
    if (nn == 1)
    {
      norm = fabs(x[0]);
    }  
    else
    {
      scale = 0.0;
      ssq = 1.0;

      /* The following loop is equivalent to this call to the LAPACK 
         auxiliary routine:   CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

      for (ix=(nn-1)*iincx; ix>=0; ix-=iincx)
      {
        if (x[ix] != 0.0)
        {
          absxi = fabs(x[ix]);
          if (scale < absxi)
          {
            temp = scale / absxi;
            ssq = ssq * (temp * temp) + 1.0;
            scale = absxi;
          }
          else
          {
            temp = absxi / scale;
            ssq += temp * temp;
          }
        }
      }
      norm = scale * sqrt(ssq);
    }
  }
  else
    norm = 0.0;

  //mexPrintf("\n -> exit dnrm2_mf mf \n");
  return norm;

} /* dnrm2_ */


void daxpy_mf(int   *   n_arg,
                      double *   da_arg,
                      double *   dx,
                      int *      incx_arg,
                      double *   dy,
                      int *      incy_arg)
{
    //mexPrintf("\n -> entering daxpy mf \n");
  int i,ix,iy;
  int n=*n_arg;
  double da=*da_arg;
  int incx = *incx_arg;
  int incy = *incy_arg;

  if (n<=0)
    return;

  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*incx;
    if(incy<0)
      iy = (1-n)*incy;
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) 
      dy[iy] += da*dx[ix];

    return;

  } else {

    /* unroll */
    
    for(i=0;i<(n-4);i+=4) {
      dy[i]   += da*dx[i];
      dy[i+1] += da*dx[i+1];
      dy[i+2] += da*dx[i+2];
      dy[i+3] += da*dx[i+3];
    }
    /* continue with current value of i */
    for(;i<n;i++)
      dy[i]   += da*dx[i];
  }
}

#define 	GMX_DOUBLE_EPS   2.2204460492503131e-16 //Double precision accuracy.
 
#define 	GMX_DOUBLE_MAX   1.7976931348623157e+308 //Maximum double precision value - reduced 1 unit in last digit for MSVC.
 
#define 	GMX_DOUBLE_MIN   2.2250738585072014e-308 //Minimum double precision value.
            
void
dgemv_mf(char *trans, 
       int *m__,
       int *n__,
       double *alpha__,
       double *a,
       int *lda__,
       double *x,
       int *incx__,
       double *beta__,
       double *y,
       int *incy__)
{
    //mexPrintf("\n -> entering dgemv mf \n");
  char ch=(*trans);
  int lenx,leny,kx,ky;
  int i,j,jx,jy,ix,iy;
  double temp;

  int m = *m__;
  int n = *n__;
  double alpha = *alpha__;
  double beta = *beta__;
  int incx = *incx__;
  int incy = *incy__;
  int lda = *lda__;

  
  if(n<=0 || m<=0 || (fabs(alpha)<GMX_DOUBLE_MIN && fabs(beta-1.0)<GMX_DOUBLE_EPS))
    return;
   

  if(ch=='N' || ch == 'n') {
    lenx = n;
    leny = m;
  } else {
    lenx = m;
    leny = n;
  }
  
   if(incx>0)
    kx = 1;
  else
    kx = 1 - (lenx -1)*(incx);

  if(incy>0)
    ky = 1;
  else
    ky = 1 - (leny -1)*(incy);
 
  if(fabs(beta-1.0)>GMX_DOUBLE_EPS) {
    if(incy==1) {
      if(fabs(beta)<GMX_DOUBLE_MIN)
	for(i=0;i<leny;i++)
	  y[i] = 0.0;
      else
	for(i=0;i<leny;i++)
	  y[i] *= beta;
    } else {
      /* non-unit incr. */
      iy = ky;
      if(fabs(beta)<GMX_DOUBLE_MIN) 
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] = 0.0;
      else
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] *= beta;
    }
  }
  
  if(fabs(alpha)<GMX_DOUBLE_MIN)
    return;
  
  if(ch=='N' || ch == 'n') {
    jx = kx;
    if(incy==1) {
      for(j=1;j<=n;j++,jx+=incx) 
	if(fabs(x[jx-1])>GMX_DOUBLE_MIN) {
	  temp = alpha * x[jx-1];
	  for(i=1;i<=m;i++)
	    y[i-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jx+=incx) 
	if(fabs(x[jx-1])>GMX_DOUBLE_MIN) {
	  temp = alpha * x[jx-1];
	  iy = ky;
	  for(i=1;i<=m;i++,iy+=incy)
	    y[iy-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    }
  } else {
    /* transpose */
    jy = ky;
    if(incx==1) {
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	for(i=1;i<=m;i++)
	  temp += a[(j-1)*(lda)+(i-1)] * x[i-1];
	y[jy-1] += alpha * temp;
      }
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	ix = kx;
	for(i=1;i<=m;i++,ix+=incx)
	  temp += a[(j-1)*(lda)+(i-1)] * x[ix-1];
	y[jy-1] += alpha * temp;
      }
    }
  }
} 


double ddot_mf(int *n_arg,
                    double *dx,
                    int *incx_arg,
                    double *dy,
                    int *incy_arg)
{
    //mexPrintf("\n -> entering ddot mf \n");
    int i,ix,iy,m;
    int n=*n_arg;
    int incx = *incx_arg;
    int incy = *incy_arg;
    double t1;
    
    if(n<=0)
        return 0.0;
    
    t1 = 0.0;
    
    if(incx!=1 || incy!=1) {
        ix = 0;
        iy = 0;
        if(incx<0)
            ix = (1-n)*incx;
        if(incy<0)
            iy = (1-n)*incy;
        
        for(i=0;i<n;i++,ix+=incx,iy+=incy) 
            t1 += dx[ix] * dy[iy];
        
        return t1;
        
    } else {
        
        m = n%5;
        
        for(i=0;i<m;i++)
            t1 += dx[i] * dy[i];
        
        /* unroll */
        for(i=m;i<n;i+=5) 
            t1  =  t1 + dx[i] * dy[i]   
                +    dx[i+1] * dy[i+1] 
                +    dx[i+2] * dy[i+2] 
                +    dx[i+3] * dy[i+3]   
                +    dx[i+4] * dy[i+4];   
        
        return t1;
    }
}




void
dgemm_mf(char *transa,
                      char *transb,
                      int *m__,
                      int *n__,
                      int *k__,
                      double *alpha__,
                      double *a,
                      int *lda__,
                      double *b,
                      int *ldb__,
                      double *beta__,
                      double *c,
                      int *ldc__)
{
    //mexPrintf("\n -> entering dgemm mf \n");
  char tra=(*transa);
  char trb=(*transb);
  double temp;
  int i,j,l;
  int nrowa,ncola,nrowb;

  int m = *m__;
  int n = *n__;
  int k = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  double alpha = *alpha__;
  double beta  = *beta__;
  
  if(tra=='N' || tra == 'n') {
    nrowa = m;
    ncola = k;
  } else {
    nrowa = k;
    ncola = m;
  }

  if(trb=='N'|| trb == 'n') 
    nrowb = k;
   else 
    nrowb = n;
  
  if(m==0 || n==0 || (( fabs(alpha)<GMX_DOUBLE_MIN || k==0) && fabs(beta-1.0)<GMX_DOUBLE_EPS))
    return;

  if(fabs(alpha)<GMX_DOUBLE_MIN) {
    if(fabs(beta)<GMX_DOUBLE_MIN) {
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] = 0.0;
    } else {
      /* nonzero beta */
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] *= beta;
    }
    return;
  }

  if(trb=='N' || trb == 'n') {
    if(tra=='N' || tra=='n') {
      
      for(j=0;j<n;j++) {
	if(fabs(beta)<GMX_DOUBLE_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(fabs(beta-1.0)>GMX_DOUBLE_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( fabs(b[ j*(ldb) + l ])>GMX_DOUBLE_MIN) {
	    temp = alpha * b[ j*(ldb) + l ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
    } else {
      /* transpose A, but not B */
      for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[j*(ldb)+l];
	  if(fabs(beta)<GMX_DOUBLE_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
      }
    }
  } else {
    /* transpose B */
    if(tra=='N' || tra=='n') {

      /* transpose B, but not A */

      for(j=0;j<n;j++) {
	if(fabs(beta)<GMX_DOUBLE_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(fabs(beta-1.0)>GMX_DOUBLE_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( fabs(b[ l*(ldb) + j ])>GMX_DOUBLE_MIN) {
	    temp = alpha * b[ l*(ldb) + j ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
 
    } else {
      /* Transpose both A and B */
       for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[l*(ldb)+j];
	  if(fabs(beta)<GMX_DOUBLE_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
       }
    }
  }
}