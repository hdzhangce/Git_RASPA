/*****************************************************************************************************
    nonuniform_cubic_spline_1d.c -  description
    -------------------------------------------
    Copyright © 2006-2011 David Dubbeldam, Sofia Calero, Donald E. Ellis 
                          and Randall Q. Snurr. All Rights Reserved.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/depa/webdex/quimfis/miembros/Web_Sofia/Sofia.htm
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/
 *****************************************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "cubic_spline_1d.h"
#include "nonuniform_cubic_spline_1d.h"
#include "matrix.h"
#include "linear_equations.h"

int diag5dec(int   n,REAL *ld2, REAL *ld1,REAL *d,
             REAL *ud1,REAL *ud2);
int diag5sol(int n,REAL *ld2,REAL *ld1,REAL *d,
             REAL *ud1,REAL *ud2,REAL *b);
int diag5pddec(int n,REAL *d,REAL *ud1,REAL *ud2);
int diag5pdsol(int n,REAL *d,REAL *ud1,REAL *ud2,REAL *b);
int diag5pd(int mod,int n,REAL *d,REAL *ud1,
            REAL *ud2,REAL *b);
int diag5(int mod,int n,REAL *ld2,REAL *ld1, REAL *d,
          REAL *ud1,REAL *ud2,REAL *b);
int fzyfsy(int n,REAL * md,REAL *ud1, REAL *ud2,
           REAL *rs,REAL *x,REAL *cmd,REAL *cld_1,
           REAL *cld_2,REAL *cld_l2,REAL *cld_l1,
           REAL *bud_1,REAL *bud_2,REAL *brs_2,REAL *brs_1);
int fzyfsz(int n,REAL *md,REAL *ud1,REAL *ud2,REAL *cmd,
           REAL *cld_1,REAL *cld_2,REAL *cld_l2,REAL *cld_l1,
           REAL *bud_1,REAL *bud_2,REAL *brs_2, REAL *brs_1);
int fzyfsl(int n,REAL *rs,REAL *x,REAL *cmd,REAL *cld_1,
           REAL *cld_2,REAL *cld_l2,REAL *cld_l1,REAL *bud_1,
           REAL *bud_2,REAL *brs_2,REAL *brs_1);



CUBIC_SPLINE CreateCubicSpline(int m,REAL *x,REAL *y,int boundary_condition,
                      REAL left_boundary,REAL right_boundary)
{
  int n,i;
  REAL *h,*lower,*diag,*upper,*lowrow=NULL,*ricol=NULL;
  CUBIC_SPLINE spline;

  spline.n=m; 
  n=m-1;
  spline.a=(REAL*)calloc(m+1,sizeof(REAL));
  spline.b=(REAL*)calloc(m+1,sizeof(REAL));
  spline.c=(REAL*)calloc(m+1,sizeof(REAL));
  spline.d=(REAL*)calloc(m+1,sizeof(REAL));

  spline.x=(REAL*)calloc(m+1,sizeof(REAL));
  spline.y=(REAL*)calloc(m+1,sizeof(REAL));

  for(i=0;i<m;i++)
  {
    spline.x[i]=x[i];
    spline.y[i]=y[i];
    spline.a[i]=y[i];
  }

  h=(REAL*)calloc(n+1,sizeof(REAL));
  lower=(REAL*)calloc(n+1,sizeof(REAL));
  diag=(REAL*)calloc(n+1,sizeof(REAL));
  upper=(REAL*)calloc(n+1,sizeof(REAL));

  if(boundary_condition==PERIODIC_SPLINE)
    if(n>2)
    {
      lowrow=(REAL*)calloc(n-2+1,sizeof(REAL));
      ricol=(REAL*)calloc(n-2+1,sizeof(REAL));
    }
  
  for(i=0;i<n;i++) 
    h[i]=x[i+1]-x[i];

  for(i=0;i<n-1;i++)    /* form linear system */
  {
    spline.c[i]=3.0*((y[i+2]-y[i+1])/h[i+1]-(y[i+1]-y[i])/h[i]);
    diag[i]=2.0*(h[i]+h[i+1]);
    lower[i+1]=upper[i]=h[i+1];
  }

  switch (boundary_condition)
  {
    case NOT_A_NODE:
      if(n==2)
      {
        spline.c[0]/=3.0;
        diag[0]*=0.5;
      }
      else
      {
        spline.c[0]*=h[1]/(h[0]+h[1]);
        spline.c[n-2]*=h[n-2]/(h[n-1]+h[n-2]);
        diag[0]-=h[0],
        diag[n-2]-=h[n-1],
        upper[0]-=h[0],
        lower[n-2]-=h[n-1];
      }
      break;
    case FIRST_DERIVATIVE:
      spline.c[0]-=1.5*((y[1]-y[0])/h[0]-left_boundary);
      spline.c[n-2]-=1.5*(right_boundary-(y[n]-y[n-1])/h[n-1]);
      diag[0]-=0.5*h[0];
      diag[n-2]-=0.5*h[n-1];
      break;
    case SECOND_DERIVATIVE:
      spline.c[0]-=h[0]*0.5*left_boundary;
      spline.c[n-2]-=h[n-1]*0.5*right_boundary;
      break;
    case THIRD_DERIVATIVE:
      spline.c[0]+=0.5*left_boundary*h[0]*h[0];
      spline.c[n-2]-=0.5*right_boundary*h[n-1]*h[n-1];
      diag[0]+=h[0];
      diag[n-2]+=h[n-1];
      break;
    case PERIODIC_SPLINE:
      spline.c[n-1]=3.0*((y[1]-y[0])/h[0]-(y[n]-y[n-1])/h[n-1]);
      if(n>2)
      {
        diag[n-1]=2.0*(h[0]+h[n-1]);
        ricol[0]=lowrow[0] = h[0];
      }
  }

  switch(n)
  {
    case 2:
      if(boundary_condition==4)
      {
        spline.c[1]=3.0*(y[0]-y[1])/(x[2]-x[1])/(x[1]-x[0]);
        spline.c[2]=-spline.c[1];
      }
      else
        spline.c[1]=spline.c[0]/diag[0];
      break;
    default:    
      if (boundary_condition==PERIODIC_SPLINE)
        tzdiag(n, lower, diag, upper, lowrow,ricol, spline.c, 0);
      else
        trdiag(n - 1, lower, diag, upper, spline.c, 0);
      for(i=n;i!=0;i--)
        spline.c[i]=spline.c[i-1];
      break;
  }  

  switch (boundary_condition)
  {
    case NOT_A_NODE:
      if(n==2)
        spline.c[0]=spline.c[2]=spline.c[1];
      else
      {
        spline.c[0]=spline.c[1]+h[0]*(spline.c[1]-spline.c[2])/h[1];
        spline.c[n]=spline.c[n-1]+h[n-1]*(spline.c[n-1]-spline.c[n-2])/h[n-2];
      }
      break;
    case FIRST_DERIVATIVE:
      spline.c[0]=1.5*((y[1]-y[0])/h[0]-left_boundary);
      spline.c[0]=(spline.c[0]-spline.c[1]*h[0]*0.5)/h[0];
      spline.c[n]=-1.5*((y[n]-y[n-1])/h[n-1]-right_boundary);
      spline.c[n]=(spline.c[n]-spline.c[n-1]*h[n-1]*0.5)/h[n-1];
      break;
    case SECOND_DERIVATIVE:
      spline.c[0]=left_boundary*0.5;
      spline.c[n]=right_boundary*0.5;
      break;
    case THIRD_DERIVATIVE:
      spline.c[0]=spline.c[1]-left_boundary*0.5*h[0];
      spline.c[n]=spline.c[n-1]+right_boundary*0.5*h[n-1];
      break;
    case PERIODIC_SPLINE:
      spline.c[0]=spline.c[n];
  }
  for(i=0;i<n;i++)
  {
    spline.b[i]=(y[i+1]-y[i])/h[i]-h[i]*(spline.c[i+1]+2.0*spline.c[i])/3.0;
    spline.d[i]=(spline.c[i+1]-spline.c[i])/(3.0*h[i]);
  }
  return spline;
}


static REAL hf[10];   /* aux vector that may not be altered after     */
                      /* the first call of a run                      */

// Compute the coefficients of a cubic fitting spline for given first
// derivatives at the end points. 
int glsp1a(int n,REAL *xn,REAL *fn,REAL *w,
           REAL marg_0,REAL marg_n,int rep,
           REAL *a,REAL *b,REAL *c,REAL *d,
           REAL *h,REAL *h1,REAL *h2,
           REAL *md,REAL *ud1,REAL *ud2)
{
  int i,k,error;
  REAL h_var_1,h_var_2;

  if(rep!=0&&rep!= 1) return (3);
  if(!rep) // First call: determine aux values and LU factorization
  {
    for(i=0;i<=n-1;++i)
    {
      h[i]=xn[i+1]-xn[i];
      h1[i]=1./h[i];
      c[i]=h1[i]*h1[i];
      b[i]=6./w[i];
    }
    b[n]=6./w[n];
    for(i=0;i<=n-2;++i)
      h2[i]=h1[i]+h1[i+1];
    // second  co-diagonal
    for(i=1;i<=n-3;++i)
      ud2[i]=b[i+1]*h1[i]*h1[i+1];

    // first co-diagonal
    for(i=1;i<=n-2;++i)
      ud1[i]=h[i]-b [i]*h1[i]*h2[i-1]-b[i+1]*h1[i]*h2[i];

    // main diagonal
    for(i=1;i<=n-1;++i)
    {
      k=i-1;
      md[i]=2.*(h[k]+h[i])+b[k]*c[k]
             +b[i]*h2[k]*h2[k]+b[i+1]*c[i];
    }

    // The global aux vector hf is used to alter the corners in the system
    // matrix. Its 10 elements must be available for repeated calls.
    hf[0]=h[0]-b[0]*c[0]-b[1]*h2[0]*h1[0];
    hf[1]=b[1]*h1[0]*h1[1];
    hf[2]=b[n-1]*h1[n-2]*h1[n-1];
    hf[3]=h[n-1]-b[n-1]*h2[n-2]*h1[n-1]-b[n]*c[n-1];
    hf[4]=h1[0]*(b[1]+b[0])+2.*h[0]*h[0];
    hf[5]=h1[n-1]*(b[n]+b[n-1])+2.*h[n-1]*h[n-1];
    hf[6]=b[1]*h2[0]+b[0]*h1[0]-h[0]*h[0];
    hf[7]=b[n-1]*h2[n-2]+b[n]*h1[n-1]-h[n-1]*h[n-1];
    hf[8]=b[1]*h1[1];
    hf[9]=b[n-1]*h1[n-2];

    md[1]+=hf[0]/hf[4]*hf[6];
    ud1[1]-=hf[0]/hf[4]*hf[8];
    md[2]-=hf[1]/hf[4]*hf[8];
    md[n-2]-=hf[2]/hf[5]*hf[9];
    ud1[n-2]+=hf[2]/hf[5]*hf[7];
    md[n-1] +=hf[3]/hf[5]*hf[7];
  }

  // Compute right hand side
  h_var_1=(fn[1]-fn[0])*h1[0];
  for(i=1;i<=n-1;++i,h_var_1=h_var_2)
  {
    h_var_2=(fn[i+1]-fn[i])*h1[i];
    c[i]=3.*(h_var_2-h_var_1);
  }

  h_var_1=3.*((fn[1]-fn[0])-marg_0*h[0]);
  h_var_2=3.*(marg_n*h[n-1]-(fn[n]-fn[n-1]));

  c[1]-=hf[0]/hf[4]*h_var_1;
  c[2]-=hf[1]/hf[4]*h_var_1;
  c[n-2]-=hf[2]/hf[5]*h_var_2;
  c[n-1]-=hf[3]/hf[5]*h_var_2;

  // Compute the  c[i], i=1, ..., n-1
  // by solving the linear system via diag5pd
  if(!rep)
  {
    error=diag5pd(0,n-1,md+1,ud1+1,ud2+1,c+1);
    if(error!=0)
    {
      if(error==1)
        return 2;
      else
        return 1;
    }
  }
  else
    diag5pd(2,n-1,md+1,ud1+1,ud2+1,c+1);

  // compute remaining spline coefficients
  c[0]=(h_var_1+c[1]*hf[6]-c[2]*hf[8])/hf[4];
  c[n]=(h_var_2+c[n-1]*hf[7]-c[n-2]*hf[9])/hf[5];

  a[0]=fn[0]+2./w[0]*h1[0]*(c[0]-c[1]);
  for(i=1;i<=n-1;++i)
  {
    k=i-1;
    a[i]=fn[i]-2./w[i]*(c[k]*h1[k]-h2[k]*c[i]+c[i+1]*h1[i]);
  }
  a[n]=fn[n]-2./w[n]*h1[n-1]*(c[n-1]-c[n]);

  for(i=0;i<=n-1;++i)
  {
    k=i+1;
    b[i]=h1[i]*(a[k]-a[i])-h[i]/3.*(c[k]+2.*c[i]);
    d[i]=h1[i]/3.*(c[k]-c[i]);
  }
  return 0;
}

// Compute the coefficients of a cubic fitting spline for given second
// derivatives at the end points.
int glsp2a(int n,REAL *xn,REAL *fn,REAL *w,
           REAL marg_0,REAL marg_n,int rep,
           REAL *a,REAL *b,REAL *c,REAL *d,
           REAL *h,REAL *h1,REAL *h2,
           REAL *md,REAL *ud1,REAL *ud2)
{
  int i,k,error;
  REAL h_var_1,h_var_2;

  if (rep!=0&&rep!=1) return (3);
  if(!rep)  // Compute aux values and three diagonals of linear system
  {
    for(i=0;i<=n-1;++i)
    {
      h[i]=xn[i+1]-xn[i];
      h1[i]=1./h[i];
      c[i]=h1[i]*h1[i];
      b[i]=6./w[i];
    }
    b[n]=6./w[n];

    for(i=0;i<=n-2;++i)
      h2[i]=h1[i]+h1[i+1];

    // second co-diagonal
    for(i=1;i<=n-3;++i)
      ud2[i]=b[i+1]*h1[i]*h1[i+1];

    // first co-diagonal
    for(i=1;i<=n-2;++i)
      ud1[i]=h[i]-b[i]*h1[i]*h2[i-1]-b[i+1]*h1[i]*h2[i];

    // main diagonal
    for(i=1;i<=n-1;++i)
    {
      k=i-1;
      md[i]=2.*(h[k]+h [i])+b[k]*c[k]
            +b[i]*h2[k]*h2[k]+b[i+1]*c[i];
    }
  }

  // Compute right hand side
  c[0]=0.5*marg_0;
  c[n]=0.5*marg_n;

  h_var_2=(fn[2]-fn[1])*h1[1];
  h_var_1=(fn[3]-fn[2])*h1[2];
  c[1]=3.*(h_var_2-(fn[1]-fn[0])*h1[0])-c[0]*
       (h[0]-6./w[0]*h1[0]*h1[0]-6./w[1]*h1[0]*h2[0]);
  c[2]=3.*(h_var_1-h_var_2)-c[0]*(6./w[1])*h1[0]*h1[1];
  for(i=3;i<=n-3;++i,h_var_1=h_var_2)
  {
    h_var_2=(fn[i+1]-fn[i])*h1[i];
    c[i]=3.*(h_var_2-h_var_1);
  }
  h_var_2=(fn[n-1]-fn[n-2])*h1[n-2];
  c[n-2]=3.*(h_var_2-h_var_1)
             -c[n]*6./w[n-1]*h1[n-2]*h1[n-1];
  c[n-1]=3.*((fn[n]-fn[n-1])*h1[n-1]-h_var_2)
         -c[n]*(h[n-1]-6./w[n-1]*h1[n-1]*h2[n-2]
                              -6./w[n]*c[n-1]);

  // compute coefficients  c[i], i=1, ..., n-1
  // via LU factorization in  diag5pd
  if(!rep)
  {
    error=diag5pd(0,n-1,md+1,ud1+1,ud2+1,c+1);
    if (error != 0)
    {
      if (error == 1)
        return 2;
      else
        return 1;
    }
  }
  else
    diag5pd(2,n-1,md+1,ud1+1,ud2+1,c+1);

  // compute remaining spline coefficients
  a[0]=fn[0]+2./w[0]*h1[0]*(c[0]-c[1]);
  for(i=1;i<=n-1;++i)
  {
    k=i-1;
    a[i]=fn[i]-2./w[i]*(c[k]*h1[k]-h2[k]*c[i]
         +c[i+1]*h1[i]);
  }
  a[n]=fn[n]-2./w[n]*h1[n-1]*(c[n-1]-c[n]);

  for(i=0;i<=n-1;++i)
  {
    k=i+1;
    b[i]=h1[i]*(a[k]-a[i])-h[i]/3.*(c[k]+2.*c[i]);
    d[i]=h1[i]/3.*(c[k]-c[i]);
  }

  return 0;
}

// Compute the coefficients of a cubic fitting spline for given third
// derivatives at the end points.
int glsp3a(int n,REAL *xn,REAL *fn,REAL *w,
           REAL marg_0,REAL marg_n,
           REAL *a,REAL *b,REAL *c,REAL *d,
           REAL *ld1,REAL *ld2,REAL *ud1,REAL *ud2,
           REAL *h1,REAL *h)
{
  int i,k,error;
  REAL h_var_1,h_var_2;

  // compute aux values
  for(i=0;i<=n-1;++i)
  {
    h[i]=xn[i+1]-xn[i];
    h1[i]=1./h[i];
    c[i]=h1[i]*h1[i];
    b[i]=6./w[i];
  }
  b[n]=6./w[n];

  for(i=0;i<=n-2;++i)
    d[i]=h1[i]+h1[i+1];

  // Compute five diagonal system matrix A and the right hand side RS
  // of the system  A * C = RS
  for(i=3;i<=n-1;++i)
  {
    ld2[i]=b[i-1]*h1[i-2]*h1[i-1];  /* 2nd co-diagonals   */
    ud2[i-2]=ld2[i];
  }

  h_var_1=h[1]-b[2]*h1[1]*d[1];        /* first co-diagonals */
  ld1[2]=h_var_1-b[1]*c[1];
  ud1[1]=h_var_1-b[1]*h1[1]*d[0];
  for(i=3;i<=n-2;++i)
  {
    k=i-1;
    ld1[i]=h[k]-b[k]*h1[k]*d[k-1]-b[i]*h1[k]*d[k];
    ud1[k]=ld1[i];
  }
  h_var_1=h[n-2]-b[n-2]*h1[n-2]*d[n-3];
  ld1[n-1]=h_var_1-b[n-1]*h1[n-2]*d[n-2];
  ud1[n-2]=h_var_1-b[n-1]*c[n-2];

  a[1]=3.*h[0]+2.*h[1]+b[1]*h1[1]*d[0]+b[2]*c[1];
  for(i=2;i<=n-2;++i)
  {
    k=i-1;                                      /*  main diagonal */
    a[i]=2.*(h[k]+h[i])+b[k]*c[k]
         +b[i]*d[k]*d[k]+b[i+1]*c[i];
  }
  a[n-1]=3.*h[n-1]+2.*h[n-2]+b[n-1]*h1[n-2]*d[n-2]+b[n-2]*c[n-2];
  c[0]=0.5*marg_0;                           /* right hand side  */
  c[n]=0.5*marg_n;

  h_var_2=(fn[2]-fn[1])*h1[1];
  h_var_1=(fn[3]-fn[2])*h1[2];
  c[1]=3.*(h_var_2-(fn[1]-fn[0])*h1[0])
       +c[0]*(h[0]*h[0]-b[0]*h1[0]-b[1]*d[0]);
  c[2]=3.*(h_var_1-h_var_2)+c[0]*b[1]*h1[1];
  for(i=3;i<=n-3;++i,h_var_1=h_var_2)
  {
    h_var_2=(fn[i+1]-fn[i])*h1[i];
    c[i] =3.*(h_var_2-h_var_1);
  }
  h_var_2=(fn[n-1]-fn[n-2])*h1[n-2];
  c[n-2]=3.*(h_var_2-h_var_1)
         -c[n]*b[n-1]*h1[n-2];
  c[n-1]=3.*((fn[n]-fn[n-1])*h1[n-1]-h_var_2)
         -c[n]*(h[n-1]*h[n-1]-b[n-1]*d[n-2]-b[n]*h1[n-1]);
  // Compute coefficients c[i], i=1, ..., n-1
  // by LU factorization in diag5
  error=diag5(0,n-1,ld2+1,ld1+1,a+1,ud1+1,ud2+1,c+1);
  if (error!=0)
  {
    if (error==2)
      return 1;
    else
      return 2;
  }

  c[0]=c[1]-c[0]*h[0];
  c[n]=c[n-1]+c[n]*h[n-1];

  // Compute remaining spline coefficients
  a[0]=fn[0]+b[0]/3.*h1[0]*(c[0]-c[1]);
  for(i=1;i<=n-1;++i)
  {
    k=i-1;
    a[i]=fn[i]-b[i]/3.*(c[k]*h1[k]-d[k]*c[i]
         +c[i+1]*h1[i]);
  }
  a[n]=fn[n]-b[n]/3.*h1[n-1]*(c[n-1]-c[n]);

  for(i=0;i<=n-1;++i)
  {
    k=i+1;
    b[i]=h1[i]*(a[k]-a[i])-h[i]/3.*(c[k]+2.*c[i]);
    d[i]=h1[i]/3.*(c[k]-c[i]);
  }
  return 0;
}

// Computes the coefficients of a periodic cubic fitting spline.
int glsppe(int n,REAL *xn,REAL *fn,REAL *w,int rep,
           REAL *a,REAL *b,REAL *c,REAL *d,
           REAL *h,REAL *h1,REAL *h2,REAL *h3,
           REAL *rs,REAL *hup)
{ 
  int i,k,error;
  REAL h_var_1,h_var_2;

  if(rep!=0&&rep!=1) return 3;

  // Check periodicity
  //if(fn[n]!=fn[0]||w[n]!=w[0]) return 4;

  if(!rep)
  // First call : i.e., we must determine the aux variables and the
  // entries of the symmetric almost cyclic five-diagonal system matrix.
  {
    for(i=0;i<=n-1;++i)                         /*  aux variables  */
    {
      h[i]=xn[i+1]-xn[i];
      h1[i]=1./h[i];
      c[i]=h1[i]*h1[i];
      h2[i]=6./w[i];
    }
    h[n]=h[0];h1[n]=h1[0];c[n]=c[0];h2[n]=h2[0];
    for(i=0;i<=n-1;++i)
      h3[i]=h1[i]+h1[i+1];

    for(i=1;i<=n-1;++i)                    /*  second co-diagonal  */
      d[i]=h2[i+1]*h1[i]*h1[i+1];
    d[n]=h2[1]*h1[0]*h1[1];

    for(i=1;i<=n-1;++i)                      /* first co-diagonal  */
      b[i]=h[i]-h2[i]*h1[i]*h3[i-1]-h2[i+1]*h1[i]*h3[i];
    b[n]=h[0]-h2[0]*h1[0]*h3[n-1]-h2[1]*h1[0]*h3[0];

    for(i=1;i<=n-1;++i)                        /*  main diagonal   */
    {
      k=i-1;
      a[i]=2.*(h[k]+h[i])+h2[k]*c[k]
           +h2[i]*h3[k]*h3[k]+h2[i+1]*c[i];
    }
    a[n]=2.*(h[n-1]+h[n])+h2[n-1]*c[n-1]+h2[n]*h3[n-1]*h3[n-1]+h2[1]*c[0];
  }
  // right hand side
  h_var_1=(fn[1]-fn[0])*h1[0];
  for(i=1;i<=n-1;++i,h_var_1=h_var_2)
  {
    h_var_2=(fn[i+1]-fn[i])*h1[i];
    rs[i]=3.*(h_var_2-h_var_1);
  }
  rs[n]=3.*((fn[1]-fn[0])*h1[0]-h_var_1);
  // Find coefficients  c[i], i=0, ..., n
  // by solving linear system
  if(!rep)
  {
    error=fzyfsy(n,a,b,d,rs,c,            /* LU factorization */
                 &hup[0],&hup[n],&hup[2*n],
                 &hup[3*n],&hup[4*n-2],&hup[5*n-5],
                 &hup[6*n-5],&hup[7*n-5],&hup[8*n-9]);
            /* the 0th entry of the vector denotes the last entry of
               previous vector. This is ok since we do not use it here.
                                                                      */
    if(error!=0) return error;
  }
  else
    fzyfsl(n,rs,c,&hup[0],&hup[n]  ,/* for repeated calls no     */
           &hup[2*n],&hup[3*n],        /* factorization !           */
           &hup[4*n-2],&hup[5*n-5],
           &hup[6*n-5],&hup[7*n-5],&hup[8*n-9]);
                                                /* see call of fzyfsy */
  c[0]=c[n];
  // Compute remaining spline coefficients
  a[0]=fn[0]-h2[0]/3.*h1[0]*(c[1]-c[0])
       +h2[n]/3.*h1[n-1]*(c[n]-c[n-1]);
  for(i=1;i<=n-1;++i)
  {
    k=i-1;
    a[i]=fn[i]-h2[i]/3.*(c[k]*h1[k]-h3[k]*c[i]+c[i+1]*h1[i]);
  }
  a[n]=a[0];

  for(i=0;i<=n-1;++i)
  {
    k=i+1;
    b[i]=h1[i]*(a[k]-a[i])-h[i]/3.*(c[k]+2.*c[i]);
    d[i]=h1[i]/3.*(c[k]-c[i]);
  }
  return 0;
}

// Compute the coefficients of a nonparametric cubic fitting spline
CUBIC_SPLINE CreateCubicFittingSpline(int n,REAL *xn,REAL *fn,REAL *w,
            int marg_cond,REAL marg_0,REAL marg_n)
{
  int error=7,i;
  REAL *h=NULL,*h1=NULL,*h2=NULL,*md=NULL,
       *ud1=NULL,*ud2=NULL,*rs=NULL,*hup=NULL;

  CUBIC_SPLINE spline;


  spline.n=n;
  spline.a=(REAL*)calloc(n,sizeof(REAL));
  spline.b=(REAL*)calloc(n,sizeof(REAL));
  spline.c=(REAL*)calloc(n,sizeof(REAL));
  spline.d=(REAL*)calloc(n,sizeof(REAL));

  spline.x=(REAL*)calloc(n,sizeof(REAL));
  spline.y=(REAL*)calloc(n,sizeof(REAL));

  for(i=0;i<n;i++)
  {
    spline.x[i]=xn[i];
    spline.y[i]=fn[i];
  }

  /* check input data */

  switch (marg_cond)
  {
    case FIRST_DERIVATIVE:
    case SECOND_DERIVATIVE:
    case THIRD_DERIVATIVE:
      h=(REAL*)calloc(n,sizeof(REAL));
      h1=(REAL*)calloc(n,sizeof(REAL));
      h2=(REAL*)calloc(n,sizeof(REAL));
      md=(REAL*)calloc(n,sizeof(REAL));
      ud1=(REAL*)calloc(n,sizeof(REAL));
      ud2=(REAL*)calloc(n,sizeof(REAL));
      break;
    case PERIODIC_SPLINE:
      h=(REAL*)calloc(n+1,sizeof(REAL));
      h1=(REAL*)calloc(n+1,sizeof(REAL));
      h2=(REAL*)calloc(n+1,sizeof(REAL));
      md=(REAL*)calloc(n+1,sizeof(REAL));
      rs=(REAL*)calloc(n+1,sizeof(REAL));
      hup=(REAL*)calloc(9*n-11,sizeof(REAL));
      break;
    default:
      break;
  }

 /* call subroutines, depending on type of end point conditions  */
  switch (marg_cond)
  {
    case FIRST_DERIVATIVE:
      error=glsp1a(n,xn,fn,w,marg_0,marg_n,0,spline.a,spline.b,spline.c,spline.d,h,h1,h2,md,ud1,ud2);
      break;
    case SECOND_DERIVATIVE:
      error=glsp2a(n,xn,fn,w,marg_0,marg_n,0,spline.a,spline.b,spline.c,spline.d,h,h1,h2,md,ud1,ud2);
      break;
    case THIRD_DERIVATIVE:
      error=glsp3a(n,xn,fn,w,marg_0,marg_n,spline.a,spline.b,spline.c,spline.d,h2,md,ud1,ud2,h,h1);
      break;
    case PERIODIC_SPLINE:
      if(n<6)
        error=-1;
      else
        error=glsppe(n,xn,fn,w,0,spline.a,spline.b,spline.c,spline.d,h,h1,h2,md,rs,hup);
      break;
  }
  free(h);
  free(h1);
  free(h2);
  free(md);
  free(ud1);
  free(ud2);
  return spline;
}
