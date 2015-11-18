/*****************************************************************************************************
    spline_1d.c -  description
    --------------------------
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
#include "spline_1d.h"
#include "utils.h"
#include "matrix.h"
#include "linear_equations.h"

#define MAX_DEGREE 50
#define SQR(x) ((x)*(x))

REAL DotProdukt(POINT a,POINT b)
{
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}

// Bezier-curves
// ================================================================

void Horner(POINT *a,int n,REAL u,POINT *c)
{
  int i;

  c->x=a[n].x;
  c->y=a[n].y;
  c->z=a[n].z;
  for(i=n-1;i>=0;i--)
  {
    c->x=c->x*u+a[i].x;
    c->y=c->y*u+a[i].y;
    c->z=c->z*u+a[i].z;
  }
}

void ComputeAllBernsteinPolynomials(int n,REAL u,REAL *B)
{
  int j,k;
  REAL u1,saved,temp;

  B[0]=1.0;
  u1=1.0-u;
  for(j=1;j<=n;j++)
  {
    saved=0.0;
    for(k=0;k<j;k++)
    {
      temp=B[k];
      B[k]=saved+u1*temp;
      saved=u*temp;
    }
    B[j]=saved;
  } 
}

void PointOnBezierCurve(CPOINT *p,int n,REAL u,POINT *c)
{
  int k;
  REAL B[MAX_DEGREE];

  ComputeAllBernsteinPolynomials(n,u,B);
  c->x=0.0;
  c->y=0.0;
  c->z=0.0;
  for(k=0;k<=n;k++)
  {
    c->x+=B[k]*p[k].x;
    c->y+=B[k]*p[k].y;
    c->z+=B[k]*p[k].z;
  }
}

// one-dimensional non-uniform parametric BSplines
// ==================================================================================

int FindSpan(int n,int p,REAL u,REAL *U)
{
  int low,mid,high;

  if(u==U[n+1]) return n;
  low=p;
  high=n+1;
  mid=(low+high)/2;
  while((u<U[mid])||(u>=U[mid+1]))
  {
    if(u<U[mid]) high=mid;
    else low=mid;
    mid=(low+high)/2;
  }
  return mid;
}

void BasisFunctions(int i,REAL u,int p,REAL *U,REAL *N)
{
  int j,r;
  REAL saved,temp;
  REAL left[MAX_DEGREE],right[MAX_DEGREE];

  N[0]=1.0;
  for(j=1;j<=p;j++)
  {
    left[j]=u-U[i+1-j];
    right[j]=U[i+j]-u;
    saved=0.0;
    for(r=0;r<j;r++)
    {
      temp=N[r]/(right[r+1]+left[j-r]);
      N[r]=saved+right[r+1]*temp;
      saved=left[j-r]*temp;
    }
    N[j]=saved;
  }
}

void BasisFunctionsDerivatives(int i,REAL u,int p,int n,REAL *U,REAL **ders)
{
  int k,j,r,s1,s2,rk,pk,j1,j2;
  REAL saved,temp,d;
  REAL left[MAX_DEGREE],right[MAX_DEGREE];
  REAL ndu[MAX_DEGREE][MAX_DEGREE],a[MAX_DEGREE][MAX_DEGREE];

  ndu[0][0]=1.0;

  for(j=1;j<=p;j++)
  {
    left[j]=u-U[i+1-j];
    right[j]=U[i+j]-u;
    saved=0.0;
    for(r=0;r<j;r++)
    {
      ndu[j][r]=right[r+1]+left[j-r];  // lower triangle
      temp=ndu[r][j-1]/ndu[j][r];
      ndu[r][j]=saved+right[r+1]*temp; // upper triangle
      saved=left[j-r]*temp;
    }
    ndu[j][j]=saved;
  }

  for(j=0;j<=p;j++)
    ders[0][j]=ndu[j][p];

  for(r=0;r<=p;r++)
  {
    s1=0;
    s2=1;
    a[0][0]=1.0;
    for(k=1;k<=n;k++)
    {
      d=0.0;
      rk=r-k;
      pk=p-k;
      if(r>=k)
      {
        a[s2][0]=a[s1][0]/ndu[pk+1][rk];
        d=a[s2][0]*ndu[rk][pk];
      }
      if(rk>=-1) j1=1;
      else j1=-rk;
      if(r-1<=pk) j2=k-1;
      else j2=p-r;
      for(j=j1;j<=j2;j++)
      {
        a[s2][j]=(a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
        d+=a[s2][j]*ndu[rk+j][pk];
      }
      if(r<=pk)
      {
        a[s2][k]=-a[s1][k-1]/ndu[pk+1][r];
        d+=a[s2][k]*ndu[r][pk];
      }
      ders[k][r]=d;
      j=s1;
      s1=s2;
      s2=j;
    } 
  }
  r=p;
  for(k=1;k<=n;k++)
  {
    for(j=0;j<=p;j++)
      ders[k][j]*=r;
    r*=(p-k);
  }
}

POINT EvaluateParametricBSpline(BSPLINE_CURVE *b,REAL u)
{
  int i,span;
  POINT c;
  CPOINT *P;
  REAL N[MAX_DEGREE];

  span=FindSpan(b->pol.n,b->p,u,b->knt.U);
  BasisFunctions(span,u,b->p,b->knt.U,N);
  c.x=0.0;
  c.y=0.0;
  c.z=0.0;
  P=b->pol.Pw;
  for(i=0;i<=b->p;i++)
  {
    c.x+=N[i]*P[span-b->p+i].x;
    c.y+=N[i]*P[span-b->p+i].y;
    c.z+=N[i]*P[span-b->p+i].z;
  }
  return c;
}


void EvaluateParametricBSplineDerivatives(BSPLINE_CURVE *b,REAL u,int d,POINT *CK)
{
  int i,j,du,k,span;
  REAL **nders;
  CPOINT *P;

  nders=(REAL**)calloc(b->pol.n+1,sizeof(REAL*));
  for(i=0;i<=b->pol.n;i++)
    nders[i]=(REAL*)calloc(b->p+1,sizeof(REAL));
  du=d<b->p?d:b->p;
  for(k=b->p+1;k<=d;k++)
  {
    CK[k].x=0.0;
    CK[k].y=0.0;
    CK[k].z=0.0;
  }
  span=FindSpan(b->pol.n,b->p,u,b->knt.U);
  BasisFunctionsDerivatives(span,u,b->p,du,b->knt.U,nders);
  P=b->pol.Pw;
  for(k=0;k<=du;k++)
  {
    CK[k].x=0.0;
    CK[k].y=0.0;
    CK[k].z=0.0;
    for(j=0;j<=b->p;j++)
    {
      CK[k].x+=nders[k][j]*P[span-b->p+j].x;
      CK[k].y+=nders[k][j]*P[span-b->p+j].y;
      CK[k].z+=nders[k][j]*P[span-b->p+j].z;
    }
  } 
  for(i=0;i<=b->pol.n;i++)
    free(nders[i]);
  free(nders);
}

void GlobalCurveInterpolation(int n,POINT *Q,int p,int m,REAL *U,CPOINT *P)
{
  int i,j;
  int span;
  REAL d,sum;
  REAL u[MAX_DEGREE],*temp;
  REAL_FORTRAN_MATRIX A,B;

  temp=(REAL*)calloc(n+1,sizeof(REAL));

  d=0.0;
  for(i=1;i<=n;i++)
    d+=fabs(sqrt(SQR(Q[i].x-Q[i-1].x)+SQR(Q[i].y-Q[i-1].y)+SQR(Q[i].z-Q[i-1].z)));
  u[0]=0.0;
  for(i=1;i<n;i++)
    u[i]=u[i-1]+fabs(sqrt(SQR(Q[i].x-Q[i-1].x)+SQR(Q[i].y-Q[i-1].y)+SQR(Q[i].z-Q[i-1].z)))/d;
  u[n]=1.0;

  for(j=0;j<=p;j++)
    U[j]=0.0;
  for(j=1;j<=n-p;j++)
  {
    sum=0.0;
    for(i=j;i<=j+p-1;i++)
      sum+=u[i];
    U[j+p]=sum/(REAL)p;
  }
  for(j=m-p;j<=m;j++)
    U[j]=1.0;

  A=CreateRealFortranMatrix(n+1,n+1);
  for(i=0;i<=n;i++)
  {
    span=FindSpan(n,p,u[i],U);
    BasisFunctions(span,u[i],p,U,temp);
    for(j=0;j<=n;j++)
      A.element[i+(n+1)*(j+span-p)]=(REAL)temp[j];
  }
  PrintRealFortranMatrix(&A);
  B=CreateRealFortranMatrix(n+1,3);
  for(i=0;i<=n;i++)
  {
    B.element[i]=(REAL)Q[i].x;
    B.element[i+1*(n+1)]=(REAL)Q[i].y;
    B.element[i+2*(n+1)]=(REAL)Q[i].z;
  }
  SolveLinearSystem(&A,&B);
  for(i=0;i<=n;i++)
  {
    P[i].x=(REAL)B.element[i];
    P[i].y=(REAL)B.element[i+1*(n+1)];
    P[i].z=(REAL)B.element[i+2*(n+1)];
  }
  DeleteRealFortranMatrix(A);
  DeleteRealFortranMatrix(B);
  free(temp);
}

void GlobalCurveInterpolationWithEndDerivatives(int n,POINT *Q,int p,int m,REAL *U,CPOINT *P,POINT start,POINT end)
{
  int i,j;
  int span;
  REAL d,sum;
  REAL u[MAX_DEGREE],*temp;
  REAL_FORTRAN_MATRIX A,B;

  temp=(REAL*)calloc(n+3,sizeof(REAL));

  d=0.0;
  for(i=1;i<=n;i++)
    d+=fabs(sqrt(SQR(Q[i].x-Q[i-1].x)+SQR(Q[i].y-Q[i-1].y)+SQR(Q[i].z-Q[i-1].z)));
  u[0]=0.0;
  for(i=1;i<n;i++)
    u[i]=u[i-1]+fabs(sqrt(SQR(Q[i].x-Q[i-1].x)+SQR(Q[i].y-Q[i-1].y)+SQR(Q[i].z-Q[i-1].z)))/d;
  u[n]=1.0;

  for(j=0;j<=p;j++)
    U[j]=0.0;
  for(j=0;j<=n-p+1;j++)
  {
    sum=0.0;
    for(i=j;i<=j+p-1;i++)
      sum+=u[i];
    U[j+p+1]=sum/(REAL)p;
  }
  for(j=m-p;j<=m;j++)
    U[j]=1.0;

  A=CreateRealFortranMatrix(n+3,n+3);
  for(i=0;i<=n;i++)
  {
    span=FindSpan(n+2,p,u[i],U);
    BasisFunctions(span,u[i],p,U,temp);
    for(j=0;j<=n+2;j++)
      A.element[i+(n+3)*(j+span-p)]=temp[j];
  }
  A.element[n+1+(n+3)*(0)]=-1.0;
  A.element[n+1+(n+3)*(1)]=1.0;
  A.element[n+2+(n+3)*(n+1)]=-1.0;
  A.element[n+2+(n+3)*(n+2)]=1.0;
  //PrintRealFortranMatrix(&A);
  //exit(0);

  B=CreateRealFortranMatrix(n+3,3);

  for(i=0;i<=n;i++)
  {
    B.element[i]=Q[i].x;
    B.element[i+1*(n+3)]=Q[i].y;
    B.element[i+2*(n+3)]=Q[i].z;
  }
  B.element[n+1]=(U[p+1]/(REAL)p)*start.x;
  B.element[n+1+1*(n+3)]=(U[p+1]/(REAL)p)*start.y;
  B.element[n+1+2*(n+3)]=(U[p+1]/(REAL)p)*start.z;

  B.element[n+2]=((1.0-U[m-p-1])/(REAL)p)*end.x;
  B.element[n+2+1*(n+3)]=((1.0-U[m-p-1])/(REAL)p)*end.y;
  B.element[n+2+2*(n+3)]=((1.0-U[m-p-1])/(REAL)p)*end.z;

  printf("Solving....\n");
  SolveLinearSystem(&A,&B);
  printf("Solve done\n");
  for(i=0;i<=n+2;i++)
  {
    P[i].x=B.element[i];
    P[i].y=B.element[i+1*(n+3)];
    P[i].z=B.element[i+2*(n+3)];
  }
  DeleteRealFortranMatrix(A);
  DeleteRealFortranMatrix(B);
  free(temp);
  printf("finished\n");
}

void CreateParametricBSpline(POLYGON *Q,int p,BSPLINE_CURVE *b)
{
  int n,m;

  n=Q->n;
  m=n+p+1;  
  b->pol.n=n;
  b->pol.Pw=(CPOINT*)calloc(n+1,sizeof(CPOINT));
  b->p=p;
  b->knt.m=m;
  b->knt.U=(REAL*)calloc(m+1,sizeof(REAL));

  GlobalCurveInterpolation(n,Q->Pw,p,m,b->knt.U,b->pol.Pw);
}

void CreateParametricBSplineWithEndDerivatives(POLYGON *Q,int p,BSPLINE_CURVE *b,POINT start,POINT end)
{
  int n,m;

  n=Q->n+2;
  m=n+p+1;
  b->pol.n=n;
  b->pol.Pw=(CPOINT*)calloc(n+1,sizeof(CPOINT));
  b->p=p;
  b->knt.m=m;
  b->knt.U=(REAL*)calloc(m+1,sizeof(REAL));
  GlobalCurveInterpolationWithEndDerivatives(n-2,Q->Pw,p,m,b->knt.U,b->pol.Pw,start,end);
}


REAL MapPointToBspline(BSPLINE_CURVE *s,POINT P)
{
  int i,n;
  int condition_1,condition_2,condition_4;
  REAL u_new,u,distance,minimal_distance,u_minimal;
  REAL f,fderiv,a,b;
  POINT c,CK[4];

  // search for minimal distance in 100 steps as u0 starting value
  n=100; 
  u_minimal=0.5;
  minimal_distance=100000.0;
  for(i=0;i<=n;i++)
  {
    u=(REAL)i/(REAL)n;
    c=EvaluateParametricBSpline(s,u);
    distance=sqrt(SQR(P.x-c.x)+SQR(P.y-c.y)+SQR(P.z-c.z));
    if(distance<minimal_distance) 
    {
      minimal_distance=distance;
      u_minimal=u;
    }
  }
  u=u_minimal;
  a=u_minimal-1.0/(REAL)n;
  if(a<0.0) a=0.0;
  b=u_minimal+1.0/(REAL)n;
  if(b>1.0) b=1.0;
  do
  {
    EvaluateParametricBSplineDerivatives(s,u,3,CK);
    f=CK[1].x*(CK[0].x-P.x)+CK[1].y*(CK[0].y-P.y)+CK[1].z*(CK[0].z-P.z);
    fderiv=((CK[2].x*(CK[0].x-P.x)+CK[2].y*(CK[0].y-P.y)+CK[2].z*(CK[0].z-P.z))+
           (SQR(CK[1].x)+SQR(CK[1].y)+SQR(CK[1].z)));
    u_new=u-f/fderiv;

    // condition 1: point coincidence
    condition_1=sqrt(SQR(P.x-CK[0].x)+SQR(P.y-CK[0].y)+SQR(P.z-CK[0].z))<=(REAL)DBL_EPSILON;
    // condition 2: zero cosine
    condition_2=fabs(f)/(sqrt(SQR(CK[1].x)+SQR(CK[1].y)+SQR(CK[1].z))*
           sqrt(SQR(CK[0].x-P.x)+SQR(CK[0].y-P.y)+SQR(CK[0].z-P.z)))<=(REAL)DBL_EPSILON;
    // condition 3: clip to interval
    if(u_new<a) u_new=a;
    if(u_new>b) u_new=b;
    // condition 4: the parameter does not change significantly (end of the curve)
    condition_4=sqrt(SQR((u_new-u)*CK[1].x)+SQR((u_new-u)*CK[1].y)+SQR((u_new-u)*CK[1].z))<=(REAL)DBL_EPSILON;

    u=u_new;
  } while((!condition_1)&&(!condition_2)&&(!condition_4));
  //printf("[%f %f] %f\n",a,b,u);
  if(u<0.0) return 0.0;
  else if(u>1.0) return 1.0;
  else return u;
}

REAL MapPointToBsplinePeriodic(BSPLINE_CURVE *s,POINT P,POINT box)
{
  int i,n;
  REAL u,distance,minimal_distance,u_minimal;
  REAL dr[3];
  POINT c;

  // search for minimal distance in 100 steps as u0 starting value
  n=200;
  u_minimal=0.5;
  minimal_distance=100000.0;
  for(i=0;i<=n;i++)
  {
    u=(REAL)i/(REAL)n;
    c=EvaluateParametricBSpline(s,u);
    dr[0]=P.x-c.x;
    dr[1]=P.y-c.y;
    dr[2]=P.z-c.z;
    dr[0]-=box.x*NINT(dr[0]/box.x);
    dr[1]-=box.y*NINT(dr[1]/box.y);
    dr[2]-=box.z*NINT(dr[2]/box.z);
    distance=sqrt(SQR(dr[0])+SQR(dr[1])+SQR(dr[2]));
    if(distance<minimal_distance)
    {
      minimal_distance=distance;
      u_minimal=u;
    }
  }
  u=u_minimal;

/*
  a=u_minimal-1.0/(REAL)n;
  if(a<0.0) a=0.0;
  b=u_minimal+1.0/(REAL)n;
  if(b>1.0) b=1.0;
  do
  {
    EvaluateParametricBSplineDerivatives(s,u,3,CK);
    f=CK[1].x*(CK[0].x-P.x)+CK[1].y*(CK[0].y-P.y)+CK[1].z*(CK[0].z-P.z);
    fderiv=((CK[2].x*(CK[0].x-P.x)+CK[2].y*(CK[0].y-P.y)+CK[2].z*(CK[0].z-P.z))+
           (SQR(CK[1].x)+SQR(CK[1].y)+SQR(CK[1].z)));
    u_new=u-f/fderiv;

    // condition 1: point coincidence
    condition_1=sqrt(SQR(P.x-CK[0].x)+SQR(P.y-CK[0].y)+SQR(P.z-CK[0].z))<=DBL_EPSILON;
    // condition 2: zero cosine
    condition_2=fabs(f)/(sqrt(SQR(CK[1].x)+SQR(CK[1].y)+SQR(CK[1].z))*
           sqrt(SQR(CK[0].x-P.x)+SQR(CK[0].y-P.y)+SQR(CK[0].z-P.z)))<=DBL_EPSILON;
    // condition 3: clip to interval
    if(u_new<a) u_new=a;
    if(u_new>b) u_new=b;
    // condition 4: the parameter does not change significantly (end of the curve)
    condition_4=sqrt(SQR((u_new-u)*CK[1].x)+SQR((u_new-u)*CK[1].y)+SQR((u_new-u)*CK[1].z))<=DBL_EPSILON;

    u=u_new;
  } while((!condition_1)&&(!condition_2)&&(!condition_4));
  //printf("[%f %f] %f\n",a,b,u);
  if(u<0) return 0;
  else if(u>1.0) return 1.0;
  else return u;
*/
  return u;
}

