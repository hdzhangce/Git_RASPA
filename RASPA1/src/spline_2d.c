/*****************************************************************************************************
    spline_2d.c -  description
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
#include "spline_2d.h"
#include "vector.h"
#include "matrix.h"
#include "linear_equations.h"

void CreateParametricSurface(POINT_MATRIX *Q,int p,int q,BSPLINE_SURFACE *b)
{
  int k,l;
  int m,n,num;
  double *uk,*ul,*cds,d,total;

  m=Q->m;
  n=Q->n;

  b->p=p;
  b->q=q;
  b->net=CreateControlNet(m,n);

  b->knu=CreateKnotVector(Q->m+p+1);
  b->knv=CreateKnotVector(Q->n+q+1);

  uk=(double *)calloc(m,sizeof(double));  
  ul=(double *)calloc(n,sizeof(double));
  cds=(double *)calloc(m>n?m+1:n+1,sizeof(double));  

  num=m+1;
  total=0.0;
  uk[0]=0.0;
  uk[n]=1.0;
  for(k=1;k<n;k++)
    uk[k]=0.0;
  for(l=0;l<n;l++)
  {
    total=0.0;
    for(k=1;k<m;k++)
    {
      cds[k]=Distance3D(Q->element[k][l],Q->element[k-1][l]);
      total+=cds[k];
    }
  } 
  if(total==0.0) num--;
  else
  {
    d=0.0;
    for(k=1;k<n;k++)
    {
      d+=cds[k];
      uk[k]+=d/total;
    }
  }
  for(k=1;k<n;k++)
    uk[k]/=num;

  num=n+1;
  total=0.0;
  ul[0]=0.0;
  ul[m]=1.0;
  for(k=1;k<m;k++)
    ul[k]=0.0;
  for(l=0;l<=n;l++)
  {
    total=0.0;
    for(k=1;k<m;k++)
    {
      cds[k]=Distance3D(Q->element[k][l],Q->element[k-1][l]);
      total+=cds[k];
    }
  }
  if(total==0.0) num--;
  else
  {
    d=0.0;
    for(k=1;k<m;k++)
    {
      d+=cds[k];
      ul[k]+=d/total;
    }
  }
  for(k=1;k<m;k++)
    ul[k]/=num;

  for(k=0;k<n;k++)
    printf("uk[%d]:%g\n",k,uk[k]);

  for(k=0;k<m;k++)
    printf("ul[%d]:%g\n",k,ul[k]);
}
