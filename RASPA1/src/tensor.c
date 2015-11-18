/*****************************************************************************************************
    tensor.c -  description
    -----------------------
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
#include "tensor.h"

REAL_TENSOR CreateRealTensor(int r,int c,int d)
{
  int i,j;
  REAL_TENSOR t;

  t.d=r;
  t.n=c;
  t.m=d;

  t.element=(REAL***)calloc(r,sizeof(REAL**));
  t.element[0]=(REAL**)calloc(r*c,sizeof(REAL*));
  t.element[0][0]=(REAL*)calloc(r*c*d,sizeof(REAL));

  for(j=1;j<c;j++) 
    t.element[0][j]=t.element[0][j-1]+d;
  for(i=1;i<r;i++) 
  {
    t.element[i]=t.element[i-1]+c;
    t.element[i][0]=t.element[i-1][0]+c*d;
    for(j=1;j<c;j++) 
      t.element[i][j]=t.element[i][j-1]+d;
  }
  return t;
}

void DeleteRealTensor(REAL_TENSOR t)
{
  free(t.element[0][0]);
  free(t.element[0]);
  free(t.element);
}
