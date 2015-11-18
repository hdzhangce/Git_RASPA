/*****************************************************************************************************
    geometry.c -  description
    -------------------------
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
#include "geometry.h"
#include "utils.h"
#include "vector.h"


KNOTVECTOR CreateKnotVector(int m)
{
  KNOTVECTOR c;

  c.m=m;
  c.U=(REAL*)calloc(m+1,sizeof(REAL));
  return c;
}

void DeleteKnotVector(KNOTVECTOR c)
{
  free(c.U);
}

void PrintKnotVector(KNOTVECTOR c)
{
  int i;
  for(i=0;i<c.m;i++)
    printf("%f ",(double)c.U[i]);
  printf("\n");
}

CNET CreateControlNet(int m,int n)
{
  int i;
  CNET c;

  printf("allocing control-net %dx%d matrix\n",m,n);
  c.m=m;
  c.n=n;
  c.Pw=(CPOINT**)calloc(m,sizeof(CPOINT*));
  c.Pw[0]=(CPOINT*)calloc(m*n,sizeof(CPOINT));

  for(i=1;i<m;i++)
    c.Pw[i]=c.Pw[i-1]+n;
  return c;
}

void DeleteControlNet(CNET m)
{
  free(m.Pw[0]);
  free(m.Pw);
}

REAL Distance3D(POINT a,POINT b)
{
  return sqrt(a.x*b.x+a.y*b.y+a.z*b.z);
}

REAL ShortestDistancePointToLine(POINT point,LINE line)
{
  VECTOR v,w;
  REAL c1,c2,b;
  POINT Pb;

  v.x=line.p2.x-line.p1.x;
  v.y=line.p2.y-line.p1.y;
  v.z=line.p2.z-line.p1.z;

  w.x=point.x-line.p1.x;
  w.y=point.y-line.p1.y;
  w.z=point.z-line.p1.z;

  c1=DotProduct(w,v);
  c2=DotProduct(v,v);
  b=c1/c2;

  Pb.x=line.p1.x+b*v.x;
  Pb.y=line.p1.y+b*v.y;
  Pb.z=line.p1.z+b*v.z;

  return sqrt(SQR(point.x-Pb.x)+SQR(point.y-Pb.y)+SQR(point.z-Pb.z));
}

POINT PointToLine(POINT point,LINE line)
{
  VECTOR v,w;
  REAL c1,c2,b;
  POINT Pb;

  v.x=line.p2.x-line.p1.x;
  v.y=line.p2.y-line.p1.y;
  v.z=line.p2.z-line.p1.z;

  w.x=point.x-line.p1.x;
  w.y=point.y-line.p1.y;
  w.z=point.z-line.p1.z;

  c1=DotProduct(w,v);
  c2=DotProduct(v,v);
  b=c1/c2;

  Pb.x=line.p1.x+b*v.x;
  Pb.y=line.p1.y+b*v.y;
  Pb.z=line.p1.z+b*v.z;

  return Pb;
}


/*
REAL ShortestLineToLineDistance(LINE L1,LINE L2)
{
  REAL a,b,c,d,e,D,sc,tc;
  VETOR w;

  w.x=L1.p.x-L2.p.x;
  w.y=L1.p.y-L2.p.y;
  w.z=L1.p.z-L2.p.z;

    Vector   u = L1.P1 - L1.P0;
    Vector   v = L2.P1 - L2.P0;
    Vector   w = L1.P0 - L2.P0;

    float    a = dot(u,u);        // always >= 0
    float    b = dot(u,v);
    float    c = dot(v,v);        // always >= 0
    float    d = dot(u,w);
    float    e = dot(v,w);
    float    D = a*c - b*b;       // always >= 0
    float    sc, tc;

  a=DotProduct(L1.v,L1.v);
  b=DotProduct(L1.v,L2.v);

   
    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) {         // the lines are almost parallel
        sc = 0.0;
        tc = (b>c ? d/b : e/c);   // use the largest denominator
    }
    else {
        sc = (b*e - c*d) / D;
        tc = (a*e - b*d) / D;
    }

    // get the difference of the two closest points
    Vector   dP = w + (sc * u) - (tc * v);  // = L1(sc) - L2(tc)

    return norm(dP);   // return the closest distance
}
*/

