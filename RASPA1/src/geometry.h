/*****************************************************************************************************
    geomtry.h -  description
    ------------------------
    Copyright © 2006-2011 David Dubbeldam, Sofia Calero, Donald E. Ellis 
                          and Randall Q. Snurr. All Rights Reserved.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/depa/webdex/quimfis/miembros/Web_Sofia/Sofia.htm
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/
 *****************************************************************************************************/


#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "complex.h"

typedef short DEGREE;
typedef int INDEX;

// a general point in three-dimensional space
typedef struct point
{
  REAL x;
  REAL y;
  REAL z;
} POINT,VECTOR;

typedef struct complex_point
{
  COMPLEX x;
  COMPLEX y;
  COMPLEX z;
} CVECTOR;


// a control point
typedef struct cpoint
{
  REAL x;
  REAL y;
  REAL z;
} CPOINT;

// a polygon
typedef struct polygon
{
  int n;
  POINT *Pw;
} POLYGON;

// a control polygon
typedef struct cpolygon
{
  int n;
  CPOINT *Pw;
} CPOLYGON;

// a knotvector
typedef struct knotvector
{
  int m;
  REAL *U;
} KNOTVECTOR;

// a curve
typedef struct curve
{
  CPOLYGON pol;
  DEGREE p;
  KNOTVECTOR knt;
} BSPLINE_CURVE;

// a control net
typedef struct cnet
{
  INDEX m;
  INDEX n;
  CPOINT **Pw;
} CNET;

// a surface
typedef struct surface
{
  CNET net;
  DEGREE p;
  DEGREE q;
  KNOTVECTOR knu;
  KNOTVECTOR knv;
} BSPLINE_SURFACE;

// a line
typedef struct line
{
  POINT p1;
  POINT p2;
  VECTOR v;
} LINE;

// a plane
typedef struct plane
{
  POINT p;
  VECTOR n;
} PLANE;

// a circle
typedef struct circle
{
  POINT c;
  REAL r;
  VECTOR n;
} CIRCLE;

KNOTVECTOR CreateKnotVector(int m);
void DeleteKnotVector(KNOTVECTOR c);
void PrintKnotVector(KNOTVECTOR c);


CNET CreateControlNet(int m,int n);
void DeleteControlNet(CNET m);

REAL Distance3D(POINT a,POINT b);

REAL ShortestDistancePointToLine(POINT point,LINE line);
POINT PointToLine(POINT point,LINE line);

#endif
