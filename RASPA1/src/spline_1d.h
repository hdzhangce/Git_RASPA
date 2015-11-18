/*****************************************************************************************************
    spline_1d.h -  description
    --------------------------
    Copyright © 2006-2011 David Dubbeldam, Sofia Calero, Donald E. Ellis 
                          and Randall Q. Snurr. All Rights Reserved.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/depa/webdex/quimfis/miembros/Web_Sofia/Sofia.htm
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/
 *****************************************************************************************************/


#ifndef SPLINE_1D_H
#define SPLINE_1D_H

#include "geometry.h"

REAL DotProdukt(POINT a,POINT b);


void PointOnBezierCurve(CPOINT *p,int n,REAL u,POINT *c);
POINT CurvePoint(int n,int p,REAL *U,POINT *P,REAL u);
void CurvePointDerivatives(int n,int p,REAL *U,POINT *P,REAL u,int d,POINT *CK);
void GlobalCurveInterpolation(int n,POINT *Q,int p,int m,REAL *U,CPOINT *P);

void CreateParametricBSpline(POLYGON *Q,int p,BSPLINE_CURVE *b);
void CreateParametricBSplineWithEndDerivatives(POLYGON *Q,int p,BSPLINE_CURVE *b,POINT start,POINT end);
POINT EvaluateParametricBSpline(BSPLINE_CURVE *b,REAL u);
void EvaluateParametricBSplineDerivatives(BSPLINE_CURVE *b,REAL u,int d,POINT *CK);
REAL MapPointToBspline(BSPLINE_CURVE *b,POINT P);
REAL MapPointToBsplinePeriodic(BSPLINE_CURVE *s,POINT P,POINT box);

#endif
