/*****************************************************************************************************
    nonuniform_cubic_spline_1d.h -  description
    -------------------------------------------
    Copyright © 2006-2011 David Dubbeldam, Sofia Calero, Donald E. Ellis 
                          and Randall Q. Snurr. All Rights Reserved.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/depa/webdex/quimfis/miembros/Web_Sofia/Sofia.htm
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/
 *****************************************************************************************************/


#ifndef NONUNIFORM_CUBIC_SPLINE_1D_H
#define NONUNIFORM_CUBIC_SPLINE_1D_H

//#include "geometry.h"
#include "cubic_spline_1d.h"

CUBIC_SPLINE CreateCubicSpline(int m,REAL *x,REAL *y,int boundary_condition,
                      REAL left_boundary,REAL right_boundary);

CUBIC_SPLINE CreateCubicFittingSpline(int n,REAL *xn,REAL *fn,REAL *w,
              int marg_cond,REAL marg_0,REAL marg_n);
#endif
