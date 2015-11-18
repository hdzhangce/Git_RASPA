/*****************************************************************************************************
    ewald_phonon.h -  description
    -----------------------------
    Copyright © 2006-2011 David Dubbeldam, Sofia Calero, Donald E. Ellis 
                          and Randall Q. Snurr. All Rights Reserved.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/depa/webdex/quimfis/miembros/Web_Sofia/Sofia.htm
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/
 *****************************************************************************************************/

#ifndef EWALD_PHONON_H
#define EWALD_PHONON_H

#include "simulation.h"
#include "vector.h"
#include "utils.h"

int CalculateEwaldFourierPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian);

#endif
