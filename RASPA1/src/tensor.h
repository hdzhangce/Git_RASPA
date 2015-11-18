/*****************************************************************************************************
    tensor.h -  description
    -----------------------
    Copyright © 2006-2011 David Dubbeldam, Sofia Calero, Donald E. Ellis 
                          and Randall Q. Snurr. All Rights Reserved.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/depa/webdex/quimfis/miembros/Web_Sofia/Sofia.htm
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/
 *****************************************************************************************************/


#ifndef TENSOR_H
#define TENSOR_H

#include "constants.h"

typedef struct real_tensor
{
  int m;
  int n;
  int d;
  REAL ***element;
} REAL_TENSOR;

REAL_TENSOR CreateRealTensor(int r,int c,int d);
void DeleteRealTensor(REAL_TENSOR t);

#endif
