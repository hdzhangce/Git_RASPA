/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'grids.c' is part of RASPA.

    RASPA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RASPA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "potentials.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "framework_hessian.h"
#include "simulation.h"
#include "integration.h"
#include "cbmc.h"
#include "ewald.h"
#include "utils.h"
#include "grids.h"
#include "output.h"
#include "grids.h"
#include "movies.h"
#include "warnings.h"

// For a framework that is kept rigid it is effecient to precompute the energy and forces.
// The amount of points is compute from
//   (a) option 'SpacingVDWGrid 0.15' for Van der Waal grids
//   (b) option 'SpacingCoulombGrid 0.15' for the electrostatic grid
// 0.15 Angstrom per point seems to be a good value
// the sizes of the grids can be large, so due to memory-restrictions the amount and type of pseudo-atoms to
// be precomputed needs to chosen with care.

// for 2D bicubic interpolation, see Numerical Recipes 3.6 Interpolation in more dimensions
// the code here is the extension to 3D using the transformation matrix of Lekien and Marsden
//   F. Lakien and J. Marsden, "Tricubic interpolation in three dimensions"
//   International Journal for Numerical Methods in Engineering, 63, 455-471, 2005
// 
// The algorithme is a local tricubic interpolation scheme in three dimensions that is both C^1 and isotropic. It uses
// a 64x64 matrix that gives the relationship between the derivatives at the corners of the elements and the coefficients of
// the trcubic interpolant for this element.
// The problem is *not* separated into 3 three-dimensional problems. This allows for a much easier and accurate computation of
// higher derivatives, e.g. the force, of the extrapolated field.
//
// Grids for triclinic frameowrks are based on a enclosing rectangular box. The tricubic interpolation algorithm operates
// on a unit cube. Notice that the derivatives are given in terms of the unit cube and must be divided by the actual size
// of the elements.

static int Coeff[64][64] = {
{  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 , 0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  9, -9, -9,  9,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0,  6, -6,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   4,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ -6,  6,  6, -6,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0, -4,  4, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -2, -2, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ -6,  6,  6, -6,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0, -3,  3, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -2, -1, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  4, -4, -4,  4,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0,  2, -2,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  9, -9, -9,  9,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0,  6, -6,  3, -3,  0,  0,  0,  0,  4,  2,  2,  1,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0, -4,  4, -2,  2,  0,  0,  0,  0, -2, -2, -1, -1,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0, -3,  3, -3,  3,  0,  0,  0,  0, -2, -1, -2, -1,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -4, -4,  4,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0,  2, -2,  2, -2,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0},
{ -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  9, -9,  0,  0, -9,  9,  0,  0,  6,  3,  0,  0, -6, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  3, -3,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ -6,  6,  0,  0,  6, -6,  0,  0, -3, -3,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  4,  0,  0, -2,  2,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  9, -9,  0,  0, -9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   6,  3,  0,  0, -6, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  3, -3,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -3, -3,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  4,  0,  0, -2,  2,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0},
{  9,  0, -9,  0, -9,  0,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  2,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  9,  0, -9,  0, -9,  0,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  2,  0,  2,  0,  1,  0},
{-27, 27, 27,-27, 27,-27,-27, 27,-18, -9, 18,  9, 18,  9,-18, -9,-18, 18, -9,  9, 18,-18,  9, -9,-18, 18, 18,-18, -9,  9,  9, -9,
 -12, -6, -6, -3, 12,  6,  6,  3,-12, -6, 12,  6, -6, -3,  6,  3,-12, 12, -6,  6, -6,  6, -3,  3, -8, -4, -4, -2, -4, -2, -2, -1},
{ 18,-18,-18, 18,-18, 18, 18,-18,  9,  9, -9, -9, -9, -9,  9,  9, 12,-12,  6, -6,-12, 12, -6,  6, 12,-12,-12, 12,  6, -6, -6,  6,
   6,  6 , 3,  3, -6, -6, -3, -3,  6,  6, -6, -6,  3,  3, -3, -3,  8, -8,  4, -4,  4, -4,  2, -2,  4,  4,  2,  2,  2,  2,  1,  1},
{ -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -2,  0, -1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -2,  0, -1,  0, -1,  0},
{ 18,-18,-18, 18,-18, 18, 18,-18, 12,  6,-12, -6,-12, -6, 12,  6,  9, -9,  9, -9, -9,  9, -9,  9, 12,-12,-12, 12,  6, -6, -6,  6,
   6,  3,  6,  3, -6, -3, -6, -3,  8,  4, -8, -4,  4,  2, -4, -2,  6, -6,  6, -6,  3, -3,  3, -3,  4,  2,  4,  2,  2,  1,  2,  1},
{-12, 12, 12,-12, 12,-12,-12, 12, -6, -6,  6,  6,  6,  6, -6, -6, -6,  6, -6,  6,  6, -6,  6, -6, -8,  8,  8, -8, -4,  4,  4, -4,
  -3, -3, -3, -3,  3,  3,  3,  3, -4, -4,  4,  4, -2, -2,  2,  2, -4,  4, -4,  4, -2,  2, -2,  2, -2, -2, -2, -2, -1, -1, -1, -1},
{  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ -6,  6,  0,  0,  6, -6,  0,  0, -4, -2,  0,  0,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -3,  3,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  4, -4,  0,  0, -4,  4,  0,  0,  2,  2,  0,  0, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  2, -2,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -4, -2,  0,  0,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   2,  2,  0,  0, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0},
{ -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0, -2,  0, -1,  0},
{ 18,-18,-18, 18,-18, 18, 18,-18, 12,  6,-12, -6,-12, -6, 12,  6, 12,-12,  6, -6,-12, 12, -6,  6,  9, -9, -9,  9,  9, -9, -9,  9,
   8,  4,  4,  2, -8, -4, -4, -2,  6,  3, -6, -3,  6,  3, -6, -3,  6, -6,  3, -3,  6, -6,  3, -3,  4,  2,  2,  1,  4,  2,  2,  1},
{-12, 12, 12,-12, 12,-12,-12, 12, -6, -6,  6,  6,  6,  6, -6, -6, -8,  8, -4,  4,  8, -8,  4, -4, -6,  6,  6, -6, -6,  6,  6, -6,
  -4, -4, -2, -2 , 4,  4,  2,  2, -3, -3,  3,  3, -3, -3,  3 , 3, -4,  4, -2,  2, -4,  4, -2,  2, -2, -2, -1, -1, -2, -2, -1, -1},
{  4,  0, -4,  0, -4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{  0,  0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0, -4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0},
{-12, 12, 12,-12, 12,-12,-12, 12, -8, -4,  8,  4,  8,  4, -8, -4, -6,  6, -6,  6,  6, -6,  6, -6, -6,  6,  6, -6, -6,  6,  6, -6,
  -4, -2, -4, -2,  4,  2,  4,  2, -4, -2,  4,  2, -4, -2,  4,  2, -3,  3, -3,  3, -3,  3, -3,  3, -2, -1, -2, -1, -2, -1, -2, -1},
{  8, -8, -8,  8, -8,  8,  8, -8,  4,  4, -4, -4, -4, -4,  4,  4,  4, -4,  4, -4, -4,  4, -4,  4,  4, -4, -4,  4,  4, -4, -4,  4,
   2,  2,  2,  2, -2, -2, -2, -2,  2,  2, -2, -2,  2,  2, -2, -2,  2, -2,  2, -2,  2, -2,  2, -2,  1,  1,  1,  1,  1,  1,  1,  1}};

int UseTabularGrid;
int UseDynamicGrid;
int UseDelaunayGrid;
int EBCBMC;
int HermiteInterpolation;

REAL SpacingDelaunayGrid;
REAL SpacingVDWGrid;
REAL SpacingCoulombGrid;
REAL SpacingDynamicGrid;

int NumberOfGrids;

int *****DelaunayGrid;
int ****DynamicGrid;
float *****VDWGrid;
float ****CoulombGrid;

INT_VECTOR3 NumberOfVDWGridPoints;
INT_VECTOR3 NumberOfDelaunayGridPoints;
INT_VECTOR3 NumberOfDynamicGridPoints;
static INT_VECTOR3 NumberOfCoulombGridPoints;
static VECTOR SizeGrid;
static VECTOR ShiftGrid;
static VECTOR DeltaVDWGrid;
static VECTOR DeltaDelaunayGrid;
static VECTOR DeltaCoulombGrid;
static VECTOR DeltaDynamicGrid;
int *GridTypeList;
char (*GridTypeListName)[256];
char DelaunayPolyhedraFileName[256];
int DelaunayPolyhedraFile;

// Blocking variables
static char ***BlockingGrid;
int BlockEnergyGrids;
int BlockGridPockets;
int BlockGridPores;
REAL BlockEnergyGridOverlapCriteria;
int NumberOfGridSeeds;
VECTOR *GridSeeds;

static INT_VECTOR3 *Queue;
static int QueueSize;

// Delaunay variables
VECTOR *DelaunayUnitCellSize;
int *PolyNumberFaces, *PolyNumberVertices, **PolyNumVertPerFace, **PolyVertices, ***PolyVertPerFace,*PolyNumberEdges,***PolyEdges;
float ***PolyPlanes,*projectedpoint,*point;
VECTOR ***Poly2DPoints,***Plane2D,**PolyVertLocation,PointOnPlane;
INT_VECTOR3 **PolyDelsVertices;
int NumberPolyhedra;

//EBCBMC Variables
float ***EBCBMCValue;
INT_VECTOR3 ***EBCBMCGridLocation;
int **SizeOfEBCBMCBins;
REAL *SumVDWGridRosenbluth;
int NumberOfEBCBMCBins;
int EBCBMCVolume;
int UseEBCBMCEnergyFirstBead;



VECTOR MapToUnitCell(VECTOR pos)
{
  VECTOR s;

  switch(BoundaryCondition[CurrentSystem])
  {
    case RECTANGULAR:
      pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
      pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
      pos.z-=UnitCellSize[CurrentSystem].z*((REAL)NINT(pos.z/UnitCellSize[CurrentSystem].z));
      if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
      if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
      if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
      break;
    case TRICLINIC:
      s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
      s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
      s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;
      // apply boundary condition
      s.x-=(REAL)NINT(s.x);
      s.y-=(REAL)NINT(s.y);
      s.z-=(REAL)NINT(s.z);
      if(s.x<0.0) s.x+=1.0;
      if(s.y<0.0) s.y+=1.0;
      if(s.z<0.0) s.z+=1.0;
      pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
  }
  return pos;
}

int CheckIfPointIsInUnitCell(VECTOR pos)
{
  VECTOR s;

  s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
  s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
  s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

  if(s.x<0.0) return FALSE;
  if(s.y<0.0) return FALSE;
  if(s.z<0.0) return FALSE;
  if(s.x>=1.0) return FALSE;
  if(s.y>=1.0) return FALSE;
  if(s.z>=1.0) return FALSE;
  return TRUE;
}

VECTOR MapZToBox(VECTOR pos)
{
  pos.z-=Box[0].cz*NINT(pos.z/Box[0].cz);
  if(pos.z<0.0) pos.z+=Box[0].cz;
  return pos;
}

void MakeASCIGrid(void)
{
  int i,j,k,l,typeA,save;
  REAL percent,teller,TailEnergy;
  POINT pos;
  REAL value,third_derivative;
  VECTOR first_derivative;
  REAL_MATRIX3x3 second_derivative;
  char buffer[1024];
  FILE *FilePtr;

  CurrentSystem=0;
  if(MIN3(BoxProperties[CurrentSystem].cx,BoxProperties[CurrentSystem].cy,
          BoxProperties[CurrentSystem].cz)<CutOffVDW)
  {
     printf("ERROR:  Cutoff smaller than half of one of the perpendicular boxlengths !!!\n");
     printf("        (Cutoff: %lf perpendicular boxlengths: %lf %lf %lf)\n",(double)CutOffVDW,
             (double)BoxProperties[CurrentSystem].cx,(double)BoxProperties[CurrentSystem].cy,(double)BoxProperties[CurrentSystem].cz);
     printf("Advice: choose more unitcells for constructing the grid.\n");
     exit(0);
  }

  NumberOfAdsorbateMolecules[CurrentSystem]=0;
  NumberOfCationMolecules[CurrentSystem]=0;

  CalculateForce();
  TailEnergy=UTailCorrection[CurrentSystem];

  save=ChargeMethod;
  ChargeMethod=NONE;

  // compute the size of the grid and the shift of the origin of the enclosing box
  SizeGrid.x=SizeGrid.y=SizeGrid.z=0.0;
  ShiftGrid.x=ShiftGrid.y=ShiftGrid.z=0.0;

  SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].ax);
  SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].ay);
  SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].az);
  if(UnitCellBox[CurrentSystem].ax<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].ax;
  if(UnitCellBox[CurrentSystem].ay<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].ay;
  if(UnitCellBox[CurrentSystem].az<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].az;

  SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].bx);
  SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].by);
  SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].bz);
  if(UnitCellBox[CurrentSystem].bx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].bx;
  if(UnitCellBox[CurrentSystem].by<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].by;
  if(UnitCellBox[CurrentSystem].bz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].bz;

  SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].cx);
  SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].cy);
  SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].cz);
  if(UnitCellBox[CurrentSystem].cx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].cx;
  if(UnitCellBox[CurrentSystem].cy<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].cy;
  if(UnitCellBox[CurrentSystem].cz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].cz;

  // compute the number of grid points
  NumberOfVDWGridPoints.x=(int)(SizeGrid.x/SpacingVDWGrid);
  NumberOfVDWGridPoints.y=(int)(SizeGrid.y/SpacingVDWGrid);
  NumberOfVDWGridPoints.z=(int)(SizeGrid.z/SpacingVDWGrid);

  DeltaVDWGrid.x=SizeGrid.x/(REAL)NumberOfVDWGridPoints.x;
  DeltaVDWGrid.y=SizeGrid.y/(REAL)NumberOfVDWGridPoints.y;
  DeltaVDWGrid.z=SizeGrid.z/(REAL)NumberOfVDWGridPoints.z;

  printf("ShiftGrid: %g %g %g\n",ShiftGrid.x,ShiftGrid.y,ShiftGrid.z);
  printf("SizeGrid: %g %g %g\n",SizeGrid.x,SizeGrid.y,SizeGrid.z);
  printf("Number of grid points: %d %d %d\n",NumberOfVDWGridPoints.x,NumberOfVDWGridPoints.y,NumberOfVDWGridPoints.z);

  percent=100.0/(REAL)((NumberOfVDWGridPoints.x+1)*NumberOfGrids);
  teller=0.0;

  fprintf(OutputFilePtr[CurrentSystem],"\n\n");
  fprintf(OutputFilePtr[CurrentSystem],"Generating an ASCI interpolation grid (%d x %d x %d)\n",NumberOfVDWGridPoints.x,NumberOfVDWGridPoints.y,NumberOfVDWGridPoints.z);
  fprintf(OutputFilePtr[CurrentSystem],"========================================================\n\n");
  fflush(OutputFilePtr[CurrentSystem]);

  mkdir("ASCI_Grids",S_IRWXU);

  for(l=0;l<NumberOfGrids;l++)
  {
    fprintf(OutputFilePtr[CurrentSystem],"Creating grid %d [%s]\n",l,PseudoAtoms[GridTypeList[l]].Name);
    fflush(OutputFilePtr[CurrentSystem]);

    sprintf(buffer,"ASCI_Grids/asci_grid_%s.grid",PseudoAtoms[GridTypeList[l]].Name);
    FilePtr=fopen(buffer,"w");

    typeA=GridTypeList[l];
    for(i=0;i<=NumberOfVDWGridPoints.x;i++)
    {
      teller=teller+1.0;
      for(j=0;j<=NumberOfVDWGridPoints.y;j++)
      {
        for(k=0;k<=NumberOfVDWGridPoints.z;k++)
        {
          switch(BoundaryCondition[CurrentSystem])
          {
            case RECTANGULAR:
            case TRICLINIC:
            default:
              pos.x=i*SizeGrid.x/NumberOfVDWGridPoints.x+ShiftGrid.x;
              pos.y=j*SizeGrid.y/NumberOfVDWGridPoints.y+ShiftGrid.y;
              pos.z=k*SizeGrid.z/NumberOfVDWGridPoints.z+ShiftGrid.z;

              // apply boundary condition
              //value=CalculateFrameworkVDWEnergyAtPosition(pos,typeA);
              CalculateDerivativesAtPositionVDW(pos,typeA,&value,&first_derivative,&second_derivative,&third_derivative);
              break;
          }

          // cap the value
          if(value<EnergyOverlapCriteria) 
            fprintf(FilePtr,"%g %g %g %g %g %g %g\n",pos.x,pos.y,pos.z,value*ENERGY_TO_KELVIN,
                    -first_derivative.x*ENERGY_TO_KELVIN,-first_derivative.y*ENERGY_TO_KELVIN,-first_derivative.z*ENERGY_TO_KELVIN);
          else
            fprintf(FilePtr,"%g %g %g %s %s %s %s\n",pos.x,pos.y,pos.z,"?","?","?","?");
        }
      }
      fprintf(OutputFilePtr[CurrentSystem],"Percentage finished                      : %d\n",(int)(teller*percent));
      fflush(OutputFilePtr[CurrentSystem]);
    }
    fclose(FilePtr);
  }
}

void MakeGrid(void)
{
	int i,j,k,l,m,tel2,typeA,save;
	int Inaccessible=0;
	REAL percent,teller,TailEnergy;
	POINT pos;
	REAL value,third_derivative;
	VECTOR first_derivative;
	REAL_MATRIX3x3 second_derivative;
	char buffer3[256];
	FILE *FilePtr2;
	FILE *FilePtr3;
	FILE *FilePtr4;
	
	
	if(MIN3(BoxProperties[CurrentSystem].cx,BoxProperties[CurrentSystem].cy,
			BoxProperties[CurrentSystem].cz)<CutOffVDW)
	{
		printf("ERROR:  Cutoff smaller than half of one of the perpendicular boxlengths !!!\n");
		printf("        (Cutoff: %lf perpendicular boxlengths: %lf %lf %lf)\n",(double)CutOffVDW,
			   (double)BoxProperties[CurrentSystem].cx,(double)BoxProperties[CurrentSystem].cy,(double)BoxProperties[CurrentSystem].cz);
		printf("Advice: choose more unitcells for constructing the grid.\n");
		exit(0);
	}
	
	NumberOfAdsorbateMolecules[CurrentSystem]=0;
	TotalAdsorbateMoleculesAdded[CurrentSystem]=0;
	NumberOfCationMolecules[CurrentSystem]=0;
	
	CalculateForce();
	TailEnergy=UTailCorrection[CurrentSystem];
	
	save=ChargeMethod;
	ChargeMethod=NONE;
	
	// compute the size of the grid and the shift of the origin of the enclosing box
	SizeGrid.x=SizeGrid.y=SizeGrid.z=0.0;
	ShiftGrid.x=ShiftGrid.y=ShiftGrid.z=0.0;
	
	SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].ax);
	SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].ay);
	SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].az);
	if(UnitCellBox[CurrentSystem].ax<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].ax;
	if(UnitCellBox[CurrentSystem].ay<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].ay;
	if(UnitCellBox[CurrentSystem].az<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].az;
	
	SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].bx);
	SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].by);
	SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].bz);
	if(UnitCellBox[CurrentSystem].bx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].bx;
	if(UnitCellBox[CurrentSystem].by<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].by;
	if(UnitCellBox[CurrentSystem].bz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].bz;
	
	SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].cx);
	SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].cy);
	SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].cz);
	if(UnitCellBox[CurrentSystem].cx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].cx;
	if(UnitCellBox[CurrentSystem].cy<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].cy;
	if(UnitCellBox[CurrentSystem].cz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].cz;
	
	// compute the number of grid points
	NumberOfVDWGridPoints.x=(int)(SizeGrid.x/SpacingVDWGrid);
	NumberOfVDWGridPoints.y=(int)(SizeGrid.y/SpacingVDWGrid);
	NumberOfVDWGridPoints.z=(int)(SizeGrid.z/SpacingVDWGrid);
	
	if(NumberOfVDWGridPoints.x%2==0) NumberOfVDWGridPoints.x++;
	if(NumberOfVDWGridPoints.y%2==0) NumberOfVDWGridPoints.y++;
	if(NumberOfVDWGridPoints.z%2==0) NumberOfVDWGridPoints.z++;
	
	DeltaVDWGrid.x=SizeGrid.x/(REAL)NumberOfVDWGridPoints.x;
	DeltaVDWGrid.y=SizeGrid.y/(REAL)NumberOfVDWGridPoints.y;
	DeltaVDWGrid.z=SizeGrid.z/(REAL)NumberOfVDWGridPoints.z;
	
	printf("ShiftGrid: %g %g %g\n",ShiftGrid.x,ShiftGrid.y,ShiftGrid.z);
	printf("SizeGrid: %g %g %g\n",SizeGrid.x,SizeGrid.y,SizeGrid.z);
	printf("Number of grid points: %d %d %d\n",NumberOfVDWGridPoints.x,NumberOfVDWGridPoints.y,NumberOfVDWGridPoints.z);
	
	percent=100.0/(REAL)((NumberOfVDWGridPoints.x+1)*NumberOfGrids);
	teller=0.0;
	tel2=1;
	
	fprintf(OutputFilePtr[CurrentSystem],"\n\n");
	fprintf(OutputFilePtr[CurrentSystem],"Generating an VDW interpolation grid (%d x %d x %d)\n",NumberOfVDWGridPoints.x,NumberOfVDWGridPoints.y,NumberOfVDWGridPoints.z);
	fprintf(OutputFilePtr[CurrentSystem],"========================================================\n\n");
	fflush(OutputFilePtr[CurrentSystem]);
	
	printf("NumberOfGrids: %d\n",NumberOfGrids);
	for(l=0;l<NumberOfGrids;l++)
	{
		fprintf(OutputFilePtr[CurrentSystem],"Creating grid %d [%s]\n",l,PseudoAtoms[GridTypeList[l]].Name);
		fflush(OutputFilePtr[CurrentSystem]);
		
		// Allocate the memory for this grid
		VDWGrid[GridTypeList[l]]=(float****)calloc(NumberOfVDWGridPoints.x+1,sizeof(float***));
		for(i=0;i<=NumberOfVDWGridPoints.x;i++)
		{
			VDWGrid[GridTypeList[l]][i]=(float***)calloc(NumberOfVDWGridPoints.y+1,sizeof(float**));
			for(j=0;j<=NumberOfVDWGridPoints.y;j++)
			{
				VDWGrid[GridTypeList[l]][i][j]=(float**)calloc(NumberOfVDWGridPoints.z+1,sizeof(float*));
				for(k=0;k<=NumberOfVDWGridPoints.z;k++)
					VDWGrid[GridTypeList[l]][i][j][k]=(float*)calloc(8,sizeof(float));
			}
		}
		
		typeA=GridTypeList[l];
		for(i=0;i<=NumberOfVDWGridPoints.x;i++)
		{
			teller=teller+1.0;
			for(j=0;j<=NumberOfVDWGridPoints.y;j++)
			{
				for(k=0;k<=NumberOfVDWGridPoints.z;k++)
				{
					switch(BoundaryCondition[CurrentSystem])
					{
						case RECTANGULAR:
						case TRICLINIC:
						default:
							pos.x=i*SizeGrid.x/NumberOfVDWGridPoints.x+ShiftGrid.x;
							pos.y=j*SizeGrid.y/NumberOfVDWGridPoints.y+ShiftGrid.y;
							pos.z=k*SizeGrid.z/NumberOfVDWGridPoints.z+ShiftGrid.z;
							
							// apply boundary condition
							CalculateDerivativesAtPositionVDW(pos,typeA,&value,&first_derivative,&second_derivative,&third_derivative);
							break;
					}
					
					// cap the value
					if(value>1.0e12) 
					{
						value=1.0e12;
						if(first_derivative.x>EnergyOverlapCriteria) first_derivative.x=EnergyOverlapCriteria;
						if(first_derivative.x<-EnergyOverlapCriteria) first_derivative.x=-EnergyOverlapCriteria;
						if(first_derivative.y>EnergyOverlapCriteria) first_derivative.y=EnergyOverlapCriteria;
						if(first_derivative.y<-EnergyOverlapCriteria) first_derivative.y=-EnergyOverlapCriteria;
						if(first_derivative.z>EnergyOverlapCriteria) first_derivative.z=EnergyOverlapCriteria;
						if(first_derivative.z<-EnergyOverlapCriteria) first_derivative.z=-EnergyOverlapCriteria;
						second_derivative.ax=second_derivative.bx=second_derivative.cx=0.0;
						second_derivative.ay=second_derivative.by=second_derivative.cy=0.0;
						second_derivative.az=second_derivative.bz=second_derivative.cz=0.0;
						third_derivative=0.0;
					}
					
					// store energy, and derivatives in "unit cube format"
					VDWGrid[GridTypeList[l]][i][j][k][0]=(float)value;
					VDWGrid[GridTypeList[l]][i][j][k][1]=(float)first_derivative.x*DeltaVDWGrid.x;
					VDWGrid[GridTypeList[l]][i][j][k][2]=(float)first_derivative.y*DeltaVDWGrid.y;
					VDWGrid[GridTypeList[l]][i][j][k][3]=(float)first_derivative.z*DeltaVDWGrid.z;
					VDWGrid[GridTypeList[l]][i][j][k][4]=(float)second_derivative.ay*(DeltaVDWGrid.x*DeltaVDWGrid.y);
					VDWGrid[GridTypeList[l]][i][j][k][5]=(float)second_derivative.az*(DeltaVDWGrid.x*DeltaVDWGrid.z);
					VDWGrid[GridTypeList[l]][i][j][k][6]=(float)second_derivative.bz*(DeltaVDWGrid.y*DeltaVDWGrid.z);
					VDWGrid[GridTypeList[l]][i][j][k][7]=(float)third_derivative*(DeltaVDWGrid.x*DeltaVDWGrid.y*DeltaVDWGrid.z);
				}
			}
			if((int)(teller*percent)>(int)((teller-1.0)*percent))
				fprintf(OutputFilePtr[CurrentSystem],"Percentage finished                      : %d\n",(int)(teller*percent));
			fflush(OutputFilePtr[CurrentSystem]);
		}
		printf("Writing Grid\n");
		WriteVDWGrid(GridTypeList[l]);
	}
	
	//Read in Delaunay polyhedra and find planes for faces
	if (DelaunayPolyhedraFile)
	{
		
		ReadPolyhedraData();
		
		//Check to see if RASPA unit cell is same size as Delaunay unit cell
		/* Code here */
		
		for (i=0; i<NumberPolyhedra; i++) {
			for (j=0; j<PolyNumberVertices[i]; j++) {
				//Set PolyVeertLocation to the value at the atom inside the box and shifted by PBCs from PolyDels
				PolyVertLocation[i][j] = PBCCorrection(Framework[0].Atoms[0][PolyVertices[i][j]-1].Position,i,j);
			}
			for (j=0; j<PolyNumberFaces[i]; j++) {
				
				PlaneFinder(i,j);
				for (k=0; k<PolyNumVertPerFace[i][j]; k++) {
					
					//Place all vertices in each face into 2D, stored in Poly2DPoints
					Poly2DPoints[i][j][k]=CoordinateSystemShift(PolyVertLocation[i][PolyVertPerFace[i][j][k]-1], i, j);
				}
			}
		}
	}
	// Creation of the Delaunay Grid points
	// compute the number of grid points
	NumberOfDelaunayGridPoints.x=(int)(SizeGrid.x/SpacingDelaunayGrid);
	NumberOfDelaunayGridPoints.y=(int)(SizeGrid.y/SpacingDelaunayGrid);
	NumberOfDelaunayGridPoints.z=(int)(SizeGrid.z/SpacingDelaunayGrid);
	
	if(NumberOfDelaunayGridPoints.x%2==0) NumberOfDelaunayGridPoints.x++;
	if(NumberOfDelaunayGridPoints.y%2==0) NumberOfDelaunayGridPoints.y++;
	if(NumberOfDelaunayGridPoints.z%2==0) NumberOfDelaunayGridPoints.z++;
	
	DeltaDelaunayGrid.x=SizeGrid.x/(REAL)NumberOfDelaunayGridPoints.x;
	DeltaDelaunayGrid.y=SizeGrid.y/(REAL)NumberOfDelaunayGridPoints.y;
	DeltaDelaunayGrid.z=SizeGrid.z/(REAL)NumberOfDelaunayGridPoints.z;
	
	printf("ShiftGrid: %g %g %g\n",ShiftGrid.x,ShiftGrid.y,ShiftGrid.z);
	printf("SizeGrid: %g %g %g\n",SizeGrid.x,SizeGrid.y,SizeGrid.z);
	printf("Number of grid points: %d %d %d\n",NumberOfDelaunayGridPoints.x,NumberOfDelaunayGridPoints.y,NumberOfDelaunayGridPoints.z);
	
	percent=100.0/(REAL)((NumberOfDelaunayGridPoints.x+1)*NumberOfGrids);
	teller=0.0;
	tel2=1;
	
	fprintf(OutputFilePtr[CurrentSystem],"\n\n");
	fprintf(OutputFilePtr[CurrentSystem],"Generating a Delaunay accessible grid (%d x %d x %d)\n",NumberOfDelaunayGridPoints.x,NumberOfDelaunayGridPoints.y,NumberOfDelaunayGridPoints.z);
	fprintf(OutputFilePtr[CurrentSystem],"========================================================\n\n");
	fflush(OutputFilePtr[CurrentSystem]);
	
	printf("NumberOfGrids: %d\n",NumberOfGrids);
	for(l=0;l<NumberOfGrids;l++)
	{
		fprintf(OutputFilePtr[CurrentSystem],"Creating delaunay grid %d [%s]\n",l,PseudoAtoms[GridTypeList[l]].Name);
		fflush(OutputFilePtr[CurrentSystem]);
		
		// Allocate the memory for this grid
		DelaunayGrid[GridTypeList[l]]=(int****)calloc(NumberOfDelaunayGridPoints.x+1,sizeof(int***));
		for(i=0;i<=NumberOfDelaunayGridPoints.x;i++)
		{
			DelaunayGrid[GridTypeList[l]][i]=(int***)calloc(NumberOfDelaunayGridPoints.y+1,sizeof(int**));
			for(j=0;j<=NumberOfDelaunayGridPoints.y;j++)
			{
				DelaunayGrid[GridTypeList[l]][i][j]=(int**)calloc(NumberOfDelaunayGridPoints.z+1,sizeof(int*));
				for(k=0;k<=NumberOfDelaunayGridPoints.z;k++)
					DelaunayGrid[GridTypeList[l]][i][j][k]=(int*)calloc(1,sizeof(int));
			}
		}
		
		sprintf(buffer3,"FrameworkGridForVMD.xyz");
		fopen(buffer3,"w");
		FilePtr2 = fopen(buffer3,"r+");
		
		sprintf(buffer3,"AccessGridForVMD.xyz");
		fopen(buffer3,"w");
		FilePtr3 = fopen(buffer3,"r+");
		
		sprintf(buffer3,"InAccessGridForVMD.xyz");
		fopen(buffer3,"w");
		FilePtr4 = fopen(buffer3,"r+");
		
		typeA=GridTypeList[l];
		
		for(i=0;i<=NumberOfDelaunayGridPoints.x;i++)
		{
			teller=teller+1.0;
			for(j=0;j<=NumberOfDelaunayGridPoints.y;j++)
			{
				for(k=0;k<=NumberOfDelaunayGridPoints.z;k++)
				{
					switch(BoundaryCondition[CurrentSystem])
					{
						case RECTANGULAR:
						case TRICLINIC:
						default:
							pos.x=i*SizeGrid.x/NumberOfDelaunayGridPoints.x+ShiftGrid.x;
							pos.y=j*SizeGrid.y/NumberOfDelaunayGridPoints.y+ShiftGrid.y;
							pos.z=k*SizeGrid.z/NumberOfDelaunayGridPoints.z+ShiftGrid.z;
							
							// apply boundary condition
							
							break;
					}
					if (SpacingDelaunayGrid == SpacingVDWGrid && VDWGrid[GridTypeList[l]][i][j][k][0] >= EnergyOverlapCriteria) {
						Inaccessible = 1;
						fprintf(FilePtr2,"g %f %f %f\n",pos.x,pos.y,pos.z);
						fflush(FilePtr2);
					}
					
					else if (DelaunayPolyhedraFile && DelaunayGrid[GridTypeList[l]][i][j][k][0] < 1 && NumberPolyhedra > 0){
						Inaccessible = 0;
						m = 0;
						while (m<NumberPolyhedra && Inaccessible == 0) {
							Inaccessible=PointInPolyhedra(m, Inaccessible, pos);								
							if (Inaccessible == 1) {
								fprintf(FilePtr4,"g %f %f %f\n",pos.x,pos.y,pos.z);
								fflush(FilePtr4);
							}
							m++;
						}
						
						if (Inaccessible == 0) {
							fprintf(FilePtr3,"g %f %f %f\n",pos.x,pos.y,pos.z);
							fflush(FilePtr3);
						}
					}
					else {
						Inaccessible = 0;
					}
					
					// store accessibility in "unit cube format"
					DelaunayGrid[GridTypeList[l]][i][j][k][0]=Inaccessible;
				}
			}
			if((int)(teller*percent)>(int)((teller-1.0)*percent))
				fprintf(OutputFilePtr[CurrentSystem],"Percentage finished                      : %d\n",(int)(teller*percent));
			fflush(OutputFilePtr[CurrentSystem]);
		}
		printf("Writing Grid\n");
		WriteDelaunayGrid(GridTypeList[l]);
	}
	//Free the memory intensive variables
	free(PolyNumberFaces);
	free(PolyNumberVertices);
	free(PolyNumVertPerFace);
	free(PolyVertices);
	free(PolyVertPerFace);
	free(PolyNumberEdges);
	free(PolyEdges);
	free(PolyPlanes);
	free(projectedpoint);
	free(Poly2DPoints);
	free(Plane2D);
	free(PolyVertLocation);
	free(PolyDelsVertices);
	
	// compute the number of grid points
	NumberOfCoulombGridPoints.x=(int)(SizeGrid.x/SpacingCoulombGrid);
	NumberOfCoulombGridPoints.y=(int)(SizeGrid.y/SpacingCoulombGrid);
	NumberOfCoulombGridPoints.z=(int)(SizeGrid.z/SpacingCoulombGrid);
	
	if(NumberOfCoulombGridPoints.x%2==0) NumberOfCoulombGridPoints.x++;
	if(NumberOfCoulombGridPoints.y%2==0) NumberOfCoulombGridPoints.y++;
	if(NumberOfCoulombGridPoints.z%2==0) NumberOfCoulombGridPoints.z++;
	
	DeltaCoulombGrid.x=SizeGrid.x/(REAL)NumberOfCoulombGridPoints.x;
	DeltaCoulombGrid.y=SizeGrid.y/(REAL)NumberOfCoulombGridPoints.y;
	DeltaCoulombGrid.z=SizeGrid.z/(REAL)NumberOfCoulombGridPoints.z;
	
	ChargeMethod=save;
	fprintf(OutputFilePtr[CurrentSystem],"\n\n");
	switch(ChargeMethod)
	{
		case EWALD:
			fprintf(OutputFilePtr[CurrentSystem],"Generating an Ewald-Coulomb interpolation grid (%d x %d x %d)\n",
					NumberOfCoulombGridPoints.x,NumberOfCoulombGridPoints.y,NumberOfCoulombGridPoints.z);
			break;
		case WOLFS_METHOD:
			fprintf(OutputFilePtr[CurrentSystem],"Generating an Coulomb interpolation grid using Wolf's method (%d x %d x %d)\n",
					NumberOfCoulombGridPoints.x,NumberOfCoulombGridPoints.y,NumberOfCoulombGridPoints.z);
			break;
		default:
			printf("Skipping the Coulombic grid\n");
			exit(0);
			break;
	}
	fprintf(OutputFilePtr[CurrentSystem],"===============================================================\n\n");
	fflush(OutputFilePtr[CurrentSystem]);
	
	percent=100.0/(REAL)(NumberOfCoulombGridPoints.x+1);
	teller=0.0;
	
	// Allocate the memory for this grid
	CoulombGrid=(float****)calloc(NumberOfCoulombGridPoints.x+1,sizeof(float***));
	for(i=0;i<=NumberOfCoulombGridPoints.x;i++)
	{
		CoulombGrid[i]=(float***)calloc(NumberOfCoulombGridPoints.y+1,sizeof(float**));
		for(j=0;j<=NumberOfCoulombGridPoints.y;j++)
		{
			CoulombGrid[i][j]=(float**)calloc(NumberOfCoulombGridPoints.z+1,sizeof(float*));
			for(k=0;k<=NumberOfCoulombGridPoints.z;k++)
				CoulombGrid[i][j][k]=(float*)calloc(8,sizeof(float));
		}
	}
	
	typeA=ReturnPseudoAtomNumber("UNIT");
	for(i=0;i<=NumberOfCoulombGridPoints.x;i++)
	{
		teller=teller+1.0;
		for(j=0;j<=NumberOfCoulombGridPoints.y;j++)
		{
			for(k=0;k<=NumberOfCoulombGridPoints.z;k++)
			{
				switch(BoundaryCondition[CurrentSystem])
				{
					case RECTANGULAR:
					case TRICLINIC:
					default:
						pos.x=i*SizeGrid.x/NumberOfCoulombGridPoints.x+ShiftGrid.x;
						pos.y=j*SizeGrid.y/NumberOfCoulombGridPoints.y+ShiftGrid.y;
						pos.z=k*SizeGrid.z/NumberOfCoulombGridPoints.z+ShiftGrid.z;
						
						CalculateDerivativesAtPositionReal(pos,typeA,&value,&first_derivative,&second_derivative,&third_derivative);
						break;
				}
				
				// cap the value
				if(value>1.0e12)
				{
					value=1.0e12;
					if(first_derivative.x>EnergyOverlapCriteria) first_derivative.x=EnergyOverlapCriteria;
					if(first_derivative.x<-EnergyOverlapCriteria) first_derivative.x=-EnergyOverlapCriteria;
					if(first_derivative.y>EnergyOverlapCriteria) first_derivative.y=EnergyOverlapCriteria;
					if(first_derivative.y<-EnergyOverlapCriteria) first_derivative.y=-EnergyOverlapCriteria;
					if(first_derivative.z>EnergyOverlapCriteria) first_derivative.z=EnergyOverlapCriteria;
					if(first_derivative.z<-EnergyOverlapCriteria) first_derivative.z=-EnergyOverlapCriteria;
					second_derivative.ax=second_derivative.bx=second_derivative.cx=0.0;
					second_derivative.ay=second_derivative.by=second_derivative.cy=0.0;
					second_derivative.az=second_derivative.bz=second_derivative.cz=0.0;
					third_derivative=0.0;
				}
				
				// store energy, and derivatives in "unit cube format"
				CoulombGrid[i][j][k][0]=(float)value;
				CoulombGrid[i][j][k][1]=(float)first_derivative.x*DeltaCoulombGrid.x;
				CoulombGrid[i][j][k][2]=(float)first_derivative.y*DeltaCoulombGrid.y;
				CoulombGrid[i][j][k][3]=(float)first_derivative.z*DeltaCoulombGrid.z;
				CoulombGrid[i][j][k][4]=(float)second_derivative.ay*(DeltaCoulombGrid.x*DeltaCoulombGrid.y);
				CoulombGrid[i][j][k][5]=(float)second_derivative.az*(DeltaCoulombGrid.x*DeltaCoulombGrid.z);
				CoulombGrid[i][j][k][6]=(float)second_derivative.bz*(DeltaCoulombGrid.y*DeltaCoulombGrid.z);
				CoulombGrid[i][j][k][7]=(float)third_derivative*(DeltaCoulombGrid.x*DeltaCoulombGrid.y*DeltaCoulombGrid.z);
			}
		}
		
		if((int)(teller*percent)>(int)((teller-1.0)*percent))
			fprintf(OutputFilePtr[CurrentSystem],"Percentage finished                      : %d\n",(int)(teller*percent));
		fflush(OutputFilePtr[CurrentSystem]);
	}
	printf("Writing Grid\n");
	WriteCoulombGrid();
}

void MakeDynamicGrid(void)
{
	int i,j,k;
	
	//	//Find smallest radius used and use to make dynamic grid (Still need to figure out if want to use this way)
	//	//Currently using user defined grid size	
	//	Components[comp].NumberOfAtoms
	//	Components[comp].Type[temp]
	//	for (i=0; i<NumberOfComponents; i++) {
	//		for (j=0; j<Components[i].NumberOfAtoms; j++) {
	//			PseudoAtoms[Components[i].Type[j]].Radius
	//			
	//		}
	//	}
	
	TotalAdsorbateMoleculesAdded[CurrentSystem]=0;
	
	// compute the size of the grid and the shift of the origin of the enclosing box
	SizeGrid.x=SizeGrid.y=SizeGrid.z=0.0;
	ShiftGrid.x=ShiftGrid.y=ShiftGrid.z=0.0;
	
	SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].ax);
	SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].ay);
	SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].az);
	if(UnitCellBox[CurrentSystem].ax<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].ax;
	if(UnitCellBox[CurrentSystem].ay<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].ay;
	if(UnitCellBox[CurrentSystem].az<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].az;
	
	SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].bx);
	SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].by);
	SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].bz);
	if(UnitCellBox[CurrentSystem].bx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].bx;
	if(UnitCellBox[CurrentSystem].by<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].by;
	if(UnitCellBox[CurrentSystem].bz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].bz;
	
	SizeGrid.x+=fabs(UnitCellBox[CurrentSystem].cx);
	SizeGrid.y+=fabs(UnitCellBox[CurrentSystem].cy);
	SizeGrid.z+=fabs(UnitCellBox[CurrentSystem].cz);
	if(UnitCellBox[CurrentSystem].cx<0.0) ShiftGrid.x+=UnitCellBox[CurrentSystem].cx;
	if(UnitCellBox[CurrentSystem].cy<0.0) ShiftGrid.y+=UnitCellBox[CurrentSystem].cy;
	if(UnitCellBox[CurrentSystem].cz<0.0) ShiftGrid.z+=UnitCellBox[CurrentSystem].cz;
	
	// compute the number of grid points
	NumberOfDynamicGridPoints.x=(int)(SizeGrid.x/SpacingDynamicGrid);
	NumberOfDynamicGridPoints.y=(int)(SizeGrid.y/SpacingDynamicGrid);
	NumberOfDynamicGridPoints.z=(int)(SizeGrid.z/SpacingDynamicGrid);
	
	if(NumberOfDynamicGridPoints.x%2==0) NumberOfVDWGridPoints.x++;
	if(NumberOfDynamicGridPoints.y%2==0) NumberOfVDWGridPoints.y++;
	if(NumberOfDynamicGridPoints.z%2==0) NumberOfVDWGridPoints.z++;
	
	DeltaDynamicGrid.x=SizeGrid.x/(REAL)NumberOfDynamicGridPoints.x;
	DeltaDynamicGrid.y=SizeGrid.y/(REAL)NumberOfDynamicGridPoints.y;
	DeltaDynamicGrid.z=SizeGrid.z/(REAL)NumberOfDynamicGridPoints.z;
	
	printf("ShiftGrid: %g %g %g\n",ShiftGrid.x,ShiftGrid.y,ShiftGrid.z);
	printf("SizeGrid: %g %g %g\n",SizeGrid.x,SizeGrid.y,SizeGrid.z);
	printf("Number of grid points: %d %d %d\n",NumberOfDynamicGridPoints.x,NumberOfDynamicGridPoints.y,NumberOfDynamicGridPoints.z);
	
	// Allocate the memory for this dynamic grid
	DynamicGrid=(int****)calloc(NumberOfDynamicGridPoints.x+1,sizeof(int***));
	for(i=0;i<=NumberOfDynamicGridPoints.x;i++)
	{
		DynamicGrid[i]=(int***)calloc(NumberOfDynamicGridPoints.y+1,sizeof(int**));
		for(j=0;j<=NumberOfDynamicGridPoints.y;j++)
		{
			DynamicGrid[i][j]=(int**)calloc(NumberOfDynamicGridPoints.z+1,sizeof(int*));
			for(k=0;k<=NumberOfDynamicGridPoints.z;k++)
				DynamicGrid[i][j][k]=(int*)calloc(1,sizeof(int));
		}
	}
	
}

void AdjustDynamicGrid(int GridValue, int typeA, VECTOR pos)
{
	int x0,y0,z0,x1,y1,z1;
	VECTOR s;
	
	switch(BoundaryCondition[CurrentSystem])
	{
		case RECTANGULAR:
			// the position has to be moved back to the main unit cell using the rectangular boundary condition
			pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
			pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
			pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
			if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
			if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
			if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
			break;
		case TRICLINIC:
		default:
			// the position first has to be moved back to the main unit cell using the triclinic boundary condition
			s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
			s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
			s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;
			
			s.x-=(REAL)NINT(s.x);
			s.y-=(REAL)NINT(s.y);
			s.z-=(REAL)NINT(s.z);
			
			if(s.x<0.0) s.x+=1.0;
			if(s.y<0.0) s.y+=1.0;
			if(s.z<0.0) s.z+=1.0;
			
			pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
			pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
			pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
			break;
	}
	
	
	s.x=((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)(NumberOfDynamicGridPoints.x);
	s.y=((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)(NumberOfDynamicGridPoints.y);
	s.z=((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)(NumberOfDynamicGridPoints.z);
	
	// find the corresponding cell, between 0 and NumberOfVDWGridPoints
	// determine lower boundary
	x0=(int)s.x;
	y0=(int)s.y;
	z0=(int)s.z;
	
	// determine upper boundary, apply periodic boundary condition
	x1=x0+1;
	y1=y0+1;
	z1=z0+1;
	
	if (x1 == NumberOfDynamicGridPoints.x) 
		x1=0;
	if (y1 == NumberOfDynamicGridPoints.y) 
		y1=0;
	if (z1 == NumberOfDynamicGridPoints.z) 
		z1=0;
	
	// Change Grid Value
	
	DynamicGrid[x0][y0][z0][0]=GridValue;		
	DynamicGrid[x1][y0][z0][0]=GridValue;
	DynamicGrid[x0][y1][z0][0]=GridValue;
	DynamicGrid[x0][y0][z1][0]=GridValue;
	DynamicGrid[x1][y1][z0][0]=GridValue;	
	DynamicGrid[x1][y0][z1][0]=GridValue;	
	DynamicGrid[x0][y1][z1][0]=GridValue;
	DynamicGrid[x1][y1][z1][0]=GridValue;
	
	
}

void PrintDynamicGrid(void)
{
	FILE *FilePtr2;
	char buffer3[256];
	int i,j,k,ii,Total,AdsorbNumber;
	VECTOR pos;
	
	Total=TotalAdsorbateMoleculesAdded[CurrentSystem];
	sprintf(buffer3,"DynamicGrid.xyz");
	fopen(buffer3,"a+");
	FilePtr2 = fopen(buffer3,"a+");
	fprintf(FilePtr2,"New Grid\n");
	for (ii=0; ii<=Total; ii++) {
		for(i=0;i<=NumberOfDynamicGridPoints.x;i++)
		{
			for(j=0;j<=NumberOfDynamicGridPoints.y;j++)
			{
				for(k=0;k<=NumberOfDynamicGridPoints.z;k++)
				{
					switch(BoundaryCondition[CurrentSystem])
					{
						case RECTANGULAR:
						case TRICLINIC:
						default:
							pos.x=i*SizeGrid.x/NumberOfDynamicGridPoints.x+ShiftGrid.x;
							pos.y=j*SizeGrid.y/NumberOfDynamicGridPoints.y+ShiftGrid.y;
							pos.z=k*SizeGrid.z/NumberOfDynamicGridPoints.z+ShiftGrid.z;
							
							break;
					}
					AdsorbNumber=DynamicGrid[i][j][k][0];
					if (AdsorbNumber==ii) {
						fprintf(FilePtr2,"%i %f %f %f\n",DynamicGrid[i][j][k][0], pos.x,pos.y,pos.z);
					}
				}
			}
		}
	}
	fflush(FilePtr2);
}

int WriteVDWGrid(int l)
{
  int i,j,k,m;
  int ngrid;
  FILE *FilePtr;
  char buffer[256];

  ngrid=NumberOfVDWGridPoints.x*NumberOfVDWGridPoints.y*NumberOfVDWGridPoints.z;

  sprintf(buffer,"%s/share/raspa/grids",RASPA_DIRECTORY);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s",RASPA_DIRECTORY,ForceField);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s",RASPA_DIRECTORY,ForceField,Framework[CurrentSystem].Name[0]);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf",RASPA_DIRECTORY,ForceField,
    Framework[CurrentSystem].Name[0],(double)SpacingVDWGrid);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%s_%s_%s.grid",
      RASPA_DIRECTORY,
      ForceField,
      Framework[CurrentSystem].Name[0],
      (double)SpacingVDWGrid,
      Framework[CurrentSystem].Name[0],
      PseudoAtoms[l].Name,
      ShiftPotentials?"shifted":"truncated");
  FilePtr=fopen(buffer,"w");

  fwrite(&SpacingVDWGrid,1,sizeof(REAL),FilePtr);
  fwrite(&NumberOfVDWGridPoints,1,sizeof(INT_VECTOR3),FilePtr);
  fwrite(&SizeGrid,1,sizeof(VECTOR),FilePtr);
  fwrite(&ShiftGrid,1,sizeof(VECTOR),FilePtr);
  fwrite(&DeltaVDWGrid,1,sizeof(VECTOR),FilePtr);
  fwrite(&UnitCellSize[CurrentSystem],1,sizeof(VECTOR),FilePtr);
  fwrite(&NumberOfUnitCells[CurrentSystem],1,sizeof(INT_VECTOR3),FilePtr);

  for(m=0;m<8;m++)
    for(i=0;i<=NumberOfVDWGridPoints.x;i++)
      for(j=0;j<=NumberOfVDWGridPoints.y;j++)
        for(k=0;k<=NumberOfVDWGridPoints.z;k++)
          fwrite(&VDWGrid[l][i][j][k][m],1,sizeof(float),FilePtr);


  fclose(FilePtr);
  return 0;
}

void ReadVDWGrid(void)
{
  int i,j,k,l,m;
  VECTOR unit_cell_size;
  INT_VECTOR3 number_of_unit_cells;
  FILE *FilePtr;
  char buffer[256];

  printf("Reading VDW grid\n");

  for(i=0;i<NumberOfPseudoAtoms;i++)
    VDWGrid[i]=NULL;

  for(l=0;l<NumberOfGrids;l++)
  {
    sprintf(buffer,"%s/share/raspa/grids",RASPA_DIRECTORY);
    sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%s_%s_%s.grid",
        RASPA_DIRECTORY,
        ForceField,
        Framework[CurrentSystem].Name[0],
        (double)SpacingVDWGrid,
        Framework[CurrentSystem].Name[0],
        PseudoAtoms[GridTypeList[l]].Name,
        ShiftPotentials?"shifted":"truncated");
    printf("Opening: %s\n",buffer);
    if(!(FilePtr=fopen(buffer,"r")))
    {
      printf("Error:  file %s does not exists.\n",buffer);
      exit(1);
    }

    fread(&SpacingVDWGrid,1,sizeof(REAL),FilePtr);
    fread(&NumberOfVDWGridPoints,1,sizeof(INT_VECTOR3),FilePtr);
    fread(&SizeGrid,1,sizeof(VECTOR),FilePtr);
    fread(&ShiftGrid,1,sizeof(VECTOR),FilePtr);
    fread(&DeltaVDWGrid,1,sizeof(VECTOR),FilePtr);
    fread(&unit_cell_size,1,sizeof(VECTOR),FilePtr);
    fread(&number_of_unit_cells,1,sizeof(INT_VECTOR3),FilePtr);

    if((fabs(UnitCellSize[CurrentSystem].x-unit_cell_size.x)>1e-8)||(fabs(UnitCellSize[CurrentSystem].y-unit_cell_size.y)>1e-8)||
       (fabs(UnitCellSize[CurrentSystem].z-unit_cell_size.z)>1e-8))
    {
      printf("Error in reading the grid: unit-cell size does not match\n");
      exit(0);
    }

    // Allocate the memory for this grid
    VDWGrid[GridTypeList[l]]=(float****)calloc(NumberOfVDWGridPoints.x+1,sizeof(float***));
    for(i=0;i<=NumberOfVDWGridPoints.x;i++)
    {
      VDWGrid[GridTypeList[l]][i]=(float***)calloc(NumberOfVDWGridPoints.y+1,sizeof(float**));
      for(j=0;j<=NumberOfVDWGridPoints.y;j++)
      {
        VDWGrid[GridTypeList[l]][i][j]=(float**)calloc(NumberOfVDWGridPoints.z+1,sizeof(float*));
        for(k=0;k<=NumberOfVDWGridPoints.z;k++)
          VDWGrid[GridTypeList[l]][i][j][k]=(float*)calloc(8,sizeof(float));
      }
    }

    for(m=0;m<8;m++)
      for(i=0;i<=NumberOfVDWGridPoints.x;i++)
        for(j=0;j<=NumberOfVDWGridPoints.y;j++)
          for(k=0;k<=NumberOfVDWGridPoints.z;k++)
            fread(&VDWGrid[GridTypeList[l]][i][j][k][m],1,sizeof(float),FilePtr);


    fclose(FilePtr);
  }
}

int WriteCoulombGrid(void)
{
  int i,j,k,m,max_size;
  FILE *FilePtr;
  char buffer[256],name[256];

  switch(ChargeMethod)
  {
    case EWALD:
      strcpy(name,"Ewald");
      break;
    case WOLFS_METHOD:
      strcpy(name,"WolfsMethod");
      break;
  }

  max_size=NumberOfCoulombGridPoints.x*NumberOfCoulombGridPoints.y*NumberOfCoulombGridPoints.z;

  sprintf(buffer,"%s/share/raspa/grids",RASPA_DIRECTORY);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s",RASPA_DIRECTORY,ForceField);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s",RASPA_DIRECTORY,ForceField,Framework[CurrentSystem].Name[0]);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf",RASPA_DIRECTORY,ForceField,
    Framework[CurrentSystem].Name[0],(double)SpacingCoulombGrid);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%dx%dx%d",RASPA_DIRECTORY,ForceField,
    Framework[CurrentSystem].Name[0],(double)SpacingCoulombGrid,
    NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z);
  mkdir(buffer,S_IRWXU);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%dx%dx%d/%s_%s_%s.grid",
      RASPA_DIRECTORY,
      ForceField,
      Framework[CurrentSystem].Name[0],
      (double)SpacingCoulombGrid,
      NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z,
      Framework[CurrentSystem].Name[0],
      "Electrostatics",
      name);
  FilePtr=fopen(buffer,"w");

  fwrite(&SpacingCoulombGrid,1,sizeof(REAL),FilePtr);
  fwrite(&NumberOfCoulombGridPoints,1,sizeof(INT_VECTOR3),FilePtr);
  fwrite(&SizeGrid,1,sizeof(VECTOR),FilePtr);
  fwrite(&ShiftGrid,1,sizeof(VECTOR),FilePtr);
  fwrite(&DeltaCoulombGrid,1,sizeof(VECTOR),FilePtr);
  fwrite(&UnitCellSize[CurrentSystem],1,sizeof(VECTOR),FilePtr);
  fwrite(&NumberOfUnitCells[CurrentSystem],1,sizeof(INT_VECTOR3),FilePtr);
  fwrite(&EwaldPrecision,1,sizeof(REAL),FilePtr);

  for(m=0;m<8;m++)
    for(i=0;i<=NumberOfCoulombGridPoints.x;i++)
      for(j=0;j<=NumberOfCoulombGridPoints.y;j++)
        for(k=0;k<=NumberOfCoulombGridPoints.z;k++)
          fwrite(&CoulombGrid[i][j][k][m],1,sizeof(float),FilePtr);

  fclose(FilePtr);
  return 0;
}

void ReadCoulombGrid(void)
{
  int i,j,k,m;
  INT_VECTOR3 number_of_unit_cells;
  VECTOR unit_cell_size;
  REAL ewald_precision;
  FILE *FilePtr;
  char buffer[256],name[256];

  switch(ChargeMethod)
  {
    case EWALD:
      strcpy(name,"Ewald");
      break;
    case WOLFS_METHOD:
      strcpy(name,"WolfsMethod");
      break;
  }

  sprintf(buffer,"%s/share/raspa/grids/%s",RASPA_DIRECTORY,ForceField);
  sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%dx%dx%d/%s_%s_%s.grid",
      RASPA_DIRECTORY,
      ForceField,
      Framework[CurrentSystem].Name[0],
      (double)SpacingCoulombGrid,  
      NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z,
      Framework[CurrentSystem].Name[0],
      "Electrostatics",
      name);
  printf("Opening: %s\n",buffer);
  if(!(FilePtr=fopen(buffer,"r")))
  {
    printf("Error:  file %s does not exists.\n",buffer);
    exit(0);
  }
  else
  {
    fread(&SpacingCoulombGrid,1,sizeof(REAL),FilePtr);
    fread(&NumberOfCoulombGridPoints,1,sizeof(INT_VECTOR3),FilePtr);
    fread(&SizeGrid,1,sizeof(VECTOR),FilePtr);
    fread(&ShiftGrid,1,sizeof(VECTOR),FilePtr);
    fread(&DeltaCoulombGrid,1,sizeof(VECTOR),FilePtr);
    fread(&unit_cell_size,1,sizeof(VECTOR),FilePtr);
    fread(&number_of_unit_cells,1,sizeof(INT_VECTOR3),FilePtr);
    fread(&ewald_precision,1,sizeof(REAL),FilePtr);
  
    if((fabs(UnitCellSize[CurrentSystem].x-unit_cell_size.x)>1e-8)||(fabs(UnitCellSize[CurrentSystem].y-unit_cell_size.y)>1e-8)||
       (fabs(UnitCellSize[CurrentSystem].z-unit_cell_size.z)>1e-8))
    {
      printf("Error in reading the grid: unit-cell size does not match\n");
      exit(0);
    }
  
    if((number_of_unit_cells.x!=NumberOfUnitCells[CurrentSystem].x)||
       (number_of_unit_cells.y!=NumberOfUnitCells[CurrentSystem].y)||
       (number_of_unit_cells.z!=NumberOfUnitCells[CurrentSystem].z))
    {
      printf("Mismatch of number of unitcells, defined: %d,%d,%d, grid: %d,%d,%d\n",
          NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z,
          number_of_unit_cells.x,number_of_unit_cells.y,number_of_unit_cells.z);
      exit(0);
    }
  
    if(fabs(EwaldPrecision-ewald_precision)>1e-8)
    {
      printf("Mismatch of Ewald precision, currently: %g but the grid had been made using %g\n",
             EwaldPrecision,ewald_precision);
    }
  
    // Allocate the memory for this grid
    CoulombGrid=(float****)calloc(NumberOfCoulombGridPoints.x+1,sizeof(float***));
    for(i=0;i<=NumberOfCoulombGridPoints.x;i++)
    {
      CoulombGrid[i]=(float***)calloc(NumberOfCoulombGridPoints.y+1,sizeof(float**));
      for(j=0;j<=NumberOfCoulombGridPoints.y;j++)
      {
        CoulombGrid[i][j]=(float**)calloc(NumberOfCoulombGridPoints.z+1,sizeof(float*));
        for(k=0;k<=NumberOfCoulombGridPoints.z;k++)
          CoulombGrid[i][j][k]=(float*)calloc(8,sizeof(float));
      }
    }
  
    for(m=0;m<8;m++)
      for(i=0;i<=NumberOfCoulombGridPoints.x;i++)
        for(j=0;j<=NumberOfCoulombGridPoints.y;j++)
          for(k=0;k<=NumberOfCoulombGridPoints.z;k++)
            fread(&CoulombGrid[i][j][k][m],1,sizeof(float),FilePtr);
  
  
    fclose(FilePtr);
  }
  printf("End Reading Coulomb Grid\n");
  printf("Done !!!\n");
}

REAL InterpolateVDWGrid(int typeA,VECTOR pos)
{
	int i,j,k;
	int x0,y0,z0,x1,y1,z1;
	VECTOR s,r,H00,H01,H10,H11,h0[8],h1[8],hfinal[64];
	REAL value;
	REAL X[64],a[64];
	
	if(!VDWGrid[typeA])
	{
		printf("InterpolateVDWGrid: Pseudoatom %d [%s] not tabulated.\n",typeA,PseudoAtoms[typeA].Name);
		exit(0);
	}
	
	switch(BoundaryCondition[CurrentSystem])
	{
		case RECTANGULAR:
			// the position has to be moved back to the main unit cell using the rectangular boundary condition
			pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
			pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
			pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
			if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
			if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
			if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
			break;
		case TRICLINIC:
		default:
			// the position first has to be moved back to the main unit cell using the triclinic boundary condition
			s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
			s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
			s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;
			
			s.x-=(REAL)NINT(s.x);
			s.y-=(REAL)NINT(s.y);
			s.z-=(REAL)NINT(s.z);
			
			if(s.x<0.0) s.x+=1.0;
			if(s.y<0.0) s.y+=1.0;
			if(s.z<0.0) s.z+=1.0;
			
			pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
			pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
			pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
			break;
	}
	
	s.x=((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)(NumberOfVDWGridPoints.x);
	s.y=((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)(NumberOfVDWGridPoints.y);
	s.z=((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)(NumberOfVDWGridPoints.z);
	
	// find the corresponding cell, between 0 and NumberOfVDWGridPoints
	// determine lower boundary
	x0=(int)s.x;
	y0=(int)s.y;
	z0=(int)s.z;
	
	// determine upper boundary, apply periodic boundary condition
	x1=x0+1;
	y1=y0+1;
	z1=z0+1;
	
	// find the corresponding position within that cell (between 0.0 and 1.0)
	r.x=(s.x-x0);
	r.y=(s.y-y0);
	r.z=(s.z-z0);
	
	// retrieve the grid value and the derivatives
	for(i=0;i<8;i++)
	{
		X[0+i*8]=VDWGrid[typeA][x0][y0][z0][i];
		X[1+i*8]=VDWGrid[typeA][x1][y0][z0][i];
		X[2+i*8]=VDWGrid[typeA][x0][y1][z0][i];
		X[3+i*8]=VDWGrid[typeA][x1][y1][z0][i];
		X[4+i*8]=VDWGrid[typeA][x0][y0][z1][i];
		X[5+i*8]=VDWGrid[typeA][x1][y0][z1][i];
		X[6+i*8]=VDWGrid[typeA][x0][y1][z1][i];
		X[7+i*8]=VDWGrid[typeA][x1][y1][z1][i];
	}
	
	// cap values to avoid a negative interpolation result
	for(i=0;i<8;i++)
		if(X[i]>1e10) return 1e8;
	
	if (HermiteInterpolation) {
		H00.x = SQR(1.0-r.x)*(3.0-2.0*(1.0-r.x));
		H01.x = SQR(r.x)*(3.0-2.0*(r.x));
		H00.y = SQR(1.0-r.y)*(3.0-2.0*(1.0-r.y));
		H01.y = SQR(r.y)*(3.0-2.0*(r.y));
		H00.z = SQR(1.0-r.z)*(3.0-2.0*(1.0-r.z));
		H01.z = SQR(r.z)*(3.0-2.0*(r.z));
		
		H10.x = -DeltaVDWGrid.x*SQR(1.0-r.x)*(-r.x);
		H11.x = DeltaVDWGrid.x*SQR(r.x)*(r.x-1.0);
		H10.y = -DeltaVDWGrid.y*SQR(1.0-r.y)*(-r.y);
		H11.y = DeltaVDWGrid.y*SQR(r.y)*(r.y-1.0);
		H10.z = -DeltaVDWGrid.z*SQR(1.0-r.z)*(-r.z);
		H11.z = DeltaVDWGrid.z*SQR(r.z)*(r.z-1.0);
		
		h0[0].x = H00.x;
		h0[1].x = H01.x;
		h0[2].x = H00.x;
		h0[3].x = H00.x;
		h0[4].x = H01.x;
		h0[5].x = H01.x;
		h0[6].x = H00.x;
		h0[7].x = H01.x;
			
		h1[0].x = H10.x;
		h1[1].x = H11.x;
		h1[2].x = H10.x;
		h1[3].x = H10.x;
		h1[4].x = H11.x;
		h1[5].x = H11.x;
		h1[6].x = H10.x;
		h1[7].x = H11.x;
		
		for (i=0; i<8; i++) {
			hfinal[0+i*8].x=h0[i].x;
			hfinal[1+i*8].x=h1[i].x;
			hfinal[2+i*8].x=h0[i].x;
			hfinal[3+i*8].x=h0[i].x;
			hfinal[4+i*8].x=h1[i].x;
			hfinal[5+i*8].x=h1[i].x;
			hfinal[6+i*8].x=h0[i].x;
			hfinal[7+i*8].x=h1[i].x;
		}
		
		h0[0].y = H00.y;
		h0[1].y = H00.y;
		h0[2].y = H01.y;
		h0[3].y = H00.y;
		h0[4].y = H01.y;
		h0[5].y = H00.y;
		h0[6].y = H01.y;
		h0[7].y = H01.y;
		
		h1[0].y = H10.y;
		h1[1].y = H10.y;
		h1[2].y = H11.y;
		h1[3].y = H10.y;
		h1[4].y = H11.y;
		h1[5].y = H10.y;
		h1[6].y = H11.y;
		h1[7].y = H11.y;

		for (i=0; i<8; i++) {
			hfinal[0+i*8].y=h0[i].y;
			hfinal[1+i*8].y=h0[i].y;
			hfinal[2+i*8].y=h1[i].y;
			hfinal[3+i*8].y=h0[i].y;
			hfinal[4+i*8].y=h1[i].y;
			hfinal[5+i*8].y=h0[i].y;
			hfinal[6+i*8].y=h1[i].y;
			hfinal[7+i*8].y=h1[i].y;
		}
		
		h0[0].z = H00.z;
		h0[1].z = H00.z;
		h0[2].z = H00.z;
		h0[3].z = H01.z;
		h0[4].z = H00.z;
		h0[5].z = H01.z;
		h0[6].z = H01.z;
		h0[7].z = H01.z;
		
		h1[0].z = H10.z;
		h1[1].z = H10.z;
		h1[2].z = H10.z;
		h1[3].z = H11.z;
		h1[4].z = H10.z;
		h1[5].z = H11.z;
		h1[6].z = H11.z;
		h1[7].z = H11.z;
		
		for (i=0; i<8; i++) {
			hfinal[0+i*8].z=h0[i].z;
			hfinal[1+i*8].z=h0[i].z;
			hfinal[2+i*8].z=h0[i].z;
			hfinal[3+i*8].z=h1[i].z;
			hfinal[4+i*8].z=h0[i].z;
			hfinal[5+i*8].z=h1[i].z;
			hfinal[6+i*8].z=h1[i].z;
			hfinal[7+i*8].z=h1[i].z;
		}
		
		value=0.0;
		
		for (i=0; i<64; i++) {
			value+=X[i]*hfinal[i].x*hfinal[i].y*hfinal[i].z;
		}
		
		return value;
		
	}
	for (i=0;i<64;i++)
	{
		a[i]=(double)(0.0);
		for (j=0;j<64;j++)
			a[i]+=Coeff[i][j]*X[j];
	}
	
	// energy
	value=0.0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			for (k=0;k<4;k++)
				value+=a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j)*pow(r.z,k);
	
	return value;
}

int InterpolateDelaunayGrid(int typeA,VECTOR pos)
{
	int x0,y0,z0,x1,y1,z1,i;
	int X[8], InAccess;
	
	VECTOR s;
	
	if(!DelaunayGrid[typeA])
	{
		printf("InterpolateDelaunayGrid: Pseudoatom %d [%s] not tabulated.\n",typeA,PseudoAtoms[typeA].Name);
		exit(0);
	}
	
	switch(BoundaryCondition[CurrentSystem])
	{
		case RECTANGULAR:
			// the position has to be moved back to the main unit cell using the rectangular boundary condition
			pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
			pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
			pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
			if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
			if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
			if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
			break;
		case TRICLINIC:
		default:
			// the position first has to be moved back to the main unit cell using the triclinic boundary condition
			s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
			s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
			s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;
			
			s.x-=(REAL)NINT(s.x);
			s.y-=(REAL)NINT(s.y);
			s.z-=(REAL)NINT(s.z);
			
			if(s.x<0.0) s.x+=1.0;
			if(s.y<0.0) s.y+=1.0;
			if(s.z<0.0) s.z+=1.0;
			
			pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
			pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
			pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
			break;
	}
	
	s.x=((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)(NumberOfDelaunayGridPoints.x);
	s.y=((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)(NumberOfDelaunayGridPoints.y);
	s.z=((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)(NumberOfDelaunayGridPoints.z);
	
	// find the corresponding cell, between 0 and NumberOfVDWGridPoints
	// determine lower boundary
	x0=(int)s.x;
	y0=(int)s.y;
	z0=(int)s.z;
	
	// determine upper boundary, apply periodic boundary condition
	x1=x0+1;
	y1=y0+1;
	z1=z0+1;
	
	// retrieve the grid value
	X[0]=DelaunayGrid[typeA][x0][y0][z0][0];
	X[1]=DelaunayGrid[typeA][x1][y0][z0][0];
	X[2]=DelaunayGrid[typeA][x0][y1][z0][0];
	X[3]=DelaunayGrid[typeA][x1][y1][z0][0];
	X[4]=DelaunayGrid[typeA][x0][y0][z1][0];
	X[5]=DelaunayGrid[typeA][x1][y0][z1][0];
	X[6]=DelaunayGrid[typeA][x0][y1][z1][0];
	X[7]=DelaunayGrid[typeA][x1][y1][z1][0];
	
	i = 0;
	do {
		InAccess =X[i];
		i++;
	} while (InAccess > 0.5 && i<8);
	
	return InAccess;
}

int InterpolateDynamicGrid(VECTOR pos)
{
	int x0,y0,z0,x1,y1,z1,i;
	int X[8], InAccess;
	
	VECTOR s;
	
	switch(BoundaryCondition[CurrentSystem])
	{
		case RECTANGULAR:
			// the position has to be moved back to the main unit cell using the rectangular boundary condition
			pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
			pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
			pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
			if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
			if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
			if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
			break;
		case TRICLINIC:
		default:
			// the position first has to be moved back to the main unit cell using the triclinic boundary condition
			s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
			s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
			s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;
			
			s.x-=(REAL)NINT(s.x);
			s.y-=(REAL)NINT(s.y);
			s.z-=(REAL)NINT(s.z);
			
			if(s.x<0.0) s.x+=1.0;
			if(s.y<0.0) s.y+=1.0;
			if(s.z<0.0) s.z+=1.0;
			
			pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
			pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
			pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
			break;
	}
	
	s.x=((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)(NumberOfDynamicGridPoints.x);
	s.y=((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)(NumberOfDynamicGridPoints.y);
	s.z=((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)(NumberOfDynamicGridPoints.z);
	
	// find the corresponding cell, between 0 and NumberOfVDWGridPoints
	// determine lower boundary
	x0=(int)s.x;
	y0=(int)s.y;
	z0=(int)s.z;
	
	// determine upper boundary, apply periodic boundary condition
	x1=x0+1;
	y1=y0+1;
	z1=z0+1;
	
	// retrieve the cube grid values
	X[0]=DynamicGrid[x0][y0][z0][0];
	X[1]=DynamicGrid[x1][y0][z0][0];
	X[2]=DynamicGrid[x0][y1][z0][0];
	X[3]=DynamicGrid[x1][y1][z0][0];
	X[4]=DynamicGrid[x0][y0][z1][0];
	X[5]=DynamicGrid[x1][y0][z1][0];
	X[6]=DynamicGrid[x0][y1][z1][0];
	X[7]=DynamicGrid[x1][y1][z1][0];
	
	i = 0;
	while (i<8)
	{
		InAccess = X[i];
		i++;
		if (InAccess != 0 && InAccess != Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Grid) {
			return InAccess;
		}
	}
	
	return InAccess;
}

REAL InterpolateVDWForceGrid(int typeA,VECTOR pos,VECTOR *Force)
{
  int i,j,k;
  int x0,y0,z0,x1,y1,z1;
  VECTOR s,r;
  REAL value;
  REAL X[64],a[64];

  if(!VDWGrid[typeA])
  {
    printf("InterpolateVDWGrid: Pseudoatom %d [%s] not tabulated.\n",typeA,PseudoAtoms[typeA].Name);
    exit(0);
  }

  switch(BoundaryCondition[CurrentSystem])
  {
    case RECTANGULAR:
      // the position has to be moved back to the main unit cell using the rectangular boundary condition
      pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
      pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
      pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
      if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
      if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
      if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
      break;
    case TRICLINIC:
    default:
      // the position first has to be moved back to the main unit cell using the triclinic boundary condition
      s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
      s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
      s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

      s.x-=(REAL)NINT(s.x);
      s.y-=(REAL)NINT(s.y);
      s.z-=(REAL)NINT(s.z);

      if(s.x<0.0) s.x+=1.0;
      if(s.y<0.0) s.y+=1.0;
      if(s.z<0.0) s.z+=1.0;

      pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
  }

  s.x=((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)NumberOfVDWGridPoints.x;
  s.y=((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)NumberOfVDWGridPoints.y;
  s.z=((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)NumberOfVDWGridPoints.z;

  // find the corresponding cell, between 0 and NumberOfVDWGridPoints
  // determine lower boundary
  x0=(int)s.x;
  y0=(int)s.y;
  z0=(int)s.z;

  // determine upper boundary, apply periodic boundary condition
  x1=x0+1;
  y1=y0+1;
  z1=z0+1;

  // find the corresponding position within that cell (between 0.0 and 1.0)
  r.x=(s.x-x0);
  r.y=(s.y-y0);
  r.z=(s.z-z0);

  // retrieve the grid value and the derivatives
  for(i=0;i<8;i++)
  {
    X[0+i*8]=VDWGrid[typeA][x0][y0][z0][i];
    X[1+i*8]=VDWGrid[typeA][x1][y0][z0][i];
    X[2+i*8]=VDWGrid[typeA][x0][y1][z0][i];
    X[3+i*8]=VDWGrid[typeA][x1][y1][z0][i];
    X[4+i*8]=VDWGrid[typeA][x0][y0][z1][i];
    X[5+i*8]=VDWGrid[typeA][x1][y0][z1][i];
    X[6+i*8]=VDWGrid[typeA][x0][y1][z1][i];
    X[7+i*8]=VDWGrid[typeA][x1][y1][z1][i];
  }

  // cap values to avoid a negative interpolation result
  for(i=0;i<8;i++)
    if(X[i]>1e10) return 1e8;

  for (i=0;i<64;i++)
  {
    a[i]=(double)(0.0);
    for (j=0;j<64;j++)
      a[i]+=Coeff[i][j]*X[j];
  }

  value=0.0;
  for (i=1;i<4;i++)
    for (j=0;j<4;j++)
      for (k=0;k<4;k++) 
        value+=i*a[i+4*j+16*k]*pow(r.x,i-1)*pow(r.y,j)*pow(r.z,k);
  Force->x=-value/DeltaVDWGrid.x;

  value=0.0;
  for (i=0;i<4;i++)
    for (j=1;j<4;j++)
      for (k=0;k<4;k++)
        value+=j*a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j-1)*pow(r.z,k);
  Force->y=-value/DeltaVDWGrid.y;

  value=0.0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      for (k=1;k<4;k++)
        value+=k*a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j)*pow(r.z,k-1);
  Force->z=-value/DeltaVDWGrid.z;

  // energy
  value=0.0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      for (k=0;k<4;k++)
        value+=a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j)*pow(r.z,k);

  return value;
}

REAL InterpolateCoulombGrid(int typeA,VECTOR pos)
{
  int i,j,k;
  int x0,y0,z0,x1,y1,z1;
  VECTOR s,r;
  REAL value;
  REAL X[64],a[64];

  if(!CoulombGrid)
  {
    printf("InterpolateCoulombGrid: Pseudoatom %d [%s] not tabulated.\n",typeA,PseudoAtoms[typeA].Name);
    exit(0);
  }

  switch(BoundaryCondition[CurrentSystem])
  {
    case RECTANGULAR:
      // the position has to be moved back to the main unit cell using the rectangular boundary condition
      pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
      pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
      pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
      if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
      if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
      if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
      break;
    case TRICLINIC:
    default:
      // the position first has to be moved back to the main unit cell using the triclinic boundary condition
      s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
      s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
      s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

      s.x-=(REAL)NINT(s.x);
      s.y-=(REAL)NINT(s.y);
      s.z-=(REAL)NINT(s.z);

      if(s.x<0.0) s.x+=1.0;
      if(s.y<0.0) s.y+=1.0;
      if(s.z<0.0) s.z+=1.0;

      pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
  }

  s.x=((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)NumberOfCoulombGridPoints.x;
  s.y=((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)NumberOfCoulombGridPoints.y;
  s.z=((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)NumberOfCoulombGridPoints.z;


  // find the corresponding cell, between 0 and NumberOfCoulombGridPoints
  // determine lower boundary
  x0=(int)s.x;
  y0=(int)s.y;
  z0=(int)s.z;

  // determine upper boundary, apply periodic boundary condition
  x1=x0+1;
  y1=y0+1;
  z1=z0+1;

  // find the corresponding position within that cell (between 0.0 and 1.0)
  r.x=(s.x-x0);
  r.y=(s.y-y0);
  r.z=(s.z-z0);

  // retrieve the grid value and the derivatives
  for(i=0;i<8;i++)
  {
    X[0+i*8]=CoulombGrid[x0][y0][z0][i];
    X[1+i*8]=CoulombGrid[x1][y0][z0][i];
    X[2+i*8]=CoulombGrid[x0][y1][z0][i];
    X[3+i*8]=CoulombGrid[x1][y1][z0][i];
    X[4+i*8]=CoulombGrid[x0][y0][z1][i];
    X[5+i*8]=CoulombGrid[x1][y0][z1][i];
    X[6+i*8]=CoulombGrid[x0][y1][z1][i];
    X[7+i*8]=CoulombGrid[x1][y1][z1][i];
  }

  for (i=0;i<64;i++)
  {
    a[i]=(double)(0.0);
    for (j=0;j<64;j++)
      a[i]+=Coeff[i][j]*X[j];
  }

  // energy
  value=0.0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      for (k=0;k<4;k++)
        value+=a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j)*pow(r.z,k);

  return PseudoAtoms[typeA].Charge1*value;
}

REAL InterpolateCoulombForceGrid(int typeA,VECTOR pos,VECTOR *Force)
{
  int i,j,k;
  int x0,y0,z0,x1,y1,z1;
  VECTOR s,r;
  REAL value;
  REAL X[64],a[64];

  if(!CoulombGrid)
  {
    printf("InterpolateCoulombGrid: Pseudoatom %d [%s] not tabulated.\n",typeA,PseudoAtoms[typeA].Name);
    exit(0);
  }

  switch(BoundaryCondition[CurrentSystem])
  {
    case RECTANGULAR:
      // the position has to be moved back to the main unit cell using the rectangular boundary condition
      pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
      pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
      pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
      if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
      if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
      if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
      break;
    case TRICLINIC:
    default:
      // the position first has to be moved back to the main unit cell using the triclinic boundary condition
      s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
      s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
      s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

      s.x-=(REAL)NINT(s.x);
      s.y-=(REAL)NINT(s.y);
      s.z-=(REAL)NINT(s.z);

      if(s.x<0.0) s.x+=1.0;
      if(s.y<0.0) s.y+=1.0;
      if(s.z<0.0) s.z+=1.0;

      pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
  } 

  s.x=((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)NumberOfCoulombGridPoints.x;
  s.y=((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)NumberOfCoulombGridPoints.y;
  s.z=((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)NumberOfCoulombGridPoints.z;

  // find the corresponding cell, between 0 and NumberOfCoulombGridPoints
  // determine lower boundary
  x0=(int)s.x;
  y0=(int)s.y;
  z0=(int)s.z;

  // determine upper boundary, apply periodic boundary condition
  x1=x0+1;
  y1=y0+1;
  z1=z0+1;

  // find the corresponding position within that cell (between 0.0 and 1.0)
  r.x=(s.x-x0);
  r.y=(s.y-y0);
  r.z=(s.z-z0);

  // retrieve the grid value and the derivatives
  for(i=0;i<8;i++)
  {
    X[0+i*8]=CoulombGrid[x0][y0][z0][i];
    X[1+i*8]=CoulombGrid[x1][y0][z0][i];
    X[2+i*8]=CoulombGrid[x0][y1][z0][i];
    X[3+i*8]=CoulombGrid[x1][y1][z0][i];
    X[4+i*8]=CoulombGrid[x0][y0][z1][i];
    X[5+i*8]=CoulombGrid[x1][y0][z1][i];
    X[6+i*8]=CoulombGrid[x0][y1][z1][i];
    X[7+i*8]=CoulombGrid[x1][y1][z1][i];
  }

  for (i=0;i<64;i++)
  {
    a[i]=(double)(0.0);
    for (j=0;j<64;j++)
      a[i]+=Coeff[i][j]*X[j];
  }

  value=0.0;
  for (i=1;i<4;i++)
    for (j=0;j<4;j++)
      for (k=0;k<4;k++) 
        value+=i*a[i+4*j+16*k]*pow(r.x,i-1)*pow(r.y,j)*pow(r.z,k);
  Force->x=-PseudoAtoms[typeA].Charge1*value/DeltaCoulombGrid.x;

  value=0.0;
  for (i=0;i<4;i++)
    for (j=1;j<4;j++)
      for (k=0;k<4;k++)
        value+=j*a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j-1)*pow(r.z,k);
  Force->y=-PseudoAtoms[typeA].Charge1*value/DeltaCoulombGrid.y;

  value=0.0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      for (k=1;k<4;k++)
        value+=k*a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j)*pow(r.z,k-1);
  Force->z=-PseudoAtoms[typeA].Charge1*value/DeltaCoulombGrid.z;

  // energy
  value=0.0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      for (k=0;k<4;k++)
        value+=a[i+4*j+16*k]*pow(r.x,i)*pow(r.y,j)*pow(r.z,k);

  return PseudoAtoms[typeA].Charge1*value;
}

void TestForceGrid(FILE *FilePtr)
{
  int i,j,k,typeA;
  REAL resf[8],rest[8],ebolt[8],bolt[8],bolf[8],ebolf[8],errb[8],eeb[8];
  REAL Vlj,Vcoul;
  VECTOR Force,CoulombForce;
  POINT pos,s;
  int StoreBiasingMethod,already_present;

  StoreBiasingMethod=BiasingMethod;
  BiasingMethod=LJ_AND_REAL_BIASING;

  // compute the modified VDW sites for anisotropic models
  CalculateAnisotropicSites();

  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Number Of Gridpoints for VDW-interactions: %d\n",
    NumberOfVDWGridPoints.x*NumberOfVDWGridPoints.y*NumberOfVDWGridPoints.z);
  fprintf(FilePtr,"=============================================================\n");
  fprintf(FilePtr,"NumberOfGridPoints [x]: %d\n",NumberOfVDWGridPoints.x);
  fprintf(FilePtr,"NumberOfGridPoints [y]: %d\n",NumberOfVDWGridPoints.y);
  fprintf(FilePtr,"NumberOfGridPoints [z]: %d\n",NumberOfVDWGridPoints.z);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Number Of Gridpoints for Coulomb-interactions: %d\n",
    NumberOfCoulombGridPoints.x*NumberOfCoulombGridPoints.y*NumberOfCoulombGridPoints.z);
  fprintf(FilePtr,"=============================================================\n");
  fprintf(FilePtr,"NumberOfGridPoints [x]: %d\n",NumberOfCoulombGridPoints.x);
  fprintf(FilePtr,"NumberOfGridPoints [y]: %d\n",NumberOfCoulombGridPoints.y);
  fprintf(FilePtr,"NumberOfGridPoints [z]: %d\n",NumberOfCoulombGridPoints.z);
  fprintf(FilePtr,"\n");

  for(typeA=0;typeA<NumberOfPseudoAtoms;typeA++)
  {
    if(VDWGrid[typeA])
    {
      for(j=0;j<8;j++)
      {
        bolt[j]=0.0;
        ebolt[j]=0.0;
        bolf[j]=0.0; 
        ebolf[j]=0.0; 
        errb[j]=0.0;
        eeb[j]=0.0;
      }

      for(i=0;i<5000;i++)
      {
        switch(BoundaryCondition[CurrentSystem])
        {
          case RECTANGULAR:
            pos.x=Box[CurrentSystem].ax*RandomNumber();
            pos.y=Box[CurrentSystem].by*RandomNumber();
            pos.z=Box[CurrentSystem].cz*RandomNumber();
            break;
          case TRICLINIC:
            s.x=RandomNumber();
            s.y=RandomNumber();
            s.z=RandomNumber();
            pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
            pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
            pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;
            break;
        }

        Framework[CurrentSystem].FrameworkModel=FULL;
        CalculateFrameworkForceAtPosition(pos,typeA,&Vlj,&Force,&Vcoul,&CoulombForce);

        // problem: to which component belongs this type?
        //if(BlockedPocketGridTest(pos))
        //  Vlj=Force.x=Force.y=Force.z=EnergyOverlapCriteria;
 
        resf[0]=Vlj;
        resf[1]=Force.x;
        resf[2]=Force.y;
        resf[3]=Force.z;

        resf[4]=Vcoul;
        resf[5]=CoulombForce.x;
        resf[6]=CoulombForce.y;
        resf[7]=CoulombForce.z;

        Framework[CurrentSystem].FrameworkModel=GRID;
        CalculateFrameworkForceAtPosition(pos,typeA,&Vlj,&Force,&Vcoul,&CoulombForce);

        rest[0]=Vlj;
        rest[1]=Force.x;
        rest[2]=Force.y;
        rest[3]=Force.z;

        rest[4]=Vcoul;
        rest[5]=CoulombForce.x;
        rest[6]=CoulombForce.y;
        rest[7]=CoulombForce.z;

        // VDW
        for(j=0;j<8;j++)
        {
          bolt[j]+=exp(-rest[0]*Beta[CurrentSystem]);
          ebolt[j]+=rest[j]*exp(-rest[0]*Beta[CurrentSystem]);

          bolf[j]+=exp(-resf[0]*Beta[CurrentSystem]);
          ebolf[j]+=resf[j]*exp(-resf[0]*Beta[CurrentSystem]);

          errb[j]+=exp(-resf[0]*Beta[CurrentSystem])*SQR(resf[j]-rest[j]);
          eeb[j]+=exp(-resf[0]*Beta[CurrentSystem])*SQR(resf[j]);
        }
      }

      fprintf(FilePtr,"\n\n");
      fprintf(FilePtr,"PseudoAtom %d Framework-[%s]\n",typeA,PseudoAtoms[typeA].Name);

      fprintf(FilePtr,"=========================================================================================\n");
      fprintf(FilePtr,"\tBoltzmann average energy VDW (table)                 : %18.12lf\n",(double)(ebolt[0]/bolt[0]));
      fprintf(FilePtr,"\tBoltzmann average energy VDW (full)                  : %18.12lf\n",(double)(ebolf[0]/bolf[0]));
      fprintf(FilePtr,"\tBoltzmann relative error VDW                         : %18.12lf\n",(double)(fabs(ebolt[0])<1e-8?0.0:sqrt(errb[0]/eeb[0])));
      if(fabs(ebolt[0])<1e-8?0.0:sqrt(errb[0]/eeb[0])>0.01) fprintf(FilePtr,"\tError: VDW energy interpolation table is wrong...\n");
      if(PseudoAtoms[typeA].HasCharges)
      {
        fprintf(FilePtr,"\tBoltzmann average energy Coulomb (table)             : %18.12lf\n",(double)(ebolt[4]/bolt[4]));
        fprintf(FilePtr,"\tBoltzmann average energy Coulomb (full)              : %18.12lf\n",(double)(ebolf[4]/bolf[4]));
        fprintf(FilePtr,"\tBoltzmann relative error Coulomb                     : %18.12lf\n",(double)(fabs(ebolt[4])<1e-8?0.0:sqrt(errb[4]/eeb[4])));
        if(fabs(ebolt[4])<1e-8?0.0:sqrt(errb[4]/eeb[4])>0.01) fprintf(FilePtr,"\tError: Coulomb energy interpolation table is wrong...\n");
      }

      if((fabs(ebolt[0])<1e-8?0.0:sqrt(errb[0]/eeb[0])>0.01)||
         (PseudoAtoms[typeA].HasCharges&&(fabs(ebolt[4])<1e-8?0.0:sqrt(errb[4]/eeb[4])>0.01)))
      {
        already_present=FALSE;
        for(k=0;k<NumberOfWarnings[CurrentSystem];k++)
          if(Warnings[CurrentSystem][k]==GRID_ERROR_ENERGY) already_present=TRUE;
        if((!ContinueAfterCrash)&&(!already_present))
        {
          if(NumberOfWarnings[k]<MAX_NUMBER_OF_WARNINGS)
          {
            Warnings[CurrentSystem][NumberOfWarnings[CurrentSystem]]=GRID_ERROR_ENERGY;
            NumberOfWarnings[CurrentSystem]++;
          }
        }
      }


      fprintf(FilePtr,"=========================================================================================\n");
      fprintf(FilePtr,"\tBoltzmann average Force[x] VDW (table)               : %18.12lf\n",(double)(ebolt[1]/bolt[1]));
      fprintf(FilePtr,"\tBoltzmann average Force[x] VDW (full)                : %18.12lf\n",(double)(ebolf[1]/bolf[1]));
      fprintf(FilePtr,"\tBoltzmann relative error VDW                         : %18.12lf\n",(double)(fabs(ebolt[1])<1e-8?0.0:sqrt(errb[1]/eeb[1])));
      if(fabs(ebolt[1])<1e-8?0.0:sqrt(errb[1]/eeb[1])>0.01) fprintf(FilePtr,"\tError: VDW Force[x] interpolation table is wrong...\n");
      if(PseudoAtoms[typeA].HasCharges)
      {
        fprintf(FilePtr,"\tBoltzmann average Force[x] Coulomb (table)           : %18.12lf\n",(double)(ebolt[5]/bolt[5]));
        fprintf(FilePtr,"\tBoltzmann average Force[x] Coulomb (full)            : %18.12lf\n",(double)(ebolf[5]/bolf[5]));
        fprintf(FilePtr,"\tBoltzmann relative error Coulomb                     : %18.12lf\n",(double)(fabs(ebolt[5])<1e-8?0.0:sqrt(errb[5]/eeb[5])));
        if(fabs(ebolt[5])<1e-8?0.0:sqrt(errb[5]/eeb[5])>0.01) fprintf(FilePtr,"\tError: Coulomb Force[z] interpolation table is wrong...\n");
      }
      fprintf(FilePtr,"=========================================================================================\n");
      fprintf(FilePtr,"\tBoltzmann average Force[y] VDW (table)               : %18.12lf\n",(double)(ebolt[2]/bolt[2]));
      fprintf(FilePtr,"\tBoltzmann average Force[y] VDW (full)                : %18.12lf\n",(double)(ebolf[2]/bolf[2]));
      fprintf(FilePtr,"\tBoltzmann relative error VDW                         : %18.12lf\n",(double)(fabs(ebolt[2])<1e-8?0.0:sqrt(errb[2]/eeb[2])));
      if(fabs(ebolt[2])<1e-8?0.0:sqrt(errb[2]/eeb[2])>0.01) fprintf(FilePtr,"\tError: VDW Force[y] interpolation table is wrong...\n");
      if(PseudoAtoms[typeA].HasCharges)
      {
        fprintf(FilePtr,"\tBoltzmann average Force[y] Coulomb (table)           : %18.12lf\n",(double)(ebolt[6]/bolt[6]));
        fprintf(FilePtr,"\tBoltzmann average Force[y] Coulomb (full)            : %18.12lf\n",(double)(ebolf[6]/bolf[6]));
        fprintf(FilePtr,"\tBoltzmann relative error Coulomb                     : %18.12lf\n",(double)(fabs(ebolt[6])<1e-8?0.0:sqrt(errb[6]/eeb[6])));
        if(fabs(ebolt[6])<1e-8?0.0:sqrt(errb[6]/eeb[6])>0.01) fprintf(FilePtr,"\tError: Coulomb Force[y] interpolation table is wrong...\n");
      }
      fprintf(FilePtr,"=========================================================================================\n");
      fprintf(FilePtr,"\tBoltzmann average Force[z] VDW (table)               : %18.12lf\n",(double)(ebolt[3]/bolt[3]));
      fprintf(FilePtr,"\tBoltzmann average Force[z] VDW (full)                : %18.12lf\n",(double)(ebolf[3]/bolf[3]));
      fprintf(FilePtr,"\tBoltzmann relative error VDW                         : %18.12lf\n",(double)(fabs(ebolt[3])<1e-8?0.0:sqrt(errb[3]/eeb[3])));
      if(fabs(ebolt[3])<1e-8?0.0:sqrt(errb[3]/eeb[3])>0.01) fprintf(FilePtr,"\tError: VDW Force[z] interpolation table is wrong...\n");
      if(PseudoAtoms[typeA].HasCharges)
      {
        fprintf(FilePtr,"\tBoltzmann average Force[z] Coulomb (table)           : %18.12lf\n",(double)(ebolt[7]/bolt[7]));
        fprintf(FilePtr,"\tBoltzmann average Force[z] Coulomb (full)            : %18.12lf\n",(double)(ebolf[7]/bolf[7]));
        fprintf(FilePtr,"\tBoltzmann relative error Coulomb                     : %18.12lf\n",(double)(fabs(ebolt[7])<1e-8?0.0:sqrt(errb[7]/eeb[7])));
        if(fabs(ebolt[7])<1e-8?0.0:sqrt(errb[7]/eeb[7])>0.01) fprintf(FilePtr,"\tError: VDW Force[z] interpolation table is wrong...\n");
      }

      if((fabs(ebolt[1])<1e-8?0.0:sqrt(errb[1]/eeb[1])>0.01)||
         (fabs(ebolt[2])<1e-8?0.0:sqrt(errb[2]/eeb[2])>0.01)||
         (fabs(ebolt[3])<1e-8?0.0:sqrt(errb[3]/eeb[3])>0.01)||
         (PseudoAtoms[typeA].HasCharges&&(
         (fabs(ebolt[5])<1e-8?0.0:sqrt(errb[5]/eeb[5])>0.01)||
         (fabs(ebolt[6])<1e-8?0.0:sqrt(errb[6]/eeb[6])>0.01)||
         (fabs(ebolt[7])<1e-8?0.0:sqrt(errb[7]/eeb[7])>0.01))))
      {
        already_present=FALSE;
        for(k=0;k<NumberOfWarnings[CurrentSystem];k++)
          if(Warnings[CurrentSystem][k]==GRID_ERROR_FORCE) already_present=TRUE;
        if((!ContinueAfterCrash)&&(!already_present))
        {
          if(NumberOfWarnings[k]<MAX_NUMBER_OF_WARNINGS)
          {
            Warnings[CurrentSystem][NumberOfWarnings[CurrentSystem]]=GRID_ERROR_FORCE;
            NumberOfWarnings[CurrentSystem]++;
          }
        }
      }
    }
  }
  BiasingMethod=StoreBiasingMethod;
  Framework[CurrentSystem].FrameworkModel=GRID;
}

INT_VECTOR3 ConvertXYZPositionToGridIndex(VECTOR pos)
{
  VECTOR s;
  INT_VECTOR3 t;

  switch(BoundaryCondition[CurrentSystem])
  {
    case RECTANGULAR:
      // the position has to be moved back to the main unit cell using the rectangular boundary condition
      pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
      pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
      pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
      if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
      if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
      if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
      break;
    case TRICLINIC:
    default:
      // the position first has to be moved back to the main unit cell using the triclinic boundary condition
      s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
      s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
      s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;

      s.x-=(REAL)NINT(s.x);
      s.y-=(REAL)NINT(s.y);
      s.z-=(REAL)NINT(s.z);

      if(s.x<0.0) s.x+=1.0;
      if(s.y<0.0) s.y+=1.0;
      if(s.z<0.0) s.z+=1.0;

      pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
      pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
      pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
      break;
  }

  t.x=rint(((pos.x-ShiftGrid.x)/SizeGrid.x)*(REAL)(NumberOfVDWGridPoints.x));
  t.y=rint(((pos.y-ShiftGrid.y)/SizeGrid.y)*(REAL)(NumberOfVDWGridPoints.y));
  t.z=rint(((pos.z-ShiftGrid.z)/SizeGrid.z)*(REAL)(NumberOfVDWGridPoints.z));
  return t;
}

VECTOR ConvertGridIndexToXYZIndex(INT_VECTOR3 GridIndex)
{
  VECTOR pos;

  pos.x=GridIndex.x*SizeGrid.x/NumberOfVDWGridPoints.x+ShiftGrid.x;
  pos.y=GridIndex.y*SizeGrid.y/NumberOfVDWGridPoints.y+ShiftGrid.y;
  pos.z=GridIndex.z*SizeGrid.z/NumberOfVDWGridPoints.z+ShiftGrid.z;

  return pos;
}

VECTOR ConvertCBMCGridPosToXYZ(VECTOR pos)
{
	INT_VECTOR3 t;
	VECTOR s;
	
	pos.x=pos.x*SizeGrid.x/(REAL)NumberOfVDWGridPoints.x+ShiftGrid.x;
	pos.y=pos.y*SizeGrid.y/(REAL)NumberOfVDWGridPoints.y+ShiftGrid.y;
	pos.z=pos.z*SizeGrid.z/(REAL)NumberOfVDWGridPoints.z+ShiftGrid.z;
	
	if (pos.x <0.0 || pos.y < 0.0 || pos.z < 0.0) {
		
		switch(BoundaryCondition[CurrentSystem])
		{
			case RECTANGULAR:
				// the position has to be moved back to the main unit cell using the rectangular boundary condition
				pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
				pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
				pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
				if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
				if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
				if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
				break;
			case TRICLINIC:
			default:
				// the position first has to be moved back to the main unit cell using the triclinic boundary condition
				s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
				s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
				s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;
				
				s.x-=(REAL)NINT(s.x);
				s.y-=(REAL)NINT(s.y);
				s.z-=(REAL)NINT(s.z);
				
				if(s.x<0.0) s.x+=1.0;
				if(s.y<0.0) s.y+=1.0;
				if(s.z<0.0) s.z+=1.0;
				
				pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
				pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
				pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
				break;
		}
	}

	if (NumberOfUnitCells[CurrentSystem].x>1 || NumberOfUnitCells[CurrentSystem].y>1 || NumberOfUnitCells[CurrentSystem].z>1 ) {
		
		t.x = (int)(NumberOfUnitCells[CurrentSystem].x*RandomNumber());
		t.y = (int)(NumberOfUnitCells[CurrentSystem].y*RandomNumber());
		t.z = (int)(NumberOfUnitCells[CurrentSystem].z*RandomNumber());
		
		if (t.x == NumberOfUnitCells[CurrentSystem].x || t.y == NumberOfUnitCells[CurrentSystem].y || t.z ==NumberOfUnitCells[CurrentSystem].z)
			printf("t = %i %i %i\n",t.x,t.y,t.z);
			
		pos.x+=UnitCellBox[CurrentSystem].ax*t.x+UnitCellBox[CurrentSystem].bx*t.y+UnitCellBox[CurrentSystem].cx*t.z;
		pos.y+=UnitCellBox[CurrentSystem].ay*t.x+UnitCellBox[CurrentSystem].by*t.y+UnitCellBox[CurrentSystem].cy*t.z;
		pos.z+=UnitCellBox[CurrentSystem].az*t.x+UnitCellBox[CurrentSystem].bz*t.y+UnitCellBox[CurrentSystem].cz*t.z;
		
	}
	return pos;
}

void FloodFillNonRecursive(int seedx,int seedy,int seedz,int BlockingValue)
{
  int i,j,k;
  int w,e;
  int index;

  QueueSize=1;
  Queue[0].x=seedx;
  Queue[0].y=seedy;
  Queue[0].z=seedz;

  do
  {
    i=Queue[0].x;
    j=Queue[0].y;
    k=Queue[0].z;

    Queue[0]=Queue[QueueSize-1];
    QueueSize--;

    if(BlockingGrid[i][j][k]==0)
    {
      for(w=i;w>=0;w--)
        if(BlockingGrid[w][j][k]!=0) break;
      for(e=i;e<=NumberOfVDWGridPoints.x;e++)
        if(BlockingGrid[e][j][k]!=0) break;

      for(i=w+1;i<e;i++)
      {
        BlockingGrid[i][j][k]=BlockingValue;

        index=j+1;
        if(index>NumberOfVDWGridPoints.y) index=0;
        //if(index<=NumberOfVDWGridPoints.y)
        {
          if(BlockingGrid[i][index][k]==0)
          {
            Queue[QueueSize].x=i;
            Queue[QueueSize].y=index;
            Queue[QueueSize].z=k;
            QueueSize++;
          }
        }
        index=j-1;
        if(index<0) index=NumberOfVDWGridPoints.y;
        //if(index>=0)
        {
          if(BlockingGrid[i][index][k]==0)
          {
            Queue[QueueSize].x=i;
            Queue[QueueSize].y=index;
            Queue[QueueSize].z=k;
            QueueSize++;
          }
        }

        index=k+1;
        if(index>NumberOfVDWGridPoints.z) index=0;
        //if(index<=NumberOfVDWGridPoints.z)
        {
          if(BlockingGrid[i][j][index]==0)
          {
            Queue[QueueSize].x=i;
            Queue[QueueSize].y=j;
            Queue[QueueSize].z=index;
            QueueSize++;
          }
        }
        index=k-1;
        if(index<0) index=NumberOfVDWGridPoints.z;
        //if(index>=0)
        {
          if(BlockingGrid[i][j][index]==0)
          {
            Queue[QueueSize].x=i;
            Queue[QueueSize].y=j;
            Queue[QueueSize].z=index;
            QueueSize++;
          }
        }

      }

      // handle periodic boundaries x
      index=w+1;
      if(index>NumberOfVDWGridPoints.x)
      {
        index=0;
        if(BlockingGrid[index][j][k]==0)
        {
          Queue[QueueSize].x=index;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
      }
      index=w-1;
      if(index<0) 
      {
        index=NumberOfVDWGridPoints.x;
        if(BlockingGrid[index][j][k]==0)
        {
          Queue[QueueSize].x=index;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
      }
    }
  }
  while(QueueSize>0);
}

int CheckNumber(VECTOR SeedPoint,int ID)
{
  VECTOR C,vec,NewSeedPoint;
  INT_VECTOR3 GridPoint; 
  int x,y,z;
  int i,j,k;

  C.x=1.0;
  C.y=0.0;
  C.z=0.0;
  vec.x=UnitCellBox[0].ax*C.x+UnitCellBox[0].bx*C.y+UnitCellBox[0].cx*C.z;
  vec.y=UnitCellBox[0].ay*C.x+UnitCellBox[0].by*C.y+UnitCellBox[0].cy*C.z;
  vec.z=UnitCellBox[0].az*C.x+UnitCellBox[0].bz*C.y+UnitCellBox[0].cz*C.z;

  NewSeedPoint.x=SeedPoint.x+vec.x;
  NewSeedPoint.y=SeedPoint.y+vec.y;
  NewSeedPoint.z=SeedPoint.z+vec.z;
  GridPoint=ConvertXYZPositionToGridIndex(NewSeedPoint);

  for(i=-1;i<=1;i++)
    for(j=-1;j<=1;j++)
       for(k=-1;k<=1;k++)
       {
         x=GridPoint.x+i;
         y=GridPoint.x+j;
         z=GridPoint.x+k;
         if((x>=0&&x<=NumberOfVDWGridPoints.x)&&(y>=0&&y<=NumberOfVDWGridPoints.y)&&(z>=0&&z<=NumberOfVDWGridPoints.z))
             if(BlockingGrid[x][y][z]>2) return BlockingGrid[x][y][z];
       }

  NewSeedPoint.x=SeedPoint.x-vec.x;
  NewSeedPoint.y=SeedPoint.y-vec.y;
  NewSeedPoint.z=SeedPoint.z-vec.z;
  GridPoint=ConvertXYZPositionToGridIndex(NewSeedPoint);
  for(i=-1;i<=1;i++)
    for(j=-1;j<=1;j++)
       for(k=-1;k<=1;k++)
       {
         x=GridPoint.x+i;
         y=GridPoint.x+j;
         z=GridPoint.x+k;
         if((x>=0&&x<=NumberOfVDWGridPoints.x)&&(y>=0&&y<=NumberOfVDWGridPoints.y)&&(z>=0&&z<=NumberOfVDWGridPoints.z))
             if(BlockingGrid[x][y][z]>2) return BlockingGrid[x][y][z];
       }

  C.x=0.0;
  C.y=1.0;
  C.z=0.0;
  vec.x=UnitCellBox[0].ax*C.x+UnitCellBox[0].bx*C.y+UnitCellBox[0].cx*C.z;
  vec.y=UnitCellBox[0].ay*C.x+UnitCellBox[0].by*C.y+UnitCellBox[0].cy*C.z;
  vec.z=UnitCellBox[0].az*C.x+UnitCellBox[0].bz*C.y+UnitCellBox[0].cz*C.z;

  NewSeedPoint.x=SeedPoint.x+vec.x;
  NewSeedPoint.y=SeedPoint.y+vec.y;
  NewSeedPoint.z=SeedPoint.z+vec.z;
  GridPoint=ConvertXYZPositionToGridIndex(NewSeedPoint);
  for(i=-1;i<=1;i++)
    for(j=-1;j<=1;j++)
       for(k=-1;k<=1;k++)
       {
         x=GridPoint.x+i;
         y=GridPoint.x+j;
         z=GridPoint.x+k;
         if((x>=0&&x<=NumberOfVDWGridPoints.x)&&(y>=0&&y<=NumberOfVDWGridPoints.y)&&(z>=0&&z<=NumberOfVDWGridPoints.z))
             if(BlockingGrid[x][y][z]>2) return BlockingGrid[x][y][z];
       }


  NewSeedPoint.x=SeedPoint.x-vec.x;
  NewSeedPoint.y=SeedPoint.y-vec.y;
  NewSeedPoint.z=SeedPoint.z-vec.z;
  GridPoint=ConvertXYZPositionToGridIndex(NewSeedPoint);
  for(i=-1;i<=1;i++)
    for(j=-1;j<=1;j++)
       for(k=-1;k<=1;k++)
       {
         x=GridPoint.x+i;
         y=GridPoint.x+j;
         z=GridPoint.x+k;
         if((x>=0&&x<=NumberOfVDWGridPoints.x)&&(y>=0&&y<=NumberOfVDWGridPoints.y)&&(z>=0&&z<=NumberOfVDWGridPoints.z))
             if(BlockingGrid[x][y][z]>2) return BlockingGrid[x][y][z];
       }


  C.x=0.0;
  C.y=0.0;
  C.z=1.0;
  vec.x=UnitCellBox[0].ax*C.x+UnitCellBox[0].bx*C.y+UnitCellBox[0].cx*C.z;
  vec.y=UnitCellBox[0].ay*C.x+UnitCellBox[0].by*C.y+UnitCellBox[0].cy*C.z;
  vec.z=UnitCellBox[0].az*C.x+UnitCellBox[0].bz*C.y+UnitCellBox[0].cz*C.z;

  NewSeedPoint.x=SeedPoint.x+vec.x;
  NewSeedPoint.y=SeedPoint.y+vec.y;
  NewSeedPoint.z=SeedPoint.z+vec.z;
  GridPoint=ConvertXYZPositionToGridIndex(NewSeedPoint);
  for(i=-1;i<=1;i++)
    for(j=-1;j<=1;j++)
       for(k=-1;k<=1;k++)
       {
         x=GridPoint.x+i;
         y=GridPoint.x+j;
         z=GridPoint.x+k;
         if((x>=0&&x<=NumberOfVDWGridPoints.x)&&(y>=0&&y<=NumberOfVDWGridPoints.y)&&(z>=0&&z<=NumberOfVDWGridPoints.z))
             if(BlockingGrid[x][y][z]>2) return BlockingGrid[x][y][z];
       }


  NewSeedPoint.x=SeedPoint.x-vec.x;
  NewSeedPoint.y=SeedPoint.y-vec.y;
  NewSeedPoint.z=SeedPoint.z-vec.z;
  GridPoint=ConvertXYZPositionToGridIndex(NewSeedPoint);
  for(i=-1;i<=1;i++)
    for(j=-1;j<=1;j++)
       for(k=-1;k<=1;k++)
       {
         x=GridPoint.x+i;
         y=GridPoint.x+j;
         z=GridPoint.x+k;
         if((x>=0&&x<=NumberOfVDWGridPoints.x)&&(y>=0&&y<=NumberOfVDWGridPoints.y)&&(z>=0&&z<=NumberOfVDWGridPoints.z))
             if(BlockingGrid[x][y][z]>2) return BlockingGrid[x][y][z];
       }

  return ID;
}

VECTOR FloodFillNonRecursivePockets(int seedx,int seedy,int seedz,int PocketID,int Offset)
{
  int i,j,k;
  int w,e;
  int index;
  VECTOR dr;
  REAL r;
  int NumberOfVoxels;
  VECTOR PocketCenter,SeedCenter,CurrentCenter;
  INT_VECTOR3 temp;

  temp.x=seedx;
  temp.y=seedy;
  temp.z=seedz;
  SeedCenter=ConvertGridIndexToXYZIndex(temp);

  NumberOfVoxels=0;
  PocketCenter.x=0;
  PocketCenter.y=0;
  PocketCenter.z=0;

  QueueSize=1;
  Queue[0].x=seedx;
  Queue[0].y=seedy;
  Queue[0].z=seedz;

  do
  {
    i=Queue[0].x;
    j=Queue[0].y;
    k=Queue[0].z;

    Queue[0]=Queue[QueueSize-1];
    QueueSize--;

    if(BlockingGrid[i][j][k]==PocketID)
    {
      for(w=i;w>=0;w--)
        if(BlockingGrid[w][j][k]!=PocketID) break;
      for(e=i;e<=NumberOfVDWGridPoints.x;e++)
        if(BlockingGrid[e][j][k]!=PocketID) break;

      for(i=w+1;i<e;i++)
      {
        BlockingGrid[i][j][k]=PocketID+Offset;

        temp.x=i;
        temp.y=j;
        temp.z=k;
        CurrentCenter=ConvertGridIndexToXYZIndex(temp);
        dr.x=CurrentCenter.x-SeedCenter.x;
        dr.y=CurrentCenter.y-SeedCenter.y;
        dr.z=CurrentCenter.z-SeedCenter.z;
        dr=ApplyBoundaryConditionUnitCell(dr);
        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
        

        PocketCenter.x+=(SeedCenter.x+dr.x);
        PocketCenter.y+=(SeedCenter.y+dr.y);
        PocketCenter.z+=(SeedCenter.z+dr.z);
        NumberOfVoxels++;

        index=j+1;
        if(index>NumberOfVDWGridPoints.y) index=0;
        if(BlockingGrid[i][index][k]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=index;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
        index=j-1;
        if(index<0) index=NumberOfVDWGridPoints.y;
        if(BlockingGrid[i][index][k]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=index;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
        index=k+1;
        if(index>NumberOfVDWGridPoints.z) index=0;
        if(BlockingGrid[i][j][index]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=index;
          QueueSize++;
        }
        index=k-1;
        if(index<0) index=NumberOfVDWGridPoints.z;
        if(BlockingGrid[i][j][index]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=index;
          QueueSize++;
        }

      }

      // handle periodic boundaries x
      index=w+1;
      if(index>NumberOfVDWGridPoints.x)
      {
        index=0;
        if(BlockingGrid[index][j][k]==PocketID)
        {
          Queue[QueueSize].x=index;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
      }
      index=w-1;
      if(index<0)
      {
        index=NumberOfVDWGridPoints.x;
        if(BlockingGrid[index][j][k]==PocketID)
        {
          Queue[QueueSize].x=index;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
      }
    }
  }
  while(QueueSize>0);

  PocketCenter.x/=NumberOfVoxels;
  PocketCenter.y/=NumberOfVoxels;
  PocketCenter.z/=NumberOfVoxels;
  return PocketCenter;
}

REAL FloodFillNonRecursivePocketsRadius(int seedx,int seedy,int seedz,int PocketID,VECTOR SeedCenter)
{
  int i,j,k;
  int w,e;
  int index;
  VECTOR dr;
  REAL r,MaxRadius;
  INT_VECTOR3 temp;
  VECTOR CurrentCenter;

  MaxRadius=0; 
  QueueSize=1;
  Queue[0].x=seedx;
  Queue[0].y=seedy;
  Queue[0].z=seedz;

  do
  {
    i=Queue[0].x;
    j=Queue[0].y;
    k=Queue[0].z;

    Queue[0]=Queue[QueueSize-1];
    QueueSize--;

    if(BlockingGrid[i][j][k]==PocketID)
    {
      for(w=i;w>=0;w--)
        if(BlockingGrid[w][j][k]!=PocketID) break;
      for(e=i;e<=NumberOfVDWGridPoints.x;e++)
        if(BlockingGrid[e][j][k]!=PocketID) break;

      for(i=w+1;i<e;i++)
      {
        BlockingGrid[i][j][k]=PocketID+100;

        temp.x=i;
        temp.y=j;
        temp.z=k;
        CurrentCenter=ConvertGridIndexToXYZIndex(temp);
        dr.x=CurrentCenter.x-SeedCenter.x;
        dr.y=CurrentCenter.y-SeedCenter.y;
        dr.z=CurrentCenter.z-SeedCenter.z;
        dr=ApplyBoundaryConditionUnitCell(dr);
        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
        if(r>MaxRadius) MaxRadius=r;

        index=j+1;
        if(index>NumberOfVDWGridPoints.y) index=0;
        if(BlockingGrid[i][index][k]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=index;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
        index=j-1;
        if(index<0) index=NumberOfVDWGridPoints.y;
        if(BlockingGrid[i][index][k]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=index;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
        index=k+1;
        if(index>NumberOfVDWGridPoints.z) index=0;
        if(BlockingGrid[i][j][index]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=index;
          QueueSize++;
        }
        index=k-1;
        if(index<0) index=NumberOfVDWGridPoints.z;
        if(BlockingGrid[i][j][index]==PocketID)
        {
          Queue[QueueSize].x=i;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=index;
          QueueSize++;
        }
      }

      // handle periodic boundaries x
      index=w+1;
      if(index>NumberOfVDWGridPoints.x)
      {
        index=0;
        if(BlockingGrid[index][j][k]==PocketID)
        {
          Queue[QueueSize].x=index;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
      }
      index=w-1;
      if(index<0)
      {
        index=NumberOfVDWGridPoints.x;
        if(BlockingGrid[index][j][k]==PocketID)
        {
          Queue[QueueSize].x=index;
          Queue[QueueSize].y=j;
          Queue[QueueSize].z=k;
          QueueSize++;
        }
      }
    }
  }
  while(QueueSize>0);
  return MaxRadius;
}



int AnyRemainingPockets(int *seedx,int *seedy,int *seedz,int PocketID)
{
  int i,j,k;

  for(i=0;i<=NumberOfVDWGridPoints.x;i++)
    for(j=0;j<=NumberOfVDWGridPoints.y;j++)
      for(k=0;k<=NumberOfVDWGridPoints.z;k++)
        if(BlockingGrid[i][j][k]==PocketID)
        {
          *seedx=i;
          *seedy=j;
          *seedz=k;
          return TRUE;
        }
  return FALSE;
}

VECTOR *PocketsCenterOfMass;
REAL *PocketsRadius;

VECTOR DoubleCenters[1024];
int NumberOfDoubleCenters;
VECTOR StoredCenters[1024];
VECTOR FinalCenters[1024];
int CountFinalCenters[1024];
int NumberOfFinalCenters;

void ClassifyPockets(void)
{
  int i,j;
  int seedx,seedy,seedz;
  int PocketID;
  VECTOR center,dr,v1,v2,s;
  REAL radius,r;
  int CheckOverlap,ID;
  INT_VECTOR3 bl;

  seedx=seedy=seedz=0;

  PocketID=11;
  while(AnyRemainingPockets(&seedx,&seedy,&seedz,0))
  {
    printf("Pocket ID: %d\n",PocketID);
    FloodFillNonRecursive(seedx,seedy,seedz,PocketID++);
  };

  PocketsCenterOfMass=(VECTOR*)calloc(PocketID,sizeof(VECTOR));
  PocketsRadius=(REAL*)calloc(PocketID,sizeof(REAL));

  NumberOfDoubleCenters=0;
  for(i=11;i<PocketID;i++)
  {
    AnyRemainingPockets(&seedx,&seedy,&seedz,i);
    center=FloodFillNonRecursivePockets(seedx,seedy,seedz,i,50);

    if(CheckIfPointIsInUnitCell(center))
      DoubleCenters[NumberOfDoubleCenters++]=center;

  }

  NumberOfFinalCenters=0;
  printf("NumberOfDoubleCentersL %d\n",NumberOfDoubleCenters);
  for(i=0;i<NumberOfDoubleCenters;i++)
  {
    v1=DoubleCenters[i];

    CheckOverlap=FALSE;
    for(j=0;j<NumberOfFinalCenters;j++)
    {
      v2=StoredCenters[j];
      dr.x=v1.x-v2.x;
      dr.y=v1.y-v2.y;
      dr.z=v1.z-v2.z;
      dr=ApplyBoundaryConditionUnitCell(dr);
      r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
      if(r<2.0) 
      {
        CheckOverlap=TRUE;
        break;
      }
    }
    if(CheckOverlap==FALSE)
    {
       StoredCenters[NumberOfFinalCenters]=DoubleCenters[i];
       FinalCenters[NumberOfFinalCenters]=DoubleCenters[i];
       CountFinalCenters[NumberOfFinalCenters]++; 
       NumberOfFinalCenters++;
    }
    else 
    {
       FinalCenters[j].x+=DoubleCenters[j].x+dr.x;
       FinalCenters[j].y+=DoubleCenters[j].y+dr.y;
       FinalCenters[j].z+=DoubleCenters[j].z+dr.z;
       CountFinalCenters[j]++; 
    }
  }
  printf("NumberOfFinalCenters: %d\n",NumberOfFinalCenters);
  for(i=0;i<NumberOfFinalCenters;i++)
  {
    center.x=FinalCenters[i].x/CountFinalCenters[i];
    center.y=FinalCenters[i].y/CountFinalCenters[i];
    center.z=FinalCenters[i].z/CountFinalCenters[i];
/*
    s.x=InverseUnitCellBox[CurrentSystem].ax*center.x+InverseUnitCellBox[CurrentSystem].bx*center.y+InverseUnitCellBox[CurrentSystem].cx*center.z;
    s.y=InverseUnitCellBox[CurrentSystem].ay*center.x+InverseUnitCellBox[CurrentSystem].by*center.y+InverseUnitCellBox[CurrentSystem].cy*center.z;
    s.z=InverseUnitCellBox[CurrentSystem].az*center.x+InverseUnitCellBox[CurrentSystem].bz*center.y+InverseUnitCellBox[CurrentSystem].cz*center.z;
    printf("%g %g %g\n",s.x,s.y,s.z);
*/
    center=ConvertPositionToVTKPosition(center);
    printf("%g %g %g\n",center.x,center.y,center.z);
  }


  printf("\n\n");

  for(i=0;i<NumberOfFinalCenters;i++)
  {
    center.x=FinalCenters[i].x/CountFinalCenters[i];
    center.y=FinalCenters[i].y/CountFinalCenters[i];
    center.z=FinalCenters[i].z/CountFinalCenters[i];

    bl=ConvertXYZPositionToGridIndex(center);
    if((bl.x>=0)&&(bl.x<=NumberOfVDWGridPoints.x)&&
       (bl.y>=0)&&(bl.y<=NumberOfVDWGridPoints.y)&&
       (bl.z>=0)&&(bl.z<=NumberOfVDWGridPoints.z))
    {
    ID=BlockingGrid[bl.x][bl.y][bl.z];
    if(ID>2)
    {

    s.x=InverseUnitCellBox[CurrentSystem].ax*center.x+InverseUnitCellBox[CurrentSystem].bx*center.y+InverseUnitCellBox[CurrentSystem].cx*center.z;
    s.y=InverseUnitCellBox[CurrentSystem].ay*center.x+InverseUnitCellBox[CurrentSystem].by*center.y+InverseUnitCellBox[CurrentSystem].cy*center.z;
    s.z=InverseUnitCellBox[CurrentSystem].az*center.x+InverseUnitCellBox[CurrentSystem].bz*center.y+InverseUnitCellBox[CurrentSystem].cz*center.z;
    printf("%g %g %g ",s.x,s.y,s.z);

      radius=FloodFillNonRecursivePocketsRadius(bl.x,bl.y,bl.z,ID,center);
      printf("%d %g\n",i,radius);
    }
    
    }
  }
  exit(0);
}


REAL GetGridEnergy(int typeA,int x0,int y0,int z0)
{
  int i,j,k;
  REAL value;
  REAL X[64],a[64];

  // retrieve the grid value and the derivatives
  for(i=0;i<8;i++)
  {
    X[0+i*8]=VDWGrid[typeA][x0][y0][z0][i];
    X[1+i*8]=VDWGrid[typeA][x0+1][y0][z0][i];
    X[2+i*8]=VDWGrid[typeA][x0][y0+1][z0][i];
    X[3+i*8]=VDWGrid[typeA][x0+1][y0+1][z0][i];
    X[4+i*8]=VDWGrid[typeA][x0][y0][z0+1][i];
    X[5+i*8]=VDWGrid[typeA][x0+0][y0][z0+1][i];
    X[6+i*8]=VDWGrid[typeA][x0][y0+1][z0+1][i];
    X[7+i*8]=VDWGrid[typeA][x0+1][y0+1][z0+1][i];
  }

  for (i=0;i<64;i++)
  {
    a[i]=(double)(0.0);
    for (j=0;j<64;j++)
      a[i]+=Coeff[i][j]*X[j];
  }

  // energy
  value=0.0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      for (k=0;k<4;k++)
        value+=a[i+4*j+16*k]*pow(0.0,i)*pow(0.0,j)*pow(0.0,k);

  return value;
}

void BlockingVDWGrid(void)
{
  int i,j,k,l;
  REAL energy;

  // return if blocking the grids is not wanted
  if(!BlockEnergyGrids) return;

  // return if no grids are loaded
  if(!UseTabularGrid) return;


  // Allocate the memory for this grid
  BlockingGrid=(char***)calloc(NumberOfVDWGridPoints.x+1,sizeof(char**));
  for(i=0;i<=NumberOfVDWGridPoints.x;i++)
  {
    BlockingGrid[i]=(char**)calloc(NumberOfVDWGridPoints.y+1,sizeof(char*));
    for(j=0;j<=NumberOfVDWGridPoints.y;j++)
      BlockingGrid[i][j]=(char*)calloc(NumberOfVDWGridPoints.z+1,sizeof(char));
  }

  // initialize the queue for the flood-fill algorithm
  QueueSize=0;
  Queue=(INT_VECTOR3*)calloc((NumberOfVDWGridPoints.x+1)*(NumberOfVDWGridPoints.y+1)*(NumberOfVDWGridPoints.z+1),
                             sizeof(INT_VECTOR3));

  printf("BlockEnergyGridOverlapCriteria: %g\n",BlockEnergyGridOverlapCriteria);

  // patch grid: set blocking zones to 'overlap'
  for(l=0;l<NumberOfPseudoAtoms;l++)
  {
    if(VDWGrid[l])
    {
      // set energy-favorable regions to '0', overlap regions to '1'
      for(i=0;i<=NumberOfVDWGridPoints.x;i++)
        for(j=0;j<=NumberOfVDWGridPoints.y;j++)
          for(k=0;k<=NumberOfVDWGridPoints.z;k++)
          {
            energy=VDWGrid[l][i][j][k][0];
            if(energy<BlockEnergyGridOverlapCriteria)
              BlockingGrid[i][j][k]=0;
            else
              BlockingGrid[i][j][k]=1;
          }

       // before: low energy regions set to '0', high energy regions set to '1'
       // after: low energy regions of pores set to '2', high energy regions to '1',
       //        therefore the remaining values of '0' are considered as 'pockets'
       for(i=0;i<NumberOfGridSeeds;i++)
         FloodFillNonRecursive(GridSeeds[i].x,GridSeeds[i].y,GridSeeds[i].z,2);

       for(i=0;i<=NumberOfVDWGridPoints.x;i++)
         for(j=0;j<=NumberOfVDWGridPoints.y;j++)
           for(k=0;k<=NumberOfVDWGridPoints.z;k++)
           {
             // if the current voxel is a pocket and pockets are set to 'blocked' or
             // if the current voxel is a pore and pores are set to 'blocked', then block the voxel
             // (pores can also be blocked to visualize just the pockets in VTK)
             if(((BlockingGrid[i][j][k]==0)&&(BlockGridPockets))||((BlockingGrid[i][j][k]==2)&&(BlockGridPores)))
             {
               VDWGrid[l][i][j][k][0]=EnergyOverlapCriteria;
               VDWGrid[l][i][j][k][1]=EnergyOverlapCriteria;
               VDWGrid[l][i][j][k][2]=EnergyOverlapCriteria;
               VDWGrid[l][i][j][k][3]=EnergyOverlapCriteria;
               VDWGrid[l][i][j][k][4]=0.0;
               VDWGrid[l][i][j][k][5]=0.0;
               VDWGrid[l][i][j][k][6]=0.0;
               VDWGrid[l][i][j][k][7]=0.0;
             }
           }
    }
  }

  ClassifyPockets();   

  int count[200];
  int total_count=0;

  for(l=0;l<50;l++)
  {
    count[l]=0;
    for(i=0;i<=NumberOfVDWGridPoints.x;i++)
      for(j=0;j<=NumberOfVDWGridPoints.y;j++)
        for(k=0;k<=NumberOfVDWGridPoints.z;k++)
          count[(int)BlockingGrid[i][j][k]]++;
    total_count+=count[l];
    if(count[l]>0) printf("count: %d -> %d\n",l,count[l]);
  }
  printf("total: %d %d\n",(1+NumberOfVDWGridPoints.x)*(1+NumberOfVDWGridPoints.y)*(1+NumberOfVDWGridPoints.z),total_count);


  printf("exit\n");

  free(Queue);
}

void WriteRestartGrids(FILE *FilePtr)
{
  fwrite(&SpacingVDWGrid,1,sizeof(REAL),FilePtr);
  fwrite(&SpacingCoulombGrid,1,sizeof(REAL),FilePtr);
  fwrite(&NumberOfGrids,1,sizeof(int),FilePtr);
  fwrite(GridTypeList,NumberOfGrids,sizeof(int),FilePtr);
}

void AllocateGridMemory(void)
{

  GridTypeList=(int*)calloc(NumberOfGrids,sizeof(int));
  VDWGrid=(float*****)calloc(NumberOfPseudoAtoms,sizeof(float****));
  DelaunayGrid=(int*****)calloc(NumberOfPseudoAtoms,sizeof(int****));
}

void ReadRestartGrids(FILE *FilePtr)
{
  fread(&SpacingVDWGrid,1,sizeof(REAL),FilePtr);
  fread(&SpacingCoulombGrid,1,sizeof(REAL),FilePtr);
  fread(&NumberOfGrids,1,sizeof(int),FilePtr);

  AllocateGridMemory();
  fread(GridTypeList,NumberOfGrids,sizeof(int),FilePtr);

  if(Framework[0].FrameworkModel==GRID)
  {
    CurrentSystem=0;
    ReadVDWGrid();
    if(ChargeMethod!=NONE)
      ReadCoulombGrid();
  }
}

int WriteDelaunayGrid(int l)
{
	int i,j,k,m;
	int ngrid;
	FILE *FilePtr;
	char buffer[256];
	
	ngrid=NumberOfDelaunayGridPoints.x*NumberOfDelaunayGridPoints.y*NumberOfDelaunayGridPoints.z;
	
	sprintf(buffer,"%s/share/raspa/grids",RASPA_DIRECTORY);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s",RASPA_DIRECTORY,ForceField);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s/%s",RASPA_DIRECTORY,ForceField,Framework[CurrentSystem].Name[0]);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf",RASPA_DIRECTORY,ForceField,
			Framework[CurrentSystem].Name[0],(double)SpacingDelaunayGrid);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%s_%s_Delaunay.grid",
			RASPA_DIRECTORY,
			ForceField,
			Framework[CurrentSystem].Name[0],
			(double)SpacingDelaunayGrid,
			Framework[CurrentSystem].Name[0],
			PseudoAtoms[l].Name);
	FilePtr=fopen(buffer,"w");
	
	fwrite(&SpacingDelaunayGrid,1,sizeof(REAL),FilePtr);
	fwrite(&NumberOfDelaunayGridPoints,1,sizeof(INT_VECTOR3),FilePtr);
	fwrite(&SizeGrid,1,sizeof(VECTOR),FilePtr);
	fwrite(&ShiftGrid,1,sizeof(VECTOR),FilePtr);
	fwrite(&DeltaDelaunayGrid,1,sizeof(VECTOR),FilePtr);
	fwrite(&UnitCellSize[CurrentSystem],1,sizeof(VECTOR),FilePtr);
	fwrite(&NumberOfUnitCells[CurrentSystem],1,sizeof(INT_VECTOR3),FilePtr);
	
	for(m=0;m<8;m++)
		for(i=0;i<=NumberOfDelaunayGridPoints.x;i++)
			for(j=0;j<=NumberOfDelaunayGridPoints.y;j++)
				for(k=0;k<=NumberOfDelaunayGridPoints.z;k++)
					fwrite(&DelaunayGrid[l][i][j][k][m],1,sizeof(int),FilePtr);
	
	
	fclose(FilePtr);
	return 0;
}

void ReadDelaunayGrid(void)
{
	int i,j,k,l,m;
	VECTOR unit_cell_size;
	INT_VECTOR3 number_of_unit_cells;
	FILE *FilePtr;
	char buffer[256];
	
	printf("Reading VDW grid\n");
	
	for(i=0;i<NumberOfPseudoAtoms;i++)
		DelaunayGrid[i]=NULL;
	
	for(l=0;l<NumberOfGrids;l++)
	{
		sprintf(buffer,"%s/share/raspa/grids",RASPA_DIRECTORY);
		sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%s_%s_Delaunay.grid",
				RASPA_DIRECTORY,
				ForceField,
				Framework[CurrentSystem].Name[0],
				(double)SpacingVDWGrid,
				Framework[CurrentSystem].Name[0],
				PseudoAtoms[GridTypeList[l]].Name);
		printf("Opening: %s\n",buffer);
		if(!(FilePtr=fopen(buffer,"r")))
		{
			printf("Error:  file %s does not exists.\n",buffer);
			exit(1);
		}
		
		fread(&SpacingDelaunayGrid,1,sizeof(REAL),FilePtr);
		fread(&NumberOfDelaunayGridPoints,1,sizeof(INT_VECTOR3),FilePtr);
		fread(&SizeGrid,1,sizeof(VECTOR),FilePtr);
		fread(&ShiftGrid,1,sizeof(VECTOR),FilePtr);
		fread(&DeltaDelaunayGrid,1,sizeof(VECTOR),FilePtr);
		fread(&unit_cell_size,1,sizeof(VECTOR),FilePtr);
		fread(&number_of_unit_cells,1,sizeof(INT_VECTOR3),FilePtr);
		
		if((fabs(UnitCellSize[CurrentSystem].x-unit_cell_size.x)>1e-8)||(fabs(UnitCellSize[CurrentSystem].y-unit_cell_size.y)>1e-8)||
		   (fabs(UnitCellSize[CurrentSystem].z-unit_cell_size.z)>1e-8))
		{
			printf("Error in reading the grid: unit-cell size does not match\n");
			exit(0);
		}
		
		// Allocate the memory for this grid
		DelaunayGrid[GridTypeList[l]]=(int****)calloc(NumberOfDelaunayGridPoints.x+1,sizeof(int***));
		for(i=0;i<=NumberOfDelaunayGridPoints.x;i++)
		{
			DelaunayGrid[GridTypeList[l]][i]=(int***)calloc(NumberOfDelaunayGridPoints.y+1,sizeof(int**));
			for(j=0;j<=NumberOfDelaunayGridPoints.y;j++)
			{
				DelaunayGrid[GridTypeList[l]][i][j]=(int**)calloc(NumberOfDelaunayGridPoints.z+1,sizeof(int*));
				for(k=0;k<=NumberOfDelaunayGridPoints.z;k++)
					DelaunayGrid[GridTypeList[l]][i][j][k]=(int*)calloc(1,sizeof(int));
			}
		}
		
		for(m=0;m<8;m++)
			for(i=0;i<=NumberOfDelaunayGridPoints.x;i++)
				for(j=0;j<=NumberOfDelaunayGridPoints.y;j++)
					for(k=0;k<=NumberOfDelaunayGridPoints.z;k++)
						fread(&DelaunayGrid[GridTypeList[l]][i][j][k][m],1,sizeof(float),FilePtr);
		
		
		fclose(FilePtr);
	}
}

int PointInPolygon(VECTOR pos,VECTOR Ray, int i,int j)
{
	int LinePlanePoint, k, l, c = 0;
	VECTOR TestPoint,TP;
	
	//	printf("PointOnPlane %f %f %f \n",PointOnPlane.x,PointOnPlane.y,PointOnPlane.z);
	LinePlanePoint=LinePlaneIntersect(i, j, pos, Ray);
	//	printf("PointOnPlane %f %f %f LinePlane Point %i pos %f %f %f Ray %f %f %f\n",PointOnPlane.x,PointOnPlane.y,PointOnPlane.z,LinePlanePoint,pos.x,pos.y,pos.z,Ray.x,Ray.y,Ray.z);
	
	switch (LinePlanePoint) {
		case 0:
			return 0;
			break;
		case 1:
			//Convert PointOnPlane to 2D Coordinates with plane of interest
			TP=PointOnPlane;
			TestPoint=CoordinateSystemShift(TP, i, j);
			
			for (k = 0, l = PolyNumVertPerFace[i][j]-1; k < PolyNumVertPerFace[i][j]; l = k++) {	
				if ( ((Poly2DPoints[i][j][k].y>TestPoint.y) != (Poly2DPoints[i][j][l].y>TestPoint.y)) &&
					(TestPoint.x < (Poly2DPoints[i][j][l].x-Poly2DPoints[i][j][k].x) * (TestPoint.y-Poly2DPoints[i][j][k].y) / (Poly2DPoints[i][j][l].y-Poly2DPoints[i][j][k].y) + Poly2DPoints[i][j][k].x) )
					c = !c;
			}
			return c;
			break;
		case 2:
			//Convert PointOnPlane to 2D Coordinates with plane of interest
			TP=PointOnPlane;
			TestPoint=CoordinateSystemShift(TP, i, j);
			
			for (k = 0, l = PolyNumVertPerFace[i][j]-1; k < PolyNumVertPerFace[i][j]; l = k++) {	
				if ( ((Poly2DPoints[i][j][k].y>TestPoint.y) != (Poly2DPoints[i][j][l].y>TestPoint.y)) &&
					(TestPoint.x < (Poly2DPoints[i][j][l].x-Poly2DPoints[i][j][k].x) * (TestPoint.y-Poly2DPoints[i][j][k].y) / (Poly2DPoints[i][j][l].y-Poly2DPoints[i][j][k].y) + Poly2DPoints[i][j][k].x) )
					c = !c;
			}
			if (c==1) {
				c=10;
			}
			return c;
			break;
		default:
			break;
	}
	
	//Convert PointOnPlane to 2D Coordinates with plane of interest
	TP=PointOnPlane;
	TestPoint=CoordinateSystemShift(TP, i, j);
	
	for (k = 0, l = PolyNumVertPerFace[i][j]-1; k < PolyNumVertPerFace[i][j]; l = k++) {	
		if ( ((Poly2DPoints[i][j][k].y>TestPoint.y) != (Poly2DPoints[i][j][l].y>TestPoint.y)) &&
			(TestPoint.x < (Poly2DPoints[i][j][l].x-Poly2DPoints[i][j][k].x) * (TestPoint.y-Poly2DPoints[i][j][k].y) / (Poly2DPoints[i][j][l].y-Poly2DPoints[i][j][k].y) + Poly2DPoints[i][j][k].x) )
			c = !c;
	}
	if (LinePlanePoint==2 && c==1) {
		c=10;
	}
	return c;
}

VECTOR PointProject(VECTOR PointToProject,int i, int j)
{
	float alpha = 0;
	VECTOR projectedpoint;
	
	//	PolyPlanes[i][j][0]=A
	//	PolyPlanes[i][j][1]=B
	//	PolyPlanes[i][j][2]=C
	//	PolyPlanes[i][j][3]=D
	
	alpha = (PolyPlanes[i][j][0]*PointToProject.x+PolyPlanes[i][j][1]*PointToProject.y+PolyPlanes[i][j][2]*PointToProject.z+PolyPlanes[i][j][3])/(PolyPlanes[i][j][0]*PolyPlanes[i][j][0]+PolyPlanes[i][j][1]*PolyPlanes[i][j][1]+PolyPlanes[i][j][2]*PolyPlanes[i][j][2]);
	
	projectedpoint.x=PointToProject.x-PolyPlanes[i][j][0]*alpha;
	projectedpoint.y=PointToProject.y-PolyPlanes[i][j][1]*alpha;
	projectedpoint.z=PointToProject.z-PolyPlanes[i][j][2]*alpha;
	
	return projectedpoint;
	
}

VECTOR CoordinateSystemShift(VECTOR PointToShift, int i, int j)
{
	
	VECTOR TempPoint;
	
	TempPoint.x = PointToShift.x;
	TempPoint.y = PointToShift.y;
	TempPoint.z = PointToShift.z;
	//Correct for when plane lies in YZ coordinates only (B and C are 0)
	if (PolyPlanes[i][j][1]>-0.000001 && PolyPlanes[i][j][1]<0.000001 && PolyPlanes[i][j][2]>-0.000001 && PolyPlanes[i][j][2]<0.000001) {
		PointToShift.z = DotProduct(TempPoint, Plane2D[i][j][0]);
		PointToShift.y = DotProduct(TempPoint, Plane2D[i][j][1]);
		PointToShift.x = DotProduct(TempPoint, Plane2D[i][j][2]);
	}
	
	else {
		PointToShift.x = DotProduct(TempPoint, Plane2D[i][j][0]);
		PointToShift.y = DotProduct(TempPoint, Plane2D[i][j][1]);
		PointToShift.z = DotProduct(TempPoint, Plane2D[i][j][2]);
	}
	
	//	printf("Temp point = %f %f %f PointToShift = %f %f %f\n",TempPoint.x,TempPoint.y,TempPoint.z,PointToShift.x,PointToShift.y,PointToShift.z);
	return PointToShift;
}

void PlaneFinder(int i, int j)
{
	float A=0,B=0,C=0,D=0;
	
	VECTOR Xnew, Ynew, Znew, UnitXnew, UnitYnew, UnitZnew, Xold;
	
	Xold.x = 1;
	Xold.y = 0;
	Xold.z = 0;
	
	//  Use three Vertices to find plane values A,B,C,D.
	A = PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].y*(PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].z-PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].z)+PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].y*(PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].z-PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].z)+PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].y*(PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].z-PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].z);
	B = PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].z*(PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].x-PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].x)+PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].z*(PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].x-PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].x)+PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].z*(PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].x-PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].x);
	C = PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].x*(PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].y-PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].y)+PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].x*(PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].y-PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].y)+PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].x*(PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].y-PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].y);
	D = -1*(PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].x*(PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].y*PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].z-PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].y*PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].z)+PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].x*(PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].y*PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].z-PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].y*PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].z)+PolyVertLocation[i][PolyVertPerFace[i][j][2]-1].x*(PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].y*PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].z-PolyVertLocation[i][PolyVertPerFace[i][j][1]-1].y*PolyVertLocation[i][PolyVertPerFace[i][j][0]-1].z));
	
	PolyPlanes[i][j][0]=A;
	PolyPlanes[i][j][1]=B;
	PolyPlanes[i][j][2]=C;
	PolyPlanes[i][j][3]=D;
	
	//printf("Plane = %f %f %f %f i = %d j= %d\n",PolyPlanes[i][j][0],PolyPlanes[i][j][1],PolyPlanes[i][j][2],PolyPlanes[i][j][3],i,j);
	//Correct for when plane lies in YZ coordinates only (B and C are 0)
	if (PolyPlanes[i][j][1]>-0.000001 && PolyPlanes[i][j][1]<0.000001 && PolyPlanes[i][j][2]>-0.000001 && PolyPlanes[i][j][2]<0.000001) {
		Znew.z = A;
		Znew.y = B;
		Znew.x = C;
	}
	else {
		Znew.x = A;
		Znew.y = B;
		Znew.z = C;
	}
	
	Ynew = CrossProduct(Xold, Znew);
	Xnew = CrossProduct(Znew, Ynew);
	
	UnitXnew = NormalizeVector(Xnew);
	UnitYnew = NormalizeVector(Ynew);
	UnitZnew = NormalizeVector(Znew);
	
	NormalizeVector(Xnew);
	NormalizeVector(Ynew);
	NormalizeVector(Znew);
	
	Plane2D[i][j][0] = UnitXnew;
	Plane2D[i][j][1] = UnitYnew;
	Plane2D[i][j][2] = UnitZnew; 
	
}	

void ReadPolyhedraData(void)
{
	
	char buffer[256];
	char str [30];
	float DelaunayUnitCellSizeX,DelaunayUnitCellSizeY,DelaunayUnitCellSizeZ;
	int i,j,k,l,m,c,lines,TempEdgeVert1,TempEdgeVert2,DuplicateEdge;
	FILE *pFile;
	lines=0;
	
	sprintf(buffer,"%s/share/raspa/grids/Delaunay/%s",RASPA_DIRECTORY,DelaunayPolyhedraFileName);
	printf("Opening: %s\n",buffer);
	
	pFile= fopen (buffer,"r");
	if (pFile!=NULL)
	{
		while ((c=fgetc(pFile)) != EOF) {
			if (c== '\n') {
				++lines;
			}
		}
		
		rewind(pFile);
		
		for (c=0; c<lines; c++) {
			fgets(str, 30, pFile);
			sscanf(str, "%d",&NumberPolyhedra);
		}
		
		rewind(pFile);
		
		fgets(str, 30, pFile);
		sscanf(str, "%f %f %f",&DelaunayUnitCellSizeX,&DelaunayUnitCellSizeY,&DelaunayUnitCellSizeZ);
		
		PolyNumberFaces=(int*)calloc(NumberPolyhedra,sizeof(int));
		PolyNumberVertices=(int*)calloc(NumberPolyhedra,sizeof(int));
		PolyNumberEdges=(int*)calloc(NumberPolyhedra,sizeof(int));
		PolyEdges=(int***)calloc(NumberPolyhedra,sizeof(int**));
		PolyNumVertPerFace=(int**)calloc(NumberPolyhedra,sizeof(int*));
		PolyVertices=(int**)calloc(NumberPolyhedra,sizeof(int*));
		PolyVertPerFace=(int***)calloc(NumberPolyhedra,sizeof(int**));
		PolyDelsVertices=(INT_VECTOR3**)calloc(NumberPolyhedra,sizeof(INT_VECTOR3*));
		PolyPlanes=(float***)calloc(NumberPolyhedra, sizeof(float**));
		Plane2D=(VECTOR***)calloc(NumberPolyhedra, sizeof(VECTOR**));
		Poly2DPoints=(VECTOR***)calloc(NumberPolyhedra,sizeof(VECTOR**));
		PolyVertLocation=(VECTOR**)calloc(NumberPolyhedra,sizeof(VECTOR*));
		
		for (i=0; i<NumberPolyhedra; i++) {
			fgets(str, 30, pFile);
			sscanf(str, "%d",&PolyNumberVertices[i]);
			fgets(str, 30, pFile);
			sscanf(str, "%d",&PolyNumberFaces[i]);
			
			PolyNumVertPerFace[i]=(int*)calloc(PolyNumberFaces[i],sizeof(int));
			PolyVertPerFace[i]=(int**)calloc(PolyNumberFaces[i],sizeof(int*));
			PolyEdges[i]=(int**)calloc((PolyNumberFaces[i]+PolyNumberVertices[i]-2),sizeof(int*));
			PolyPlanes[i]=(float**)calloc(PolyNumberFaces[i], sizeof(float*));
			Poly2DPoints[i]=(VECTOR**)calloc(PolyNumberFaces[i], sizeof(VECTOR*));
			Plane2D[i]=(VECTOR**)calloc(PolyNumberFaces[i], sizeof(VECTOR*));
			PolyDelsVertices[i]=(INT_VECTOR3*)calloc(PolyNumberVertices[i],sizeof(INT_VECTOR3));
			PolyVertices[i]=(int*)calloc(PolyNumberVertices[i],sizeof(int));
			PolyVertLocation[i]=(VECTOR*)calloc(PolyNumberVertices[i],sizeof(VECTOR));
			
			for (j=0; j<(PolyNumberFaces[i]+PolyNumberVertices[i]-2); j++) {
				PolyEdges[i][j]=(int*)calloc(2,sizeof(int));
			}
			
			for (j=0; j<PolyNumberFaces[i]; j++) {
				fgets(str, 30, pFile);
				sscanf(str, "%d",&PolyNumVertPerFace[i][j]);
				PolyVertPerFace[i][j]=(int*)calloc(PolyNumVertPerFace[i][j],sizeof(int));
				PolyPlanes[i][j]=(float*)calloc(4,sizeof(float));
				Poly2DPoints[i][j]=(VECTOR*)calloc(PolyNumVertPerFace[i][j], sizeof(VECTOR));
				Plane2D[i][j]=(VECTOR*)calloc(3, sizeof(VECTOR));
			}
			
			PolyNumberEdges[i] = 0;
			for (j=0; j<PolyNumberFaces[i]; j++) {
				for (k=0; k<PolyNumVertPerFace[i][j]; k++) {
					fgets(str, 30, pFile);
					sscanf(str, "%d",&PolyVertPerFace[i][j][k]);
					
				}
				
				for (k=0, l=PolyNumVertPerFace[i][j]-1; k<PolyNumVertPerFace[i][j]; l=k++) {
					TempEdgeVert1 = PolyVertPerFace[i][j][k];
					TempEdgeVert2 = PolyVertPerFace[i][j][l];
					
					DuplicateEdge = 0;
					m=0;
					do {
						if ((PolyEdges[i][m][0] == TempEdgeVert1 && PolyEdges[i][m][1]== TempEdgeVert2) || (PolyEdges[i][m][1] == TempEdgeVert1 && PolyEdges[i][m][0]== TempEdgeVert2))  {
							DuplicateEdge = 1;
						}
						m++;
					}while (m<PolyNumberEdges[i] && DuplicateEdge == 0 );
					switch (DuplicateEdge) {
						case 0:
							PolyEdges[i][PolyNumberEdges[i]][0] = TempEdgeVert1;
							PolyEdges[i][PolyNumberEdges[i]][1] = TempEdgeVert2;
							//							printf("Edge %i of poly %i is %i %i\n",PolyNumberEdges[i],i,PolyEdges[i][PolyNumberEdges[i]][0],PolyEdges[i][PolyNumberEdges[i]][1]);
							PolyNumberEdges[i]++;
							break;
						default:
							break;
					}
				}	
				
			}
			//			printf("Number of edges = %i, should be = %i\n",PolyNumberEdges[i],(PolyNumberFaces[i]+PolyNumberVertices[i]-2));
			
			for (j=0; j<PolyNumberVertices[i]; j++) {
				fgets(str, 30, pFile);
				sscanf(str, "%d",&PolyVertices[i][j]);
				
			}
			
			
			for (j=0; j<PolyNumberVertices[i]; j++) {
				fgets(str, 30, pFile);
				sscanf(str,"%d",&PolyDelsVertices[i][j].x);
				fgets(str, 30, pFile);
				sscanf(str,"%d",&PolyDelsVertices[i][j].y);
				fgets(str, 30, pFile);
				sscanf(str,"%d",&PolyDelsVertices[i][j].z);
				
			}
		}
		
		
	}
	else {
		printf("Error:  file  does not exists.\n");
		exit(1);
	}
	
	
	
	
}

VECTOR PBCCorrection(VECTOR pos,int i, int j)
{
	
	VECTOR s;
	
	switch(BoundaryCondition[CurrentSystem])
	{
		case RECTANGULAR:
			// the position has to be moved back to the main unit cell using the rectangular boundary condition
			pos.x-=UnitCellSize[CurrentSystem].x*(REAL)(NINT(pos.x/UnitCellSize[CurrentSystem].x));
			pos.y-=UnitCellSize[CurrentSystem].y*(REAL)(NINT(pos.y/UnitCellSize[CurrentSystem].y));
			pos.z-=UnitCellSize[CurrentSystem].z*(REAL)(NINT(pos.z/UnitCellSize[CurrentSystem].z));
			if(pos.x<0.0) pos.x+=UnitCellSize[CurrentSystem].x;
			if(pos.y<0.0) pos.y+=UnitCellSize[CurrentSystem].y;
			if(pos.z<0.0) pos.z+=UnitCellSize[CurrentSystem].z;
			
			pos.x+=UnitCellSize[CurrentSystem].x*PolyDelsVertices[i][j].x;
			pos.y+=UnitCellSize[CurrentSystem].y*PolyDelsVertices[i][j].y;
			pos.z+=UnitCellSize[CurrentSystem].z*PolyDelsVertices[i][j].z;
			
			break;
		case TRICLINIC:
		default:
			// the position first has to be moved back to the main unit cell using the triclinic boundary condition
			s.x=InverseUnitCellBox[CurrentSystem].ax*pos.x+InverseUnitCellBox[CurrentSystem].bx*pos.y+InverseUnitCellBox[CurrentSystem].cx*pos.z;
			s.y=InverseUnitCellBox[CurrentSystem].ay*pos.x+InverseUnitCellBox[CurrentSystem].by*pos.y+InverseUnitCellBox[CurrentSystem].cy*pos.z;
			s.z=InverseUnitCellBox[CurrentSystem].az*pos.x+InverseUnitCellBox[CurrentSystem].bz*pos.y+InverseUnitCellBox[CurrentSystem].cz*pos.z;
			
			s.x-=(REAL)NINT(s.x);
			s.y-=(REAL)NINT(s.y);
			s.z-=(REAL)NINT(s.z);
			
			if(s.x<0.0) s.x+=1.0;
			if(s.y<0.0) s.y+=1.0;
			if(s.z<0.0) s.z+=1.0;
			
			pos.x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
			pos.y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
			pos.z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
			
			pos.x+=UnitCellBox[CurrentSystem].ax*PolyDelsVertices[i][j].x+UnitCellBox[CurrentSystem].bx*PolyDelsVertices[i][j].y+UnitCellBox[CurrentSystem].cx*PolyDelsVertices[i][j].z;
			pos.y+=UnitCellBox[CurrentSystem].ay*PolyDelsVertices[i][j].x+UnitCellBox[CurrentSystem].by*PolyDelsVertices[i][j].y+UnitCellBox[CurrentSystem].cy*PolyDelsVertices[i][j].z;
			pos.z+=UnitCellBox[CurrentSystem].az*PolyDelsVertices[i][j].x+UnitCellBox[CurrentSystem].bz*PolyDelsVertices[i][j].y+UnitCellBox[CurrentSystem].cz*PolyDelsVertices[i][j].z;
			
			break;
	}
	
	pos.x=(pos.x-ShiftGrid.x);
	pos.y=(pos.y-ShiftGrid.y);
	pos.z=(pos.z-ShiftGrid.z);
	
	//	printf("i = %i, j = %i, PBCFrameworkpos %f %f %f\n",i,j,pos.x,pos.y,pos.z);
	
	return pos;
	
}

VECTOR RandomRay(void)
{
	VECTOR Ray;
	REAL randnum;
	
	//Need to correlate 200 to some number in the unit cell...
	randnum = (RandomNumber()-0.5);
	if (randnum<=0) {
		randnum = -1.0;
	}
	else if (randnum>0){
		randnum = 1.0;
	}
	
	Ray.x=(RandomNumber()+1)*200*randnum;
	Ray.y=(RandomNumber()+1)*200*randnum;
	Ray.z=(RandomNumber()+1)*200*randnum;
	
	return Ray;
	
}

int CollinearVertex(VECTOR Vertex,VECTOR Ray, VECTOR pos)
{
	
	VECTOR aminb,bminc,cross;
	double svalue = 0.0000001;
	
	aminb.x = pos.x - Ray.x;
	aminb.y = pos.y - Ray.y;
	aminb.z = pos.z - Ray.z;
	bminc.x = Ray.x - Vertex.x;
	bminc.y = Ray.y - Vertex.y;
	bminc.z = Ray.z - Vertex.z;
	cross = CrossProduct(aminb, bminc);
	
	if (cross.x <= svalue && cross.x >= -svalue && cross.y <= svalue && cross.y >= -svalue && cross.z <= svalue && cross.z >= -svalue) {
		return 1;
	}
	else {
		return 0;
	}
	
}

int EdgeIntersection(VECTOR VertexOne, VECTOR VertexTwo,VECTOR Ray, VECTOR pos)
{
	REAL r,svalue;
	VECTOR P1minP3,P2minP1,P4minP3,a,b,d,cross,t;
	svalue = 0.000001;
	
	//Check if skew
	P1minP3.x = VertexOne.x - Ray.x;
	P1minP3.y = VertexOne.y - Ray.y;
	P1minP3.z = VertexOne.z - Ray.z;
	P2minP1.x = VertexTwo.x - VertexOne.x;
	P2minP1.y = VertexTwo.y - VertexOne.y;
	P2minP1.z = VertexTwo.z - VertexOne.z;
	P4minP3.x = pos.x - Ray.x;
	P4minP3.y = pos.y - Ray.y;
	P4minP3.z = pos.z - Ray.z;
	d = CrossProduct(P2minP1, P4minP3);
	r = DotProduct(P1minP3, d);
	
	if (r >= svalue || r<= -svalue) {
		return 0;
	}
	//check if on edge
	else {
		a.x=VertexOne.x-VertexTwo.x;
		a.y=VertexOne.y-VertexTwo.y;
		a.z=VertexOne.z-VertexTwo.z;
		b.x=VertexTwo.x - pos.x;
		b.y=VertexTwo.y - pos.y;
		b.z=VertexTwo.z - pos.z;
		cross= CrossProduct(a, b);
		t.x = (pos.x-VertexOne.x)/(P2minP1.x);
		if (cross.x <= svalue && cross.x >= -svalue && cross.y <= svalue && cross.y >= -svalue && cross.z <= svalue && cross.z >= -svalue) {
			//Check if between two vertices
			if (t.x<1 && t.x>0) {
				return 2;
			}
			else {
				return 3;
			}
		}
		//Otherwise edge intersection
		else {
			return 1;
		}
	}
}

int LinePlaneIntersect(int i,int j, VECTOR pos, VECTOR Ray)
{
	VECTOR line;
	REAL Parallel, OnPlane,u;
	
	line.x = Ray.x - pos.x;
	line.y = Ray.y - pos.y;
	line.z = Ray.z - pos.z;
	Parallel = PolyPlanes[i][j][0]*line.x + PolyPlanes[i][j][1]*line.y + PolyPlanes[i][j][2]*line.z;
	OnPlane = PolyPlanes[i][j][0]*pos.x + PolyPlanes[i][j][1]*pos.y + PolyPlanes[i][j][2]*pos.z + PolyPlanes[i][j][3];
	
	if (OnPlane  >= -0.0000001 && OnPlane <= 0.0000001) {
		//Point is on the Plane
		PointOnPlane = pos;
		return 2;
	}
	//if Parallel is 0, then ray is parallel
	else if (Parallel >= -0.0000001 && Parallel<= 0.0000001) {
		printf("Point %f %f %f is parallel to plane %i\n",pos.x,pos.y,pos.z,j);		
		return 0;
		
	}
	else {
		u=(-OnPlane)/(Parallel);
		if (u<=1 && u>=0) {
			PointOnPlane.x = pos.x +u*(line.x);
			PointOnPlane.y = pos.y +u*(line.y);
			PointOnPlane.z = pos.z +u*(line.z);
			return 1;
		}
		//Ray does not intersect
		else {
			return 0;
		}
		
	}
	
}

int PointInPolyhedra(int i,int Inaccessible,VECTOR pos)
{
	VECTOR TestRay;
	int NotUsable,counter,sumcounter,j;
	NotUsable = 1;
	while (NotUsable > 0.5 && Inaccessible==0) {
		
		TestRay=RandomRay();
		//Test for Collinear points with vertex
		for (j=0; j<PolyNumberVertices[i]; j++) {
			NotUsable=CollinearVertex(PolyVertLocation[i][j], TestRay, pos);
			if (NotUsable>0.5) {
				printf("Found Collinear Vertex = %f %f %f Ray = %f %f %f pos = %f %f %f\n",PolyVertLocation[i][j].x,PolyVertLocation[i][j].y,PolyVertLocation[i][j].z,TestRay.x,TestRay.y,TestRay.z,pos.x,pos.y,pos.z);
				fflush(stdout);
				break;
			}
		}
		
		j=0;
		while (NotUsable <0.5 && j<PolyNumberEdges[i] && Inaccessible ==0) {
			NotUsable=EdgeIntersection(PolyVertLocation[i][PolyEdges[i][j][0]-1], PolyVertLocation[i][PolyEdges[i][j][1]-1], TestRay, pos);
			if (NotUsable>0.5) {
				switch (NotUsable) {
					case 1:
						//Ray intersects edge, but point is not collinear nor on edge
						break;
					case 2:
						//Point is collinear with edge and is on the edge
						Inaccessible=1;
						break;
					case 3:
						//Point is collinear with edge but not on edge special case
						Inaccessible=3;
						break;
					default:
						break;
				}
			}
			if (Inaccessible >0) {
				break;
			}
			j++;
		}
	}
	
	if (Inaccessible==0) {
		counter=0;
		sumcounter=0;
		for (j=0; j<PolyNumberFaces[i]; j++) {
			counter=PointInPolygon(pos, TestRay, i, j);
			//Check if point is on plane and in polygon
			if (counter == 10 && PointOnPlane.x == pos.x && PointOnPlane.y == pos.y && PointOnPlane.z == pos.z) {
				sumcounter = 1;
				break;
			}
			sumcounter += counter;
		}
		//check if counter is odd or even
		if (sumcounter % 2) {
			Inaccessible = 1;
		}
		else {
			Inaccessible = 0;
		}
		
	}
	//special case if point is collinear with an edge but not on edge change inaccessible back to false
	else if (Inaccessible == 3)
	{
		Inaccessible =0;
	}
	
	
	return Inaccessible;
}

void CreateEBCBMCProbBins(void)
{
	
	int i,j,k,l,typeT,t;
	REAL TempEnergy[(NumberOfVDWGridPoints.x+1)*(NumberOfVDWGridPoints.y+1)*(NumberOfVDWGridPoints.z+1)];
	REAL Energy,Energy2;
	INT_VECTOR3 GridLocation[(NumberOfVDWGridPoints.x)*(NumberOfVDWGridPoints.y)*(NumberOfVDWGridPoints.z)];
	
	k=0;
	
	printf("Start EBCBMC Allocation, NumberOfGrids = %i\n",NumberOfGrids);
	//Allocate Memory for EBCBMCValues
	EBCBMCValue=(float***)calloc(NumberOfPseudoAtoms,sizeof(float**));
	EBCBMCGridLocation=(INT_VECTOR3***)calloc(NumberOfPseudoAtoms,sizeof(INT_VECTOR3**));
	SumVDWGridRosenbluth=(REAL*)calloc(NumberOfPseudoAtoms,sizeof(REAL));
	SizeOfEBCBMCBins=(int**)calloc(NumberOfPseudoAtoms,sizeof(int*));
	
	EBCBMCVolume=NumberOfVDWGridPoints.x*NumberOfVDWGridPoints.y*NumberOfVDWGridPoints.z;
	
	Energy=0;
	Energy2=0;
	for (l=0; l<NumberOfGrids; l++) {
		SizeOfEBCBMCBins[GridTypeList[l]]=(int*)calloc(NumberOfEBCBMCBins,sizeof(int));
		EBCBMCValue[GridTypeList[l]]=(float**)calloc(NumberOfEBCBMCBins,sizeof(float*));
		EBCBMCGridLocation[GridTypeList[l]]=(INT_VECTOR3**)calloc(NumberOfEBCBMCBins,sizeof(INT_VECTOR3*));
	}
	
	for(l=0;l<NumberOfGrids;l++){
		t=0;
		typeT=GridTypeList[l];
		SumVDWGridRosenbluth[typeT]=0;
		Energy2=0.0;
		for (i=0;i<NumberOfVDWGridPoints.x;i++) {
			for (j=0; j<NumberOfVDWGridPoints.y; j++) {
				for (k=0; k<NumberOfVDWGridPoints.z; k++) {
					Energy=VDWGrid[typeT][i][j][k][0];
					if (Energy<Energy2) {
						Energy2=Energy;
					}
					if (Energy<(EnergyOverlapCriteria)) {
						SumVDWGridRosenbluth[typeT]+=Exp(-Beta[CurrentSystem]*Energy);
						TempEnergy[t]=SumVDWGridRosenbluth[typeT];
						GridLocation[t].x=i;
						GridLocation[t].y=j;
						GridLocation[t].z=k;
						t++;
					}
				}
			}
		}
		printf("Lowest Energy = %f\n",Energy2);
		for (i=0; i<NumberOfEBCBMCBins; i++) {
			SizeOfEBCBMCBins[typeT][i]=0;
		}
		for (i=0; i<t; i++) {
			TempEnergy[i]/=SumVDWGridRosenbluth[typeT];
			k=(int)(TempEnergy[i]*NumberOfEBCBMCBins);
			if (k==NumberOfEBCBMCBins) {
				SizeOfEBCBMCBins[typeT][k-1]++;
			}
			else {
				SizeOfEBCBMCBins[typeT][k]++;
			}
		}
		//Make Grids for last part
		
		for (i=0; i<NumberOfEBCBMCBins; i++) {
			EBCBMCValue[typeT][i]=(float*)calloc(SizeOfEBCBMCBins[typeT][i],sizeof(float));
			EBCBMCGridLocation[typeT][i]=(INT_VECTOR3*)calloc(SizeOfEBCBMCBins[typeT][i],sizeof(INT_VECTOR3));
		}
		
		for (i=0; i<NumberOfEBCBMCBins; i++) {
			SizeOfEBCBMCBins[typeT][i]=0;
		}
		for (i=0; i<t; i++) {
			k=(int)(TempEnergy[i]*NumberOfEBCBMCBins);
			if (k==NumberOfEBCBMCBins) {
				EBCBMCValue[typeT][k-1][SizeOfEBCBMCBins[typeT][k-1]]=TempEnergy[i];
				EBCBMCGridLocation[typeT][k-1][SizeOfEBCBMCBins[typeT][k-1]].x=GridLocation[i].x;
				EBCBMCGridLocation[typeT][k-1][SizeOfEBCBMCBins[typeT][k-1]].y=GridLocation[i].y;
				EBCBMCGridLocation[typeT][k-1][SizeOfEBCBMCBins[typeT][k-1]].z=GridLocation[i].z;
				SizeOfEBCBMCBins[typeT][k-1]++;
			}
			else {
				EBCBMCValue[typeT][k][SizeOfEBCBMCBins[typeT][k]]=TempEnergy[i];
				EBCBMCGridLocation[typeT][k][SizeOfEBCBMCBins[typeT][k]].x=GridLocation[i].x;
				EBCBMCGridLocation[typeT][k][SizeOfEBCBMCBins[typeT][k]].y=GridLocation[i].y;
				EBCBMCGridLocation[typeT][k][SizeOfEBCBMCBins[typeT][k]].z=GridLocation[i].z;
				SizeOfEBCBMCBins[typeT][k]++;
			}
		}
	printf("SumVDWGridRosenbluth = %f\n",SumVDWGridRosenbluth[typeT]);
	}
	printf("Finished EBCBMC Allocation\n");
}

int WriteEBCBMCBins(int l)
{
	int i,m;
	int ngrid,typeT;
	FILE *FilePtr;
	char buffer[256];
	
	typeT=GridTypeList[l];
	
	ngrid=NumberOfDelaunayGridPoints.x*NumberOfDelaunayGridPoints.y*NumberOfDelaunayGridPoints.z;
	
	sprintf(buffer,"%s/share/raspa/grids",RASPA_DIRECTORY);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s",RASPA_DIRECTORY,ForceField);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s/%s",RASPA_DIRECTORY,ForceField,Framework[CurrentSystem].Name[0]);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf",RASPA_DIRECTORY,ForceField,
			Framework[CurrentSystem].Name[0],(double)SpacingDelaunayGrid);
	mkdir(buffer,S_IRWXU);
	sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%s_%s_EBCBMC.bin",
			RASPA_DIRECTORY,
			ForceField,
			Framework[CurrentSystem].Name[0],
			(double)SpacingDelaunayGrid,
			Framework[CurrentSystem].Name[0],
			PseudoAtoms[GridTypeList[l]].Name);
	FilePtr=fopen(buffer,"w");
	fwrite(&SumVDWGridRosenbluth[typeT],1,sizeof(REAL),FilePtr);	
	
	for (m=0;m<NumberOfEBCBMCBins;m++)
	{
		fwrite(&SizeOfEBCBMCBins[typeT][m],1,sizeof(int),FilePtr);
	}
	for (m=0;m<NumberOfEBCBMCBins;m++)
	{
		for (i=0;i<SizeOfEBCBMCBins[typeT][m];i++)
		{
			fwrite(&EBCBMCValue[typeT][m][i],1,sizeof(float),FilePtr);
			fwrite(&EBCBMCGridLocation[typeT][m][i],1,sizeof(INT_VECTOR3),FilePtr);
		}
	}
	fclose(FilePtr);
	return 0;
}

void ReadEBCBMCBins(void)
{
	int i,l,m,typeT;
	FILE *FilePtr;
	char buffer[256];
	
	printf("Reading VDW grid\n");
	
	EBCBMCValue=(float***)calloc(NumberOfGrids,sizeof(float**));
	EBCBMCGridLocation=(INT_VECTOR3***)calloc(NumberOfGrids,sizeof(INT_VECTOR3**));
	SumVDWGridRosenbluth=(REAL*)calloc(NumberOfGrids,sizeof(REAL));
	SizeOfEBCBMCBins=(int**)calloc(NumberOfGrids,sizeof(int*));
	
	for(l=0;l<NumberOfGrids;l++)
	{
		typeT=GridTypeList[l];
		SizeOfEBCBMCBins[typeT]=(int*)calloc(NumberOfEBCBMCBins,sizeof(int));
		EBCBMCValue[typeT]=(float**)calloc(NumberOfEBCBMCBins,sizeof(float*));
		EBCBMCGridLocation[typeT]=(INT_VECTOR3**)calloc(NumberOfEBCBMCBins,sizeof(INT_VECTOR3*));
		
		sprintf(buffer,"%s/share/raspa/grids",RASPA_DIRECTORY);
		sprintf(buffer,"%s/share/raspa/grids/%s/%s/%lf/%s_%s_EBCBMC.bin",
				RASPA_DIRECTORY,
				ForceField,
				Framework[CurrentSystem].Name[0],
				(double)SpacingVDWGrid,
				Framework[CurrentSystem].Name[0],
				PseudoAtoms[GridTypeList[l]].Name);
		printf("Opening: %s\n",buffer);
		if(!(FilePtr=fopen(buffer,"r")))
		{
			printf("Error:  file %s does not exists.\n",buffer);
			exit(1);
		}
		fread(&SumVDWGridRosenbluth[typeT],1,sizeof(REAL),FilePtr);
		
		for (m=0;m<NumberOfEBCBMCBins;m++)
		{
			fread(&SizeOfEBCBMCBins[typeT][m],1,sizeof(int),FilePtr);
			EBCBMCValue[typeT][m]=(float*)calloc(SizeOfEBCBMCBins[typeT][m],sizeof(float));
			EBCBMCGridLocation[typeT][m]=(INT_VECTOR3*)calloc(SizeOfEBCBMCBins[typeT][m],sizeof(INT_VECTOR3));
		}	
		
		for (m=0;m<NumberOfEBCBMCBins;m++)
		{
			for (i=0;i<SizeOfEBCBMCBins[typeT][m];i++)
			{
				fread(&EBCBMCValue[typeT][m][i],1,sizeof(float),FilePtr);
				fread(&EBCBMCGridLocation[typeT][m][i],1,sizeof(INT_VECTOR3),FilePtr);
			}
		}	
		
		fclose(FilePtr);
	}
}
