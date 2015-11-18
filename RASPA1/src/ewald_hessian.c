/*****************************************************************************************************
    ewald_hessian.c -  description
    ------------------------------
    Copyright © 2006-2011 David Dubbeldam, Sofia Calero, Donald E. Ellis 
                          and Randall Q. Snurr. All Rights Reserved.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/depa/webdex/quimfis/miembros/Web_Sofia/Sofia.htm
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/
 *****************************************************************************************************/


#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "simulation.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "utils.h"
#include "ewald.h"
#include "ewald_exclusion.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "potentials.h"
#include "spectra.h"
#include "rigid.h"
#include "minimization.h"


// Hessian: Strain - Center of mass (part I)
// =========================================
static inline void HessianAtomicPositionStrain(REAL_MATRIX HessianMatrix,int index_i,REAL fac,REAL_MATRIX3x3 Theta,VECTOR Rk)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i>=0)
      {
        // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
        // ================================================================
        HessianMatrix.element[index_i][n  ]-=fac*(Theta.ax*Rk.x+Theta.by*Rk.x+Theta.cz*Rk.x);    // xx x + yy x + zz x
        HessianMatrix.element[index_i+1][n]-=fac*(Theta.ax*Rk.y+Theta.by*Rk.y+Theta.cz*Rk.y);    // xx y + yy y + zz y
        HessianMatrix.element[index_i+2][n]-=fac*(Theta.ax*Rk.z+Theta.by*Rk.z+Theta.cz*Rk.z);    // xx z + yy z + zz z
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_i>=0)
          {
            // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[index_i][n  ]-=fac*Theta.ax*Rk.x;    // xx x
            HessianMatrix.element[index_i][n+1]-=fac*Theta.by*Rk.x;    // yy x
            HessianMatrix.element[index_i][n+2]-=fac*Theta.cz*Rk.x;    // xz x

            HessianMatrix.element[index_i+1][n  ]-=fac*Theta.ax*Rk.y;  // xx y
            HessianMatrix.element[index_i+1][n+1]-=fac*Theta.by*Rk.y;  // yy y
            HessianMatrix.element[index_i+1][n+2]-=fac*Theta.cz*Rk.y;  // zz y

            HessianMatrix.element[index_i+2][n  ]-=fac*Theta.ax*Rk.z;  // xx z
            HessianMatrix.element[index_i+2][n+1]-=fac*Theta.by*Rk.z;  // yy z
            HessianMatrix.element[index_i+2][n+2]-=fac*Theta.cz*Rk.z;  // zz z
          }
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          if(index_i>=0)
          {
            // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[index_i][n  ]-=fac*Theta.ax*Rk.x;    // xx x
            HessianMatrix.element[index_i][n+1]-=fac*Theta.ay*Rk.x;    // xy x
            HessianMatrix.element[index_i][n+2]-=fac*Theta.az*Rk.x;    // xz x
            HessianMatrix.element[index_i][n+3]-=fac*Theta.by*Rk.x;    // yy x
            HessianMatrix.element[index_i][n+4]-=fac*Theta.bz*Rk.x;    // yz x
            HessianMatrix.element[index_i][n+5]-=fac*Theta.cz*Rk.x;    // xz x

            HessianMatrix.element[index_i+1][n  ]-=fac*Theta.ax*Rk.y;  // xx y
            HessianMatrix.element[index_i+1][n+1]-=fac*Theta.ay*Rk.y;  // xy y
            HessianMatrix.element[index_i+1][n+2]-=fac*Theta.az*Rk.y;  // xz y
            HessianMatrix.element[index_i+1][n+3]-=fac*Theta.by*Rk.y;  // yy y
            HessianMatrix.element[index_i+1][n+4]-=fac*Theta.bz*Rk.y;  // yz y
            HessianMatrix.element[index_i+1][n+5]-=fac*Theta.cz*Rk.y;  // yz y

            HessianMatrix.element[index_i+2][n  ]-=fac*Theta.ax*Rk.z;  // xx z
            HessianMatrix.element[index_i+2][n+1]-=fac*Theta.ay*Rk.z;  // xy z
            HessianMatrix.element[index_i+2][n+2]-=fac*Theta.az*Rk.z;  // xz z
            HessianMatrix.element[index_i+2][n+3]-=fac*Theta.by*Rk.z;  // yy z
            HessianMatrix.element[index_i+2][n+4]-=fac*Theta.bz*Rk.z;  // yz z
            HessianMatrix.element[index_i+2][n+5]-=fac*Theta.cz*Rk.z;  // zz z
          }
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i>=0)
              {
                // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i][n  ]-=fac*Theta.ax*Rk.x;    // xx x
                HessianMatrix.element[index_i][n+1]-=fac*Theta.by*Rk.x;    // yy x
                HessianMatrix.element[index_i][n+2]-=fac*Theta.bz*Rk.x;    // yz x
                HessianMatrix.element[index_i][n+3]-=fac*Theta.cz*Rk.x;    // xz x

                HessianMatrix.element[index_i+1][n  ]-=fac*Theta.ax*Rk.y;  // xx y
                HessianMatrix.element[index_i+1][n+1]-=fac*Theta.by*Rk.y;  // yy y
                HessianMatrix.element[index_i+1][n+2]-=fac*Theta.bz*Rk.y;  // yz y
                HessianMatrix.element[index_i+1][n+3]-=fac*Theta.cz*Rk.y;  // yz y

                HessianMatrix.element[index_i+2][n  ]-=fac*Theta.ax*Rk.z;  // xx z
                HessianMatrix.element[index_i+2][n+1]-=fac*Theta.by*Rk.z;  // yy z
                HessianMatrix.element[index_i+2][n+2]-=fac*Theta.bz*Rk.z;  // yz z
                HessianMatrix.element[index_i+2][n+3]-=fac*Theta.cz*Rk.z;  // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i>=0)
              {
                // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i][n  ]-=fac*Theta.ax*Rk.x;    // xx x
                HessianMatrix.element[index_i][n+1]-=fac*Theta.az*Rk.x;    // xz x
                HessianMatrix.element[index_i][n+2]-=fac*Theta.by*Rk.x;    // yy x
                HessianMatrix.element[index_i][n+3]-=fac*Theta.cz*Rk.x;    // xz x

                HessianMatrix.element[index_i+1][n  ]-=fac*Theta.ax*Rk.y;  // xx y
                HessianMatrix.element[index_i+1][n+1]-=fac*Theta.az*Rk.y;  // xz y
                HessianMatrix.element[index_i+1][n+2]-=fac*Theta.by*Rk.y;  // yy y
                HessianMatrix.element[index_i+1][n+3]-=fac*Theta.cz*Rk.y;  // yz y

                HessianMatrix.element[index_i+2][n  ]-=fac*Theta.ax*Rk.z;  // xx z
                HessianMatrix.element[index_i+2][n+1]-=fac*Theta.az*Rk.z;  // xz z
                HessianMatrix.element[index_i+2][n+2]-=fac*Theta.by*Rk.z;  // yy z
                HessianMatrix.element[index_i+2][n+3]-=fac*Theta.cz*Rk.z;  // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i>=0)
              {
                // the first term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i][n  ]-=fac*Theta.ax*Rk.x;    // xx x
                HessianMatrix.element[index_i][n+1]-=fac*Theta.ay*Rk.x;    // xy x
                HessianMatrix.element[index_i][n+2]-=fac*Theta.by*Rk.x;    // yy x
                HessianMatrix.element[index_i][n+3]-=fac*Theta.cz*Rk.x;    // xz x

                HessianMatrix.element[index_i+1][n  ]-=fac*Theta.ax*Rk.y;  // xx y
                HessianMatrix.element[index_i+1][n+1]-=fac*Theta.ay*Rk.y;  // xy y
                HessianMatrix.element[index_i+1][n+2]-=fac*Theta.by*Rk.y;  // yy y
                HessianMatrix.element[index_i+1][n+3]-=fac*Theta.cz*Rk.y;  // yz y

                HessianMatrix.element[index_i+2][n  ]-=fac*Theta.ax*Rk.z;  // xx z
                HessianMatrix.element[index_i+2][n+1]-=fac*Theta.ay*Rk.z;  // xy z
                HessianMatrix.element[index_i+2][n+2]-=fac*Theta.by*Rk.z;  // yy z
                HessianMatrix.element[index_i+2][n+3]-=fac*Theta.cz*Rk.z;  // zz z
              }
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

// Hessian: Strain - Center of mass (part II)
// ==========================================
static inline void HessianCenterOfMassStrainI(REAL_MATRIX HessianMatrix,int index_i,REAL f2_I,VECTOR posA,VECTOR comA,VECTOR Rk)
{
  int n;
  VECTOR dI;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i>=0)
      {
        // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
        // =================================================================
        HessianMatrix.element[index_i  ][n]-=f2_I*(Rk.x*dI.x*Rk.x+Rk.y*dI.y*Rk.x+Rk.z*dI.z*Rk.x); // xx x + yy x + zz x
        HessianMatrix.element[index_i+1][n]-=f2_I*(Rk.x*dI.x*Rk.y+Rk.y*dI.y*Rk.y+Rk.z*dI.z*Rk.y); // xx y + yy y + zz y
        HessianMatrix.element[index_i+2][n]-=f2_I*(Rk.x*dI.x*Rk.z+Rk.y*dI.y*Rk.z+Rk.z*dI.z*Rk.z); // xx z + yy z + zz z
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_i>=0)
          {
            // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // =================================================================
            HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
            HessianMatrix.element[index_i  ][n+1]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
            HessianMatrix.element[index_i  ][n+2]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

            HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
            HessianMatrix.element[index_i+1][n+1]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
            HessianMatrix.element[index_i+1][n+2]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

            HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
            HessianMatrix.element[index_i+2][n+1]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
            HessianMatrix.element[index_i+2][n+2]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
          }
          break;
        case REGULAR:
          if(index_i>=0)
          {
            // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // =================================================================
            HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
            HessianMatrix.element[index_i  ][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.x; // xy x
            HessianMatrix.element[index_i  ][n+2]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.x; // xz x
            HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
            HessianMatrix.element[index_i  ][n+4]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.x; // yz x
            HessianMatrix.element[index_i  ][n+5]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

            HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
            HessianMatrix.element[index_i+1][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.y; // xy y
            HessianMatrix.element[index_i+1][n+2]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.y; // xz y
            HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
            HessianMatrix.element[index_i+1][n+4]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.y; // yz y
            HessianMatrix.element[index_i+1][n+5]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

            HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
            HessianMatrix.element[index_i+2][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.z; // xy z
            HessianMatrix.element[index_i+2][n+2]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.z; // xz z
            HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
            HessianMatrix.element[index_i+2][n+4]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.z; // yz z
            HessianMatrix.element[index_i+2][n+5]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          if(index_i>=0)
          {
            // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // =================================================================
            HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
            HessianMatrix.element[index_i  ][n+1]-=f2_I*Rk.x*dI.y*Rk.x; // xy x
            HessianMatrix.element[index_i  ][n+2]-=f2_I*Rk.x*dI.z*Rk.x; // xz x
            HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
            HessianMatrix.element[index_i  ][n+4]-=f2_I*Rk.y*dI.z*Rk.x; // yz x
            HessianMatrix.element[index_i  ][n+5]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

            HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
            HessianMatrix.element[index_i+1][n+1]-=f2_I*Rk.x*dI.y*Rk.y; // xy y
            HessianMatrix.element[index_i+1][n+2]-=f2_I*Rk.x*dI.z*Rk.y; // xz y
            HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
            HessianMatrix.element[index_i+1][n+4]-=f2_I*Rk.y*dI.z*Rk.y; // yz y
            HessianMatrix.element[index_i+1][n+5]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

            HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
            HessianMatrix.element[index_i+2][n+1]-=f2_I*Rk.x*dI.y*Rk.z; // xy z
            HessianMatrix.element[index_i+2][n+2]-=f2_I*Rk.x*dI.z*Rk.z; // xz z
            HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
            HessianMatrix.element[index_i+2][n+4]-=f2_I*Rk.y*dI.z*Rk.z; // yz z
            HessianMatrix.element[index_i+2][n+5]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i>=0)
              {
                // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
                HessianMatrix.element[index_i  ][n+1]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
                HessianMatrix.element[index_i  ][n+2]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.x; // yz x
                HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

                HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
                HessianMatrix.element[index_i+1][n+1]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
                HessianMatrix.element[index_i+1][n+2]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.y; // yz y
                HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

                HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
                HessianMatrix.element[index_i+2][n+1]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
                HessianMatrix.element[index_i+2][n+2]-=0.5*f2_I*(Rk.y*dI.z+Rk.z*dI.y)*Rk.z; // yz z
                HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i>=0)
              {
                // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
                HessianMatrix.element[index_i  ][n+1]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.x; // xz x
                HessianMatrix.element[index_i  ][n+2]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
                HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

                HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
                HessianMatrix.element[index_i+1][n+1]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.y; // xz y
                HessianMatrix.element[index_i+1][n+2]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
                HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

                HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
                HessianMatrix.element[index_i+2][n+1]-=0.5*f2_I*(Rk.x*dI.z+Rk.z*dI.x)*Rk.z; // xz z
                HessianMatrix.element[index_i+2][n+2]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
                HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i>=0)
              {
                // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x;                 // xx x
                HessianMatrix.element[index_i  ][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.x; // xy x
                HessianMatrix.element[index_i  ][n+2]-=f2_I*Rk.y*dI.y*Rk.x;                 // yy x
                HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.z*dI.z*Rk.x;                 // zz x

                HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y;                 // xx y
                HessianMatrix.element[index_i+1][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.y; // xy y
                HessianMatrix.element[index_i+1][n+2]-=f2_I*Rk.y*dI.y*Rk.y;                 // yy y
                HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.z*dI.z*Rk.y;                 // zz y

                HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z;                 // xx z
                HessianMatrix.element[index_i+2][n+1]-=0.5*f2_I*(Rk.x*dI.y+Rk.y*dI.x)*Rk.z; // xy z
                HessianMatrix.element[index_i+2][n+2]-=f2_I*Rk.y*dI.y*Rk.z;                 // yy z
                HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.z*dI.z*Rk.z;                 // zz z
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i>=0)
              {
                // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
                HessianMatrix.element[index_i  ][n+1]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
                HessianMatrix.element[index_i  ][n+2]-=f2_I*Rk.y*dI.z*Rk.x; // yz x
                HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

                HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
                HessianMatrix.element[index_i+1][n+1]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
                HessianMatrix.element[index_i+1][n+2]-=f2_I*Rk.y*dI.z*Rk.y; // yz y
                HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

                HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
                HessianMatrix.element[index_i+2][n+1]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
                HessianMatrix.element[index_i+2][n+2]-=f2_I*Rk.y*dI.z*Rk.z; // yz z
                HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i>=0)
              {
                // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
                HessianMatrix.element[index_i  ][n+1]-=f2_I*Rk.x*dI.z*Rk.x; // xz x
                HessianMatrix.element[index_i  ][n+2]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
                HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

                HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
                HessianMatrix.element[index_i+1][n+1]-=f2_I*Rk.x*dI.z*Rk.y; // xz y
                HessianMatrix.element[index_i+1][n+2]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
                HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

                HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
                HessianMatrix.element[index_i+2][n+1]-=f2_I*Rk.x*dI.z*Rk.z; // xz z
                HessianMatrix.element[index_i+2][n+2]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
                HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i>=0)
              {
                // the second term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                HessianMatrix.element[index_i  ][n  ]-=f2_I*Rk.x*dI.x*Rk.x; // xx x
                HessianMatrix.element[index_i  ][n+1]-=f2_I*Rk.x*dI.y*Rk.x; // xy x
                HessianMatrix.element[index_i  ][n+2]-=f2_I*Rk.y*dI.y*Rk.x; // yy x
                HessianMatrix.element[index_i  ][n+3]-=f2_I*Rk.z*dI.z*Rk.x; // zz x

                HessianMatrix.element[index_i+1][n  ]-=f2_I*Rk.x*dI.x*Rk.y; // xx y
                HessianMatrix.element[index_i+1][n+1]-=f2_I*Rk.x*dI.y*Rk.y; // xy y
                HessianMatrix.element[index_i+1][n+2]-=f2_I*Rk.y*dI.y*Rk.y; // yy y
                HessianMatrix.element[index_i+1][n+3]-=f2_I*Rk.z*dI.z*Rk.y; // zz y

                HessianMatrix.element[index_i+2][n  ]-=f2_I*Rk.x*dI.x*Rk.z; // xx z
                HessianMatrix.element[index_i+2][n+1]-=f2_I*Rk.x*dI.y*Rk.z; // xy z
                HessianMatrix.element[index_i+2][n+2]-=f2_I*Rk.y*dI.y*Rk.z; // yy z
                HessianMatrix.element[index_i+2][n+3]-=f2_I*Rk.z*dI.z*Rk.z; // zz z
              }
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

// Hessian: Strain - Center of mass (part III)
// ===========================================
static inline void HessianCenterOfMassStrainJ(REAL_MATRIX HessianMatrix,int index_i,REAL f2_IJ,VECTOR posB,VECTOR comB,VECTOR Rk)
{
  int n;
  VECTOR dJ;

  n=NumberOfCoordinatesMinimizationVariables;

  dJ.x=posB.x-comB.x;
  dJ.y=posB.y-comB.y;
  dJ.z=posB.z-comB.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i>=0)
      {
        // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
        // ================================================================
        HessianMatrix.element[index_i  ][n]+=f2_IJ*(Rk.x*dJ.x*Rk.x+Rk.y*dJ.y*Rk.x+Rk.z*dJ.z*Rk.x); // xx x + yy x + zz x
        HessianMatrix.element[index_i+1][n]+=f2_IJ*(Rk.x*dJ.x*Rk.y+Rk.y*dJ.y*Rk.y+Rk.z*dJ.z*Rk.y); // xx y + yy y + zz y
        HessianMatrix.element[index_i+2][n]+=f2_IJ*(Rk.x*dJ.x*Rk.z+Rk.y*dJ.y*Rk.z+Rk.z*dJ.z*Rk.z); // xx z + yy z + zz z
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_i>=0)
          {
            // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[index_i  ][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
            HessianMatrix.element[index_i  ][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
            HessianMatrix.element[index_i  ][n+2]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x

            HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
            HessianMatrix.element[index_i+1][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
            HessianMatrix.element[index_i+1][n+2]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

            HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
            HessianMatrix.element[index_i+2][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
            HessianMatrix.element[index_i+2][n+2]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
          }
          break;
        case REGULAR:
          if(index_i>=0)
          {
            // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[index_i][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                   // xx x
            HessianMatrix.element[index_i][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.x;   // xy x
            HessianMatrix.element[index_i][n+2]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.x;   // xz x
            HessianMatrix.element[index_i][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.x;                   // yy x
            HessianMatrix.element[index_i][n+4]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.x;   // yz x
            HessianMatrix.element[index_i][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.x;                   // zz x

            HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
            HessianMatrix.element[index_i+1][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.y; // xy y
            HessianMatrix.element[index_i+1][n+2]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.y; // xz y
            HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
            HessianMatrix.element[index_i+1][n+4]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.y; // yz y
            HessianMatrix.element[index_i+1][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

            HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
            HessianMatrix.element[index_i+2][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.z; // xy z
            HessianMatrix.element[index_i+2][n+2]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.z; // xz z
            HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
            HessianMatrix.element[index_i+2][n+4]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.z; // yz z
            HessianMatrix.element[index_i+2][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          if(index_i>=0)
          {
            // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[index_i  ][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
            HessianMatrix.element[index_i  ][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.x; // xy x
            HessianMatrix.element[index_i  ][n+2]+=f2_IJ*Rk.x*dJ.z*Rk.x; // xz x
            HessianMatrix.element[index_i  ][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
            HessianMatrix.element[index_i  ][n+4]+=f2_IJ*Rk.y*dJ.z*Rk.x; // yz x
            HessianMatrix.element[index_i  ][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x

            HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
            HessianMatrix.element[index_i+1][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.y; // xy y
            HessianMatrix.element[index_i+1][n+2]+=f2_IJ*Rk.x*dJ.z*Rk.y; // xz y
            HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
            HessianMatrix.element[index_i+1][n+4]+=f2_IJ*Rk.y*dJ.z*Rk.y; // yz y
            HessianMatrix.element[index_i+1][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

            HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
            HessianMatrix.element[index_i+2][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.z; // xy z
            HessianMatrix.element[index_i+2][n+2]+=f2_IJ*Rk.x*dJ.z*Rk.z; // xz z
            HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
            HessianMatrix.element[index_i+2][n+4]+=f2_IJ*Rk.y*dJ.z*Rk.z; // yz z
            HessianMatrix.element[index_i+2][n+5]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i>=0)
              {
                // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                   // xx x
                HessianMatrix.element[index_i][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.x;                   // yy x
                HessianMatrix.element[index_i][n+2]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.x;   // yz x
                HessianMatrix.element[index_i][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x;                   // zz x

                HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
                HessianMatrix.element[index_i+1][n+2]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.y; // yz y
                HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
                HessianMatrix.element[index_i+2][n+2]+=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y)*Rk.z; // yz z
                HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i>=0)
              {
                // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                   // xx x
                HessianMatrix.element[index_i][n+1]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.x;   // xz x
                HessianMatrix.element[index_i][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x;                   // yy x
                HessianMatrix.element[index_i][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x;                   // zz x

                HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.y; // xz y
                HessianMatrix.element[index_i+1][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
                HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x)*Rk.z; // xz z
                HessianMatrix.element[index_i+2][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
                HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i>=0)
              {
                // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x;                   // xx x
                HessianMatrix.element[index_i][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.x;   // xy x
                HessianMatrix.element[index_i][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x;                   // yy x
                HessianMatrix.element[index_i][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x;                   // zz x

                HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y;                 // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.y; // xy y
                HessianMatrix.element[index_i+1][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y;                 // yy y
                HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y;                 // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z;                 // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x)*Rk.z; // xy z
                HessianMatrix.element[index_i+2][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z;                 // yy z
                HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z;                 // zz z
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i>=0)
              {
                // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i  ][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
                HessianMatrix.element[index_i  ][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
                HessianMatrix.element[index_i  ][n+2]+=f2_IJ*Rk.y*dJ.z*Rk.x; // yz x
                HessianMatrix.element[index_i  ][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x

                HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
                HessianMatrix.element[index_i+1][n+2]+=f2_IJ*Rk.y*dJ.z*Rk.y; // yz y
                HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
                HessianMatrix.element[index_i+2][n+2]+=f2_IJ*Rk.y*dJ.z*Rk.z; // yz z
                HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i>=0)
              {
                // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i  ][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
                HessianMatrix.element[index_i  ][n+1]+=f2_IJ*Rk.x*dJ.z*Rk.x; // xz x
                HessianMatrix.element[index_i  ][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
                HessianMatrix.element[index_i  ][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x
          
                HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2_IJ*Rk.x*dJ.z*Rk.y; // xz y
                HessianMatrix.element[index_i+1][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
                HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2_IJ*Rk.x*dJ.z*Rk.z; // xz z
                HessianMatrix.element[index_i+2][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
                HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i>=0)
              {
                // the third term of Eq. 41 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                HessianMatrix.element[index_i  ][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.x; // xx x
                HessianMatrix.element[index_i  ][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.x; // xy x
                HessianMatrix.element[index_i  ][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.x; // yy x
                HessianMatrix.element[index_i  ][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.x; // zz x
          
                HessianMatrix.element[index_i+1][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.y; // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.y; // xy y
                HessianMatrix.element[index_i+1][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.y; // yy y
                HessianMatrix.element[index_i+1][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.y; // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2_IJ*Rk.x*dJ.x*Rk.z; // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2_IJ*Rk.x*dJ.y*Rk.z; // xy z
                HessianMatrix.element[index_i+2][n+2]+=f2_IJ*Rk.y*dJ.y*Rk.z; // yy z
                HessianMatrix.element[index_i+2][n+3]+=f2_IJ*Rk.z*dJ.z*Rk.z; // zz z
              }
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}


// Hessian: Strain - Orientation (part I)
// ======================================
static inline void HessianOrientationStrainI(REAL_MATRIX HessianMatrix,int index_i2,int index1_rigid,
         REAL f1I,REAL f2I,VECTOR posA,VECTOR comA,VECTOR Rk,REAL_MATRIX3x3 Theta)
{
  int n;
  VECTOR dot_product;
  VECTOR veci1,veci2,veci3;
  REAL temp;
  VECTOR dI;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  veci1=DVecX[index1_rigid];
  veci2=DVecY[index1_rigid];
  veci3=DVecZ[index1_rigid];

  dot_product.x=Rk.x*veci1.x+Rk.y*veci1.y+Rk.z*veci1.z;
  dot_product.y=Rk.x*veci2.x+Rk.y*veci2.y+Rk.z*veci2.z;
  dot_product.z=Rk.x*veci3.x+Rk.y*veci3.y+Rk.z*veci3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i2>=0)
      {
        temp=f1I*Theta.ax+f2I*dI.x*Rk.x;
        HessianMatrix.element[index_i2][n]-=temp*dot_product.x+f1I*veci1.x*Rk.x;     // xx x
        HessianMatrix.element[index_i2+1][n]-=temp*dot_product.y+f1I*veci2.x*Rk.x;   // xx y
        HessianMatrix.element[index_i2+2][n]-=temp*dot_product.z+f1I*veci3.x*Rk.x;   // xx z

        temp=f1I*Theta.by+f2I*dI.y*Rk.y;
        HessianMatrix.element[index_i2][n]-=temp*dot_product.x+f1I*veci1.y*Rk.y;   // yy z
        HessianMatrix.element[index_i2+1][n]-=temp*dot_product.y+f1I*veci2.y*Rk.y; // yy y
        HessianMatrix.element[index_i2+2][n]-=temp*dot_product.z+f1I*veci3.y*Rk.y; // yy z

        temp=f1I*Theta.cz+f2I*dI.z*Rk.z;
        HessianMatrix.element[index_i2][n]-=temp*dot_product.x+f1I*veci1.z*Rk.z;   // zz x
        HessianMatrix.element[index_i2+1][n]-=temp*dot_product.y+f1I*veci2.z*Rk.z; // zz y
        HessianMatrix.element[index_i2+2][n]-=temp*dot_product.z+f1I*veci3.z*Rk.z; // zz z
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_i2>=0)
          {
            temp=f1I*Theta.ax+f2I*dI.x*Rk.x;
            HessianMatrix.element[index_i2][n]-=temp*dot_product.x+f1I*veci1.x*Rk.x;     // xx x
            HessianMatrix.element[index_i2+1][n]-=temp*dot_product.y+f1I*veci2.x*Rk.x;   // xx y
            HessianMatrix.element[index_i2+2][n]-=temp*dot_product.z+f1I*veci3.x*Rk.x;   // xx z

            temp=f1I*Theta.by+f2I*dI.y*Rk.y;
            HessianMatrix.element[index_i2][n+1]-=temp*dot_product.x+f1I*veci1.y*Rk.y;   // yy z
            HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*veci2.y*Rk.y; // yy y
            HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*veci3.y*Rk.y; // yy z

            temp=f1I*Theta.cz+f2I*dI.z*Rk.z;
            HessianMatrix.element[index_i2][n+2]-=temp*dot_product.x+f1I*veci1.z*Rk.z;   // zz x
            HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*veci2.z*Rk.z; // zz y
            HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*veci3.z*Rk.z; // zz z
          }
          break;
        case REGULAR:
          if(index_i2>=0)
          {
            // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ===================================================================================
            temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
            HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
            HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
            HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

            temp=f2I*0.5*(Rk.x*dI.y+Rk.y*dI.x)+f1I*Theta.ay;
            HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*0.5*(Rk.x*veci1.y+Rk.y*veci1.x); // xy x
            HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*0.5*(Rk.x*veci2.y+Rk.y*veci2.x); // xy y
            HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*0.5*(Rk.x*veci3.y+Rk.y*veci3.x); // xy z

            temp=f2I*0.5*(Rk.x*dI.z+Rk.z*dI.x)+f1I*Theta.az;
            HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*0.5*(veci1.z*Rk.x+veci1.x*Rk.z); // xz x
            HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*0.5*(veci2.z*Rk.x+veci2.x*Rk.z); // xz y
            HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*0.5*(veci3.z*Rk.x+veci3.x*Rk.z); // xz z

            temp=f2I*Rk.y*dI.y+f1I*Theta.by;
            HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
            HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
            HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

            temp=f2I*0.5*(Rk.y*dI.z+Rk.z*dI.y)+f1I*Theta.bz;
            HessianMatrix.element[index_i2  ][n+4]-=temp*dot_product.x+f1I*0.5*(Rk.y*veci1.z+Rk.z*veci1.y); // yz x
            HessianMatrix.element[index_i2+1][n+4]-=temp*dot_product.y+f1I*0.5*(Rk.y*veci2.z+Rk.z*veci2.y); // yz y
            HessianMatrix.element[index_i2+2][n+4]-=temp*dot_product.z+f1I*0.5*(Rk.y*veci3.z+Rk.z*veci3.y); // yz z

            temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
            HessianMatrix.element[index_i2  ][n+5]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
            HessianMatrix.element[index_i2+1][n+5]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
            HessianMatrix.element[index_i2+2][n+5]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          if(index_i2>=0)
          {
            // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ===================================================================================
            temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
            HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
            HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
            HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

            temp=f2I*Rk.x*dI.y+f1I*Theta.ay;
            HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*Rk.x*veci1.y; // xy x
            HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*Rk.x*veci2.y; // xy y
            HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*Rk.x*veci3.y; // xy z

            temp=f2I*Rk.x*dI.z+f1I*Theta.az;
            HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*Rk.x*veci1.z; // xz x
            HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*Rk.x*veci2.z; // xz y
            HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*Rk.x*veci3.z; // xz z

            temp=f2I*Rk.y*dI.y+f1I*Theta.by;
            HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
            HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
            HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

            temp=f2I*Rk.y*dI.z+f1I*Theta.bz;
            HessianMatrix.element[index_i2  ][n+4]-=temp*dot_product.x+f1I*Rk.y*veci1.z; // yz x
            HessianMatrix.element[index_i2+1][n+4]-=temp*dot_product.y+f1I*Rk.y*veci2.z; // yz y
            HessianMatrix.element[index_i2+2][n+4]-=temp*dot_product.z+f1I*Rk.y*veci3.z; // yz z

            temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
            HessianMatrix.element[index_i2  ][n+5]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
            HessianMatrix.element[index_i2+1][n+5]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
            HessianMatrix.element[index_i2+2][n+5]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i2>=0)
              {
                // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ===================================================================================
                temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
                HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
                HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

                temp=f2I*Rk.y*dI.y+f1I*Theta.by;
                HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
                HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
                HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

                temp=f2I*0.5*(Rk.y*dI.z+Rk.z*dI.y)+f1I*Theta.bz;
                HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*0.5*(Rk.y*veci1.z+Rk.z*veci1.y); // yz x
                HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*0.5*(Rk.y*veci2.z+Rk.z*veci2.y); // yz y
                HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*0.5*(Rk.y*veci3.z+Rk.z*veci3.y); // yz z

                temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
                HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
                HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
                HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i2>=0)
              {
                // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ===================================================================================
                temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
                HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
                HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

                temp=f2I*0.5*(Rk.x*dI.z+Rk.z*dI.x)+f1I*Theta.az;
                HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*0.5*(veci1.z*Rk.x+veci1.x*Rk.z); // xz x
                HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*0.5*(veci2.z*Rk.x+veci2.x*Rk.z); // xz y
                HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*0.5*(veci3.z*Rk.x+veci3.x*Rk.z); // xz z

                temp=f2I*Rk.y*dI.y+f1I*Theta.by;
                HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
                HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

                temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
                HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
                HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
                HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i2>=0)
              {
                // the first, second, and third term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ===================================================================================
                temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
                HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
                HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

                temp=f2I*0.5*(Rk.x*dI.y+Rk.y*dI.x)+f1I*Theta.ay;
                HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*0.5*(Rk.x*veci1.y+Rk.y*veci1.x); // xy x
                HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*0.5*(Rk.x*veci2.y+Rk.y*veci2.x); // xy y
                HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*0.5*(Rk.x*veci3.y+Rk.y*veci3.x); // xy z

                temp=f2I*Rk.y*dI.y+f1I*Theta.by;
                HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
                HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

                temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
                HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
                HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
                HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i2>=0)
              {
                // the first and second term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
                HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
                HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

                temp=f2I*Rk.y*dI.y+f1I*Theta.by;
                HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
                HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
                HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

                temp=f2I*Rk.y*dI.z+f1I*Theta.bz;
                HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.z; // yz x
                HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.z; // yz y
                HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.z; // yz z

                temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
                HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
                HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
                HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i2>=0)
              {
                // the first and second term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
                HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
                HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

                temp=f2I*Rk.x*dI.z+f1I*Theta.az;
                HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*Rk.x*veci1.z; // xz x
                HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*Rk.x*veci2.z; // xz y
                HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*Rk.x*veci3.z; // xz z

                temp=f2I*Rk.y*dI.y+f1I*Theta.by;
                HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
                HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

                temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
                HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
                HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
                HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i2>=0)
              {
                // the first and second term of Eq. 42 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // ================================================================
                temp=f2I*Rk.x*dI.x+f1I*Theta.ax;
                HessianMatrix.element[index_i2  ][n  ]-=temp*dot_product.x+f1I*Rk.x*veci1.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]-=temp*dot_product.y+f1I*Rk.x*veci2.x; // xx y
                HessianMatrix.element[index_i2+2][n  ]-=temp*dot_product.z+f1I*Rk.x*veci3.x; // xx z

                temp=f2I*Rk.x*dI.y+f1I*Theta.ay;
                HessianMatrix.element[index_i2  ][n+1]-=temp*dot_product.x+f1I*Rk.x*veci1.y; // xy x
                HessianMatrix.element[index_i2+1][n+1]-=temp*dot_product.y+f1I*Rk.x*veci2.y; // xy y
                HessianMatrix.element[index_i2+2][n+1]-=temp*dot_product.z+f1I*Rk.x*veci3.y; // xy z

                temp=f2I*Rk.y*dI.y+f1I*Theta.by;
                HessianMatrix.element[index_i2  ][n+2]-=temp*dot_product.x+f1I*Rk.y*veci1.y; // yy z
                HessianMatrix.element[index_i2+1][n+2]-=temp*dot_product.y+f1I*Rk.y*veci2.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]-=temp*dot_product.z+f1I*Rk.y*veci3.y; // yy z

                temp=f2I*Rk.z*dI.z+f1I*Theta.cz;
                HessianMatrix.element[index_i2  ][n+3]-=temp*dot_product.x+f1I*Rk.z*veci1.z; // zz x
                HessianMatrix.element[index_i2+1][n+3]-=temp*dot_product.y+f1I*Rk.z*veci2.z; // zz y
                HessianMatrix.element[index_i2+2][n+3]-=temp*dot_product.z+f1I*Rk.z*veci3.z; // zz z
              }
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

// Hessian: Strain - Orientation (part II)
// =======================================
static inline void HessianOrientationStrainJ(REAL_MATRIX HessianMatrix,int index_i2,int index2_rigid,
                                             REAL f2_IJ,VECTOR posB,VECTOR comB,VECTOR Rk)
{
  int n;
  VECTOR dot_product;
  VECTOR vecj1,vecj2,vecj3;
  REAL temp;
  VECTOR dJ;

  n=NumberOfCoordinatesMinimizationVariables;

  dJ.x=posB.x-comB.x;
  dJ.y=posB.y-comB.y;
  dJ.z=posB.z-comB.z;

  vecj1=DVecX[index2_rigid];
  vecj2=DVecY[index2_rigid];
  vecj3=DVecZ[index2_rigid];

  dot_product.x=Rk.x*vecj1.x+Rk.y*vecj1.y+Rk.z*vecj1.z;
  dot_product.y=Rk.x*vecj2.x+Rk.y*vecj2.y+Rk.z*vecj2.z;
  dot_product.z=Rk.x*vecj3.x+Rk.y*vecj3.y+Rk.z*vecj3.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i2>=0)
      {
        temp=f2_IJ*Rk.x*dJ.x;
        HessianMatrix.element[index_i2  ][n]+=temp*dot_product.x; // xx x
        HessianMatrix.element[index_i2+1][n]+=temp*dot_product.y; // xx y
        HessianMatrix.element[index_i2+2][n]+=temp*dot_product.z; // xx z

        temp=f2_IJ*Rk.y*dJ.y;
        HessianMatrix.element[index_i2  ][n]+=temp*dot_product.x; // yy x
        HessianMatrix.element[index_i2+1][n]+=temp*dot_product.y; // yy y
        HessianMatrix.element[index_i2+2][n]+=temp*dot_product.z; // yy z

        temp=f2_IJ*Rk.z*dJ.z;
        HessianMatrix.element[index_i2  ][n]+=temp*dot_product.x; // zz x
        HessianMatrix.element[index_i2+1][n]+=temp*dot_product.y; // zz y
        HessianMatrix.element[index_i2+2][n]+=temp*dot_product.z; // zz z
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_i2>=0)
          {
            temp=f2_IJ*Rk.x*dJ.x;
            HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
            HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
            HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

            temp=f2_IJ*Rk.y*dJ.y;
            HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // yy x
            HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // yy y
            HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // yy z

            temp=f2_IJ*Rk.z*dJ.z;
            HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // zz x
            HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // zz y
            HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // zz z
          }
          break;
        case REGULAR:
          if(index_i2>=0)
          {
            // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // =================================================================
            temp=f2_IJ*Rk.x*dJ.x;
            HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
            HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
            HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

            temp=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x);
            HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // xy x
            HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // xy y
            HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // xy z

            temp=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x);
            HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // xz x
            HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // xz y
            HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // xz z

            temp=f2_IJ*Rk.y*dJ.y;
            HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // yy x
            HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // yy y
            HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // yy z

            temp=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y);
            HessianMatrix.element[index_i2  ][n+4]+=temp*dot_product.x; // yz x
            HessianMatrix.element[index_i2+1][n+4]+=temp*dot_product.y; // yz y
            HessianMatrix.element[index_i2+2][n+4]+=temp*dot_product.z; // yz z

            temp=f2_IJ*Rk.z*dJ.z;
            HessianMatrix.element[index_i2  ][n+5]+=temp*dot_product.x; // zz x
            HessianMatrix.element[index_i2+1][n+5]+=temp*dot_product.y; // zz y
            HessianMatrix.element[index_i2+2][n+5]+=temp*dot_product.z; // zz z
          }
          break;
        case REGULAR_UPPER_TRIANGLE:
          if(index_i2>=0)
          {
            // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // =================================================================
            temp=f2_IJ*dJ.x*Rk.x;
            HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
            HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
            HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

            temp=f2_IJ*dJ.y*Rk.x;
            HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // xy x
            HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // xy y
            HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // xy z

            temp=f2_IJ*dJ.z*Rk.x;
            HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // xz x
            HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // xz y
            HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // xz z

            temp=f2_IJ*dJ.y*Rk.y;
            HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // yy x
            HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // yy y
            HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // yy z

            temp=f2_IJ*dJ.z*Rk.y;
            HessianMatrix.element[index_i2  ][n+4]+=temp*dot_product.x; // yz x
            HessianMatrix.element[index_i2+1][n+4]+=temp*dot_product.y; // yz y
            HessianMatrix.element[index_i2+2][n+4]+=temp*dot_product.z; // yz z

            temp=f2_IJ*dJ.z*Rk.z;
            HessianMatrix.element[index_i2  ][n+5]+=temp*dot_product.x; // zz x
            HessianMatrix.element[index_i2+1][n+5]+=temp*dot_product.y; // zz y
            HessianMatrix.element[index_i2+2][n+5]+=temp*dot_product.z; // zz z
          }
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i2>=0)
              {
                // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                temp=f2_IJ*Rk.x*dJ.x;
                HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
                HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

                temp=f2_IJ*Rk.y*dJ.y;
                HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // yy x
                HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // yy y
                HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // yy z

                temp=f2_IJ*0.5*(Rk.y*dJ.z+Rk.z*dJ.y);
                HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // yz x
                HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // yz y
                HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // yz z

                temp=f2_IJ*Rk.z*dJ.z;
                HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // zz x
                HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // zz y
                HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i2>=0)
              {
                // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                temp=f2_IJ*Rk.x*dJ.x;
                HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
                HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

                temp=f2_IJ*0.5*(Rk.x*dJ.z+Rk.z*dJ.x);
                HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // xz x
                HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // xz y
                HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // xz z

                temp=f2_IJ*Rk.y*dJ.y;
                HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // yy x
                HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // yy z

                temp=f2_IJ*Rk.z*dJ.z;
                HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // zz x
                HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // zz y
                HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i2>=0)
              {
                // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                temp=f2_IJ*Rk.x*dJ.x;
                HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
                HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

                temp=f2_IJ*0.5*(Rk.x*dJ.y+Rk.y*dJ.x);
                HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // xy x
                HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // xy y
                HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // xy z

                temp=f2_IJ*Rk.y*dJ.y;
                HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // yy x
                HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // yy z

                temp=f2_IJ*Rk.z*dJ.z;
                HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // zz x
                HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // zz y
                HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // zz z
              }
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i2>=0)
              {
                // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                temp=f2_IJ*Rk.x*dJ.x;
                HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
                HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

                temp=f2_IJ*Rk.y*dJ.y;
                HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // yy x
                HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // yy y
                HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // yy z

                temp=f2_IJ*Rk.y*dJ.z;
                HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // yz x
                HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // yz y
                HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // yz z

                temp=f2_IJ*Rk.z*dJ.z;
                HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // zz x
                HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // zz y
                HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i2>=0)
              {
                // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                temp=f2_IJ*Rk.x*dJ.x;
                HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
                HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

                temp=f2_IJ*Rk.x*dJ.z;
                HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // xz x
                HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // xz y
                HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // xz z

                temp=f2_IJ*Rk.y*dJ.y;
                HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // yy x
                HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // yy z

                temp=f2_IJ*Rk.z*dJ.z;
                HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // zz x
                HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // zz y
                HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i2>=0)
              {
                // the fourth term of Eq. 43 of Ref. Dubbeldam, Krishna, Snurr, 2009
                // =================================================================
                temp=f2_IJ*Rk.x*dJ.x;
                HessianMatrix.element[index_i2  ][n  ]+=temp*dot_product.x; // xx x
                HessianMatrix.element[index_i2+1][n  ]+=temp*dot_product.y; // xx y
                HessianMatrix.element[index_i2+2][n  ]+=temp*dot_product.z; // xx z

                temp=f2_IJ*Rk.x*dJ.y;
                HessianMatrix.element[index_i2  ][n+1]+=temp*dot_product.x; // xy x
                HessianMatrix.element[index_i2+1][n+1]+=temp*dot_product.y; // xy y
                HessianMatrix.element[index_i2+2][n+1]+=temp*dot_product.z; // xy z

                temp=f2_IJ*Rk.y*dJ.y;
                HessianMatrix.element[index_i2  ][n+2]+=temp*dot_product.x; // yy x
                HessianMatrix.element[index_i2+1][n+2]+=temp*dot_product.y; // yy y
                HessianMatrix.element[index_i2+2][n+2]+=temp*dot_product.z; // yy z

                temp=f2_IJ*Rk.z*dJ.z;
                HessianMatrix.element[index_i2  ][n+3]+=temp*dot_product.x; // zz x
                HessianMatrix.element[index_i2+1][n+3]+=temp*dot_product.y; // zz y
                HessianMatrix.element[index_i2+2][n+3]+=temp*dot_product.z; // zz z
              }
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}


// Hessian: Strain - Strain (part I)
// =================================
static inline void HessianAtomicStrainStrain(REAL_MATRIX HessianMatrix,REAL f,REAL InverseLamdaSquared,
                                             VECTOR Rk,REAL Rksq,REAL_MATRIX3x3 Theta)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT: 
      HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+(1.0+1.0)+
             4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-2.0*(4.0*Rk.x*Rk.x)*InverseLamdaSquared-2.0*Theta.ax);
      HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.by+
             4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
      HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.cz+
             4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));
      HessianMatrix.element[n][n]+=f*(Theta.by*Theta.by+(1.0+1.0)+
             4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-2.0*(4.0*Rk.y*Rk.y)*InverseLamdaSquared-2.0*Theta.by);
      HessianMatrix.element[n][n]+=2.0*f*(Theta.by*Theta.cz+
             4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));
      HessianMatrix.element[n][n]+=f*(Theta.cz*Theta.cz+(1.0+1.0)+
             4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-2.0*(4.0*Rk.z*Rk.z)*InverseLamdaSquared-2.0*Theta.cz);
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+(1.0+1.0)+
                 4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-2.0*(4.0*Rk.x*Rk.x)*InverseLamdaSquared-2.0*Theta.ax);
          HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.by+
                 4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
          HessianMatrix.element[n][n]+=2.0*f*(Theta.ax*Theta.cz+
                 4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));
          HessianMatrix.element[n][n]+=f*(Theta.by*Theta.by+(1.0+1.0)+
                 4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-2.0*(4.0*Rk.y*Rk.y)*InverseLamdaSquared-2.0*Theta.by);
          HessianMatrix.element[n][n]+=2.0*f*(Theta.by*Theta.cz+
                 4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));
          HessianMatrix.element[n][n]+=f*(Theta.cz*Theta.cz+(1.0+1.0)+
                 4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-2.0*(4.0*Rk.z*Rk.z)*InverseLamdaSquared-2.0*Theta.cz);
          break;
        case ANISOTROPIC:
            HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+(1.0+1.0)+
                   4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-2.0*(4.0*Rk.x*Rk.x)*InverseLamdaSquared-2.0*Theta.ax);
            HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.by+
                   4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.cz+
                   4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n+1][n+1]+=f*(Theta.by*Theta.by+(1.0+1.0)+
                   4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-2.0*(4.0*Rk.y*Rk.y)*InverseLamdaSquared-2.0*Theta.by);
            HessianMatrix.element[n+1][n+2]+=f*(Theta.by*Theta.cz+
                   4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n+2][n+2]+=f*(Theta.cz*Theta.cz+(1.0+1.0)+
                   4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-2.0*(4.0*Rk.z*Rk.z)*InverseLamdaSquared-2.0*Theta.cz);
          break;
        case REGULAR:
            // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
            HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
            HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
            HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n][n+4]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n][n+5]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.by));
            HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.az+4.0*Rk.x*Rk.y*Rk.x*Rk.z/(SQR(Rksq))-2.0*Rk.y*Rk.z*InverseLamdaSquared-0.5*Theta.bz);
            HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared-Theta.ay);
            HessianMatrix.element[n+1][n+4]+=f*(Theta.ay*Theta.bz+4.0*Rk.x*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.z*InverseLamdaSquared-0.5*Theta.az);
            HessianMatrix.element[n+1][n+5]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+2][n+2]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.cz));
            HessianMatrix.element[n+2][n+3]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n+2][n+4]+=f*(Theta.az*Theta.bz+4.0*Rk.x*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.y*InverseLamdaSquared-0.5*Theta.ay);
            HessianMatrix.element[n+2][n+5]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared-Theta.az);

            HessianMatrix.element[n+3][n+3]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
            HessianMatrix.element[n+3][n+4]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
            HessianMatrix.element[n+3][n+5]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+4][n+4]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-0.5*(Theta.by+Theta.cz));
            HessianMatrix.element[n+4][n+5]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);

            HessianMatrix.element[n+5][n+5]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
          break;
        case REGULAR_UPPER_TRIANGLE:
            // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
            // ================================================================
            HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
            HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
            HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
            HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n][n+4]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
            HessianMatrix.element[n][n+5]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-Theta.by);
            HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.az+4.0*Rk.x*Rk.y*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.z)*InverseLamdaSquared-Theta.bz);
            HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared);
            HessianMatrix.element[n+1][n+4]+=f*(Theta.ay*Theta.bz+4.0*Rk.x*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.z*InverseLamdaSquared);
            HessianMatrix.element[n+1][n+5]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+2][n+2]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-Theta.cz);
            HessianMatrix.element[n+2][n+3]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
            HessianMatrix.element[n+2][n+4]+=f*(Theta.az*Theta.bz+4.0*Rk.x*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*Rk.x*Rk.y*InverseLamdaSquared);
            HessianMatrix.element[n+2][n+5]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared);

            HessianMatrix.element[n+3][n+3]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
            HessianMatrix.element[n+3][n+4]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
            HessianMatrix.element[n+3][n+5]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

            HessianMatrix.element[n+4][n+4]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-Theta.cz);
            HessianMatrix.element[n+4][n+5]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared);

            HessianMatrix.element[n+5][n+5]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
            break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-0.5*(Theta.by+Theta.cz));
              HessianMatrix.element[n+2][n+3]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.cz));
              HessianMatrix.element[n+1][n+2]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n+1][n+3]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared-Theta.az);

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-0.5*(Theta.ax+Theta.by));
              HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared-Theta.ay);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.bz+4.0*Rk.x*Rk.x*Rk.y*Rk.z/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.by*Theta.bz+4.0*Rk.y*Rk.y*Rk.y*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared-Theta.bz);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.bz*Theta.bz+1.0+4.0*Rk.y*Rk.z*Rk.y*Rk.z/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.z*Rk.z)*InverseLamdaSquared-Theta.cz);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.bz*Theta.cz+4.0*Rk.y*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.y*Rk.z*InverseLamdaSquared);

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.az+4.0*Rk.x*Rk.x*Rk.x*Rk.z/(SQR(Rksq))-4.0*Rk.z*Rk.x*InverseLamdaSquared-Theta.az);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.az*Theta.az+1.0+4.0*Rk.x*Rk.z*Rk.x*Rk.z/(SQR(Rksq))-2.0*(Rk.z*Rk.z+Rk.x*Rk.x)*InverseLamdaSquared-Theta.cz);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.az*Theta.by+4.0*Rk.x*Rk.z*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n+1][n+3]+=f*(Theta.az*Theta.cz+4.0*Rk.x*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-4.0*Rk.x*Rk.z*InverseLamdaSquared);

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the first term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f*(Theta.ax*Theta.ax+2.0+4.0*Rk.x*Rk.x*Rk.x*Rk.x/(SQR(Rksq))-8.0*Rk.x*Rk.x*InverseLamdaSquared-2.0*Theta.ax);
              HessianMatrix.element[n][n+1]+=f*(Theta.ax*Theta.ay+4.0*Rk.x*Rk.x*Rk.x*Rk.y/(SQR(Rksq))-4.0*Rk.y*Rk.x*InverseLamdaSquared-Theta.ay);
              HessianMatrix.element[n][n+2]+=f*(Theta.ax*Theta.by+4.0*Rk.x*Rk.x*Rk.y*Rk.y/(SQR(Rksq)));
              HessianMatrix.element[n][n+3]+=f*(Theta.ax*Theta.cz+4.0*Rk.x*Rk.x*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+1][n+1]+=f*(Theta.ay*Theta.ay+1.0+4.0*Rk.x*Rk.y*Rk.x*Rk.y/(SQR(Rksq))-2.0*(Rk.y*Rk.y+Rk.x*Rk.x)*InverseLamdaSquared-Theta.by);
              HessianMatrix.element[n+1][n+2]+=f*(Theta.ay*Theta.by+4.0*Rk.x*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-4.0*Rk.x*Rk.y*InverseLamdaSquared);
              HessianMatrix.element[n+1][n+3]+=f*(Theta.ay*Theta.cz+4.0*Rk.x*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+2][n+2]+=f*(Theta.by*Theta.by+2.0+4.0*Rk.y*Rk.y*Rk.y*Rk.y/(SQR(Rksq))-8.0*Rk.y*Rk.y*InverseLamdaSquared-2.0*Theta.by);
              HessianMatrix.element[n+2][n+3]+=f*(Theta.by*Theta.cz+4.0*Rk.y*Rk.y*Rk.z*Rk.z/(SQR(Rksq)));

              HessianMatrix.element[n+3][n+3]+=f*(Theta.cz*Theta.cz+2.0+4.0*Rk.z*Rk.z*Rk.z*Rk.z/(SQR(Rksq))-8.0*Rk.z*Rk.z*InverseLamdaSquared-2.0*Theta.cz);
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

// Hessian: Strain - Strain (part II)
// ==================================
static inline void HessianAtomicCorrectionStrainStrainI(REAL_MATRIX HessianMatrix,REAL f1_I,REAL f2_I,VECTOR posA,VECTOR comA,VECTOR Rk,REAL_MATRIX3x3 Theta)
{
  int n;
  REAL temp1;
  VECTOR dI;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                         // xx xx
      HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x);   // xx yy
      HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x);   // xx zz
      HessianMatrix.element[n][n]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                     // yy yy
      HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y); // yy zz
      HessianMatrix.element[n][n]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                     // zz zz

      temp1=f2_I*dI.x*Rk.x;
      HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
      HessianMatrix.element[n][n]+=2.0*temp1*dI.y*Rk.y;                   // xx yy
      HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                   // xx zz
      temp1=f2_I*dI.y*Rk.y;
      HessianMatrix.element[n][n]+=temp1*dI.y*Rk.y;                 // yy yy
      HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                 // yy zz
      temp1=f2_I*(posA.z-comA.z)*Rk.z;
      HessianMatrix.element[n][n]+=temp1*dI.z*Rk.z;                 // zz zz

      HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
      HessianMatrix.element[n][n]+=f1_I*Rk.y*dI.y;                  // yy yy
      HessianMatrix.element[n][n]+=f1_I*Rk.z*dI.z;                  // zz zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                         // xx xx
          HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x);   // xx yy
          HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x);   // xx zz
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                     // yy yy
          HessianMatrix.element[n][n]+=2.0*(f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y); // yy zz
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                     // zz zz

          temp1=f2_I*dI.x*Rk.x;
          HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n]+=2.0*temp1*dI.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                   // xx zz
          temp1=f2_I*dI.y*Rk.y;
          HessianMatrix.element[n][n]+=temp1*dI.y*Rk.y;                 // yy yy
          HessianMatrix.element[n][n]+=2.0*temp1*dI.z*Rk.z;                 // yy zz
          temp1=f2_I*(posA.z-comA.z)*Rk.z;
          HessianMatrix.element[n][n]+=temp1*dI.z*Rk.z;                 // zz zz

          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
          HessianMatrix.element[n][n]+=f1_I*Rk.y*dI.y;                  // yy yy
          HessianMatrix.element[n][n]+=f1_I*Rk.z*dI.z;                  // zz zz
          break;
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                         // xx xx
          HessianMatrix.element[n][n+1]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;   // xx yy
          HessianMatrix.element[n][n+2]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;   // xx zz
          HessianMatrix.element[n+1][n+1]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                     // yy yy
          HessianMatrix.element[n+1][n+2]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y; // yy zz
          HessianMatrix.element[n+2][n+2]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                     // zz zz

          temp1=f2_I*dI.x*Rk.x;
          HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]+=temp1*dI.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+2]+=temp1*dI.z*Rk.z;                   // xx zz
          temp1=f2_I*dI.y*Rk.y;
          HessianMatrix.element[n+1][n+1]+=temp1*dI.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+1][n+2]+=temp1*dI.z*Rk.z;                 // yy zz
          temp1=f2_I*(posA.z-comA.z)*Rk.z;
          HessianMatrix.element[n+2][n+2]+=temp1*dI.z*Rk.z;                 // zz zz

          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
          HessianMatrix.element[n+1][n+1]+=f1_I*Rk.y*dI.y;                  // yy yy
          HessianMatrix.element[n+2][n+2]+=f1_I*Rk.z*dI.z;                  // zz zz
          break;
        case REGULAR:
          // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
          HessianMatrix.element[n][n+1]+=0.5*f1_I*Theta.ax*(dI.y*Rk.x+dI.x*Rk.y)+f1_I*Theta.ay*dI.x*Rk.x; // xx xy
          HessianMatrix.element[n][n+2]+=0.5*f1_I*Theta.ax*(dI.z*Rk.x+dI.x*Rk.z)+f1_I*Theta.az*dI.x*Rk.x; // xx xz
          HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
          HessianMatrix.element[n][n+4]+=0.5*f1_I*Theta.ax*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.x*Rk.x; // xx yz
          HessianMatrix.element[n][n+5]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

          HessianMatrix.element[n+1][n+1]+=f1_I*Theta.ay*(dI.y*Rk.x+dI.x*Rk.y);                                      // xy xy
          HessianMatrix.element[n+1][n+2]+=0.5*f1_I*(Theta.ay*(dI.z*Rk.x+dI.x*Rk.z)+Theta.az*(dI.y*Rk.x+dI.x*Rk.y)); // xy xz
          HessianMatrix.element[n+1][n+3]+=f1_I*Theta.ay*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.y*Rk.x+dI.x*Rk.y);          // xy yy
          HessianMatrix.element[n+1][n+4]+=0.5*f1_I*(Theta.ay*(dI.z*Rk.y+dI.y*Rk.z)+Theta.bz*(dI.y*Rk.x+dI.x*Rk.y)); // xy yz
          HessianMatrix.element[n+1][n+5]+=f1_I*Theta.ay*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.y*Rk.x+dI.x*Rk.y);          // xy zz

          HessianMatrix.element[n+2][n+2]+=f1_I*Theta.az*(dI.z*Rk.x+dI.x*Rk.z);                                      // xz xz
          HessianMatrix.element[n+2][n+3]+=f1_I*Theta.az*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.z*Rk.x+dI.x*Rk.z);          // xz yy
          HessianMatrix.element[n+2][n+4]+=0.5*f1_I*(Theta.az*(dI.z*Rk.y+dI.y*Rk.z)+Theta.bz*(dI.z*Rk.x+dI.x*Rk.z)); // xz yz
          HessianMatrix.element[n+2][n+5]+=f1_I*Theta.az*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.x+dI.x*Rk.z);          // xz zz

          HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
          HessianMatrix.element[n+3][n+4]+=0.5*f1_I*Theta.by*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.y*Rk.y; // yy yz
          HessianMatrix.element[n+3][n+5]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

          HessianMatrix.element[n+4][n+4]+=f1_I*Theta.bz*(dI.z*Rk.y+dI.y*Rk.z);                             // yz yz
          HessianMatrix.element[n+4][n+5]+=f1_I*Theta.bz*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.y+dI.y*Rk.z); // yz zz

          HessianMatrix.element[n+5][n+5]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz


          // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp1=f2_I*dI.x*Rk.x;
          HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y);   // xx xy
          HessianMatrix.element[n][n+2]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z);   // xx xz
          HessianMatrix.element[n][n+3]+=temp1*dI.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z);   // xx yz
          HessianMatrix.element[n][n+5]+=temp1*dI.z*Rk.z;                   // xx zz

          temp1=f2_I*0.5*(dI.x*Rk.y+dI.y*Rk.x);
          HessianMatrix.element[n+1][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y); // xy xy
          HessianMatrix.element[n+1][n+2]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z); // xy xz
          HessianMatrix.element[n+1][n+3]+=temp1*dI.y*Rk.y;                 // xy yy
          HessianMatrix.element[n+1][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // xy yz
          HessianMatrix.element[n+1][n+5]+=temp1*dI.z*Rk.z;                 // xy zz

          temp1=f2_I*0.5*(dI.x*Rk.z+dI.z*Rk.x);
          HessianMatrix.element[n+2][n+2]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z); // xz xz
          HessianMatrix.element[n+2][n+3]+=temp1*dI.y*Rk.y;                 // xz yy
          HessianMatrix.element[n+2][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // xz yz
          HessianMatrix.element[n+2][n+5]+=temp1*dI.z*Rk.z;                 // xz zz

          temp1=f2_I*dI.y*Rk.y;
          HessianMatrix.element[n+3][n+3]+=temp1*dI.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+3][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yy yz
          HessianMatrix.element[n+3][n+5]+=temp1*dI.z*Rk.z;                 // yy zz

          temp1=f2_I*0.5*(dI.y*Rk.z+dI.z*Rk.y);
          HessianMatrix.element[n+4][n+4]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yz yz
          HessianMatrix.element[n+4][n+5]+=temp1*dI.z*Rk.z;                 // yz zz

          temp1=f2_I*(posA.z-comA.z)*Rk.z; 
          HessianMatrix.element[n+5][n+5]+=temp1*dI.z*Rk.z;                 // zz zz


          // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
          HessianMatrix.element[n][n+1]+=0.5*f1_I*Rk.x*dI.y;                // xx xy
          HessianMatrix.element[n][n+2]+=0.5*f1_I*Rk.x*dI.z;                // xx xz

          HessianMatrix.element[n+1][n+1]+=0.25*f1_I*(Rk.y*dI.y+Rk.x*dI.x); // xy xy
          HessianMatrix.element[n+1][n+2]+=0.25*f1_I*Rk.y*dI.z;             // xy xz
          HessianMatrix.element[n+1][n+3]+=0.5*f1_I*Rk.x*dI.y;              // xy yy
          HessianMatrix.element[n+1][n+4]+=0.25*f1_I*Rk.x*dI.z;             // xy yz

          HessianMatrix.element[n+2][n+2]+=0.25*f1_I*(Rk.z*dI.z+Rk.x*dI.x); // xz xz
          HessianMatrix.element[n+2][n+4]+=0.25*f1_I*Rk.x*dI.y;             // xz yz
          HessianMatrix.element[n+2][n+5]+=0.5*f1_I*Rk.x*dI.z;              // xz zz

          HessianMatrix.element[n+3][n+3]+=f1_I*Rk.y*dI.y;                  // yy yy
          HessianMatrix.element[n+3][n+4]+=0.5*f1_I*Rk.y*dI.z;              // yy yz

          HessianMatrix.element[n+4][n+4]+=0.25*f1_I*(Rk.z*dI.z+Rk.y*dI.y); // yz yz
          HessianMatrix.element[n+4][n+5]+=0.5*f1_I*Rk.y*dI.z;              // yz zz

          HessianMatrix.element[n+5][n+5]+=f1_I*Rk.z*dI.z;                  // zz zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
          HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.y*Rk.x+Rk.x*dI.x*Theta.ay); // xx xy
          HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.z*Rk.x+Rk.x*dI.x*Theta.az); // xx xz
          HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
          HessianMatrix.element[n][n+4]+=f1_I*(Theta.ax*dI.z*Rk.y+Rk.x*dI.x*Theta.bz); // xx yz
          HessianMatrix.element[n][n+5]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

          HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.ay*Rk.x*dI.y+Rk.x*dI.y*Theta.ay); // xy xy
          HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.ay*Rk.x*dI.z+Rk.x*dI.y*Theta.az); // xy xz
          HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.ay*Rk.y*dI.y+Rk.x*dI.y*Theta.by); // xy yy
          HessianMatrix.element[n+1][n+4]+=f1_I*(Theta.ay*Rk.y*dI.z+Rk.x*dI.y*Theta.bz); // xy yz
          HessianMatrix.element[n+1][n+5]+=f1_I*(Theta.ay*Rk.z*dI.z+Rk.x*dI.y*Theta.cz); // xy zz

          HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.az*Rk.x*dI.z+Rk.x*dI.z*Theta.az); // xz xz
          HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.az*Rk.y*dI.y+Rk.x*dI.z*Theta.by); // xz yy
          HessianMatrix.element[n+2][n+4]+=f1_I*(Theta.az*Rk.y*dI.z+Rk.x*dI.z*Theta.bz); // xz yz
          HessianMatrix.element[n+2][n+5]+=f1_I*(Theta.az*Rk.z*dI.z+Rk.x*dI.z*Theta.cz); // xz zz

          HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
          HessianMatrix.element[n+3][n+4]+=f1_I*(Theta.by*Rk.y*dI.z+Rk.y*dI.y*Theta.bz); // yy yz
          HessianMatrix.element[n+3][n+5]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

          HessianMatrix.element[n+4][n+4]+=f1_I*(Theta.bz*Rk.y*dI.z+Rk.y*dI.z*Theta.bz); // yz yz
          HessianMatrix.element[n+4][n+5]+=f1_I*(Theta.bz*Rk.z*dI.z+Rk.y*dI.z*Theta.cz); // yz zz

          HessianMatrix.element[n+5][n+5]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


          // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp1=f2_I*Rk.x*dI.x;
          HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
          HessianMatrix.element[n][n+1]+=temp1*Rk.x*dI.y;   // xx xy
          HessianMatrix.element[n][n+2]+=temp1*Rk.x*dI.z;   // xx xz
          HessianMatrix.element[n][n+3]+=temp1*Rk.y*dI.y;   // xx yy
          HessianMatrix.element[n][n+4]+=temp1*Rk.y*dI.z;   // xx yz
          HessianMatrix.element[n][n+5]+=temp1*Rk.z*dI.z;   // xx zz

          temp1=f2_I*Rk.x*dI.y;
          HessianMatrix.element[n+1][n+1]+=temp1*Rk.x*dI.y; // xy xy
          HessianMatrix.element[n+1][n+2]+=temp1*Rk.x*dI.z; // xy xz
          HessianMatrix.element[n+1][n+3]+=temp1*Rk.y*dI.y; // xy yy
          HessianMatrix.element[n+1][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+1][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.x*dI.z;
          HessianMatrix.element[n+2][n+2]+=temp1*Rk.x*dI.z; // xy xz
          HessianMatrix.element[n+2][n+3]+=temp1*Rk.y*dI.y; // xy yy
          HessianMatrix.element[n+2][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+2][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.y*dI.y;
          HessianMatrix.element[n+3][n+3]+=temp1*Rk.y*dI.y; // xy yy
          HessianMatrix.element[n+3][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+3][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.y*dI.z;
          HessianMatrix.element[n+4][n+4]+=temp1*Rk.y*dI.z; // xy yz
          HessianMatrix.element[n+4][n+5]+=temp1*Rk.z*dI.z; // xy zz

          temp1=f2_I*Rk.z*dI.z;
          HessianMatrix.element[n+5][n+5]+=temp1*Rk.z*dI.z; // xy zz

          // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // ================================================================
          HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx
          HessianMatrix.element[n][n+1]+=f1_I*Rk.x*dI.y;   // xx xy
          HessianMatrix.element[n][n+2]+=f1_I*Rk.x*dI.z;   // xx xz

          HessianMatrix.element[n+1][n+3]+=f1_I*Rk.x*dI.y; // xy yy
          HessianMatrix.element[n+1][n+4]+=f1_I*Rk.x*dI.z; // xy yz

          HessianMatrix.element[n+2][n+5]+=f1_I*Rk.x*dI.z; // xz zz

          HessianMatrix.element[n+3][n+3]+=f1_I*Rk.y*dI.y; // yy yy
          HessianMatrix.element[n+3][n+4]+=f1_I*Rk.y*dI.z; // yy yz

          HessianMatrix.element[n+4][n+5]+=f1_I*Rk.y*dI.z; // yz zz

          HessianMatrix.element[n+5][n+5]+=f1_I*Rk.z*dI.z; // zz zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
              HessianMatrix.element[n][n+2]+=0.5*f1_I*Theta.ax*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.x*Rk.x; // xx yz
              HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

              HessianMatrix.element[n+1][n+1]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
              HessianMatrix.element[n+1][n+2]+=0.5*f1_I*Theta.by*(dI.z*Rk.y+dI.y*Rk.z)+f1_I*Theta.bz*dI.y*Rk.y; // yy yz
              HessianMatrix.element[n+1][n+3]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

              HessianMatrix.element[n+2][n+2]+=f1_I*Theta.bz*(dI.z*Rk.y+dI.y*Rk.z);                             // yz yz
              HessianMatrix.element[n+2][n+3]+=f1_I*Theta.bz*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.y+dI.y*Rk.z); // yz zz

              HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*dI.x*Rk.x;
              HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]+=temp1*dI.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+2]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z);   // xx yz
              HessianMatrix.element[n][n+3]+=temp1*dI.z*Rk.z;                   // xx zz

              temp1=f2_I*dI.y*Rk.y;
              HessianMatrix.element[n+1][n+1]+=temp1*dI.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+1][n+2]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yy yz
              HessianMatrix.element[n+1][n+3]+=temp1*dI.z*Rk.z;                 // yy zz

              temp1=f2_I*0.5*(dI.y*Rk.z+dI.z*Rk.y);
              HessianMatrix.element[n+2][n+2]+=temp1*0.5*(dI.z*Rk.y+dI.y*Rk.z); // yz yz
              HessianMatrix.element[n+2][n+3]+=temp1*dI.z*Rk.z;                 // yz zz

              temp1=f2_I*(posA.z-comA.z)*Rk.z; 
              HessianMatrix.element[n+3][n+3]+=temp1*dI.z*Rk.z;                 // zz zz


              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx

              HessianMatrix.element[n+1][n+1]+=f1_I*Rk.y*dI.y;                  // yy yy
              HessianMatrix.element[n+1][n+2]+=0.5*f1_I*Rk.y*dI.z;              // yy yz

              HessianMatrix.element[n+2][n+2]+=0.25*f1_I*(Rk.z*dI.z+Rk.y*dI.y); // yz yz
              HessianMatrix.element[n+2][n+3]+=0.5*f1_I*Rk.y*dI.z;              // yz zz

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z;                  // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Theta.ax*(dI.z*Rk.x+dI.x*Rk.z)+f1_I*Theta.az*dI.x*Rk.x; // xx xz
              HessianMatrix.element[n][n+2]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*Theta.az*(dI.z*Rk.x+dI.x*Rk.z);                                      // xz xz
              HessianMatrix.element[n+1][n+2]+=f1_I*Theta.az*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.z*Rk.x+dI.x*Rk.z);          // xz yy
              HessianMatrix.element[n+1][n+3]+=f1_I*Theta.az*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.z*Rk.x+dI.x*Rk.z);          // xz zz

              HessianMatrix.element[n+2][n+2]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

              HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*dI.x*Rk.x;
              HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z);   // xx xz
              HessianMatrix.element[n][n+2]+=temp1*dI.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*dI.z*Rk.z;                   // xx zz

              temp1=f2_I*0.5*(dI.x*Rk.z+dI.z*Rk.x);
              HessianMatrix.element[n+1][n+1]+=temp1*0.5*(dI.z*Rk.x+dI.x*Rk.z); // xz xz
              HessianMatrix.element[n+1][n+2]+=temp1*dI.y*Rk.y;                 // xz yy
              HessianMatrix.element[n+1][n+3]+=temp1*dI.z*Rk.z;                 // xz zz

              temp1=f2_I*dI.y*Rk.y;
              HessianMatrix.element[n+2][n+2]+=temp1*dI.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]+=temp1*dI.z*Rk.z;                 // yy zz

              temp1=f2_I*(posA.z-comA.z)*Rk.z; 
              HessianMatrix.element[n+3][n+3]+=temp1*dI.z*Rk.z;                 // zz zz


              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Rk.x*dI.z;                // xx xz

              HessianMatrix.element[n+1][n+1]+=0.25*f1_I*(Rk.z*dI.z+Rk.x*dI.x); // xz xz
              HessianMatrix.element[n+1][n+3]+=0.5*f1_I*Rk.x*dI.z;              // xz zz

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y;                  // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z;                  // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=2.0*f1_I*Theta.ax*dI.x*Rk.x;                                       // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Theta.ax*(dI.y*Rk.x+dI.x*Rk.y)+f1_I*Theta.ay*dI.x*Rk.x; // xx xy
              HessianMatrix.element[n][n+2]+=f1_I*Theta.ax*dI.y*Rk.y+f1_I*Theta.by*dI.x*Rk.x;                 // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*Theta.ax*dI.z*Rk.z+f1_I*Theta.cz*dI.x*Rk.x;                 // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*Theta.ay*(dI.y*Rk.x+dI.x*Rk.y);                                      // xy xy
              HessianMatrix.element[n+1][n+2]+=f1_I*Theta.ay*dI.y*Rk.y+0.5*f1_I*Theta.by*(dI.y*Rk.x+dI.x*Rk.y);          // xy yy
              HessianMatrix.element[n+1][n+3]+=f1_I*Theta.ay*dI.z*Rk.z+0.5*f1_I*Theta.cz*(dI.y*Rk.x+dI.x*Rk.y);          // xy zz

              HessianMatrix.element[n+2][n+2]+=2.0*f1_I*Theta.by*dI.y*Rk.y;                                     // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*Theta.by*dI.z*Rk.z+f1_I*Theta.cz*dI.y*Rk.y;                 // yy zz

              HessianMatrix.element[n+3][n+3]+=2.0*f1_I*Theta.cz*dI.z*Rk.z;                                     // zz zz

              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*dI.x*Rk.x;
              HessianMatrix.element[n][n]+=temp1*dI.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y);   // xx xy
              HessianMatrix.element[n][n+2]+=temp1*dI.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*dI.z*Rk.z;                   // xx zz

              temp1=f2_I*0.5*(dI.x*Rk.y+dI.y*Rk.x);
              HessianMatrix.element[n+1][n+1]+=temp1*0.5*(dI.y*Rk.x+dI.x*Rk.y); // xy xy
              HessianMatrix.element[n+1][n+2]+=temp1*dI.y*Rk.y;                 // xy yy
              HessianMatrix.element[n+1][n+3]+=temp1*dI.z*Rk.z;                 // xy zz

              temp1=f2_I*dI.y*Rk.y;
              HessianMatrix.element[n+2][n+2]+=temp1*dI.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]+=temp1*dI.z*Rk.z;                 // yy zz

              temp1=f2_I*(posA.z-comA.z)*Rk.z; 
              HessianMatrix.element[n+3][n+3]+=temp1*dI.z*Rk.z;                 // zz zz


              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;                      // xx xx
              HessianMatrix.element[n][n+1]+=0.5*f1_I*Rk.x*dI.y;                // xx xy

              HessianMatrix.element[n+1][n+1]+=0.25*f1_I*(Rk.y*dI.y+Rk.x*dI.x); // xy xy
              HessianMatrix.element[n+1][n+2]+=0.5*f1_I*Rk.x*dI.y;              // xy yy

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y;                  // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z;                  // zz zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
              HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.z*Rk.y+Rk.x*dI.x*Theta.bz); // xx yz
              HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
              HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.by*Rk.y*dI.z+Rk.y*dI.y*Theta.bz); // yy yz
              HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

              HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.bz*Rk.y*dI.z+Rk.y*dI.z*Theta.bz); // yz yz
              HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.bz*Rk.z*dI.z+Rk.y*dI.z*Theta.cz); // yz zz

              HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
              HessianMatrix.element[n][n+1]+=temp1*Rk.y*dI.y;   // xx yy
              HessianMatrix.element[n][n+2]+=temp1*Rk.y*dI.z;   // xx yz
              HessianMatrix.element[n][n+3]+=temp1*Rk.z*dI.z;   // xx zz

              temp1=f2_I*Rk.y*dI.y;
              HessianMatrix.element[n+1][n+1]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+1][n+2]+=temp1*Rk.y*dI.z; // xy yz
              HessianMatrix.element[n+1][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.y*dI.z;
              HessianMatrix.element[n+2][n+2]+=temp1*Rk.y*dI.z; // xy yz
              HessianMatrix.element[n+2][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]+=temp1*Rk.z*dI.z; // xy zz

              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx

              HessianMatrix.element[n+1][n+1]+=f1_I*Rk.y*dI.y; // yy yy
              HessianMatrix.element[n+1][n+2]+=f1_I*Rk.y*dI.z; // yy yz

              HessianMatrix.element[n+2][n+3]+=f1_I*Rk.y*dI.z; // yz zz

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z; // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.z*Rk.x+Rk.x*dI.x*Theta.az); // xx xz
              HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.az*Rk.x*dI.z+Rk.x*dI.z*Theta.az); // xz xz
              HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.az*Rk.y*dI.y+Rk.x*dI.z*Theta.by); // xz yy
              HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.az*Rk.z*dI.z+Rk.x*dI.z*Theta.cz); // xz zz

              HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

              HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
              HessianMatrix.element[n][n+1]+=temp1*Rk.x*dI.z;   // xx xz
              HessianMatrix.element[n][n+2]+=temp1*Rk.y*dI.y;   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*Rk.z*dI.z;   // xx zz

              temp1=f2_I*Rk.x*dI.z;
              HessianMatrix.element[n+1][n+1]+=temp1*Rk.x*dI.z; // xy xz
              HessianMatrix.element[n+1][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+1][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+2][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]+=temp1*Rk.z*dI.z; // xy zz

              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*Rk.x*dI.z;   // xx xz

              HessianMatrix.element[n+1][n+3]+=f1_I*Rk.x*dI.z; // xz zz

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y; // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z; // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the second term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n  ]+=f1_I*(Theta.ax*dI.x*Rk.x+Rk.x*dI.x*Theta.ax); // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*(Theta.ax*dI.y*Rk.x+Rk.x*dI.x*Theta.ay); // xx xy
              HessianMatrix.element[n][n+2]+=f1_I*(Theta.ax*dI.y*Rk.y+Rk.x*dI.x*Theta.by); // xx yy
              HessianMatrix.element[n][n+3]+=f1_I*(Theta.ax*dI.z*Rk.z+Rk.x*dI.x*Theta.cz); // xx zz

              HessianMatrix.element[n+1][n+1]+=f1_I*(Theta.ay*Rk.x*dI.y+Rk.x*dI.y*Theta.ay); // xy xy
              HessianMatrix.element[n+1][n+2]+=f1_I*(Theta.ay*Rk.y*dI.y+Rk.x*dI.y*Theta.by); // xy yy
              HessianMatrix.element[n+1][n+3]+=f1_I*(Theta.ay*Rk.z*dI.z+Rk.x*dI.y*Theta.cz); // xy zz

              HessianMatrix.element[n+2][n+2]+=f1_I*(Theta.by*Rk.y*dI.y+Rk.y*dI.y*Theta.by); // yy yy
              HessianMatrix.element[n+2][n+3]+=f1_I*(Theta.by*Rk.z*dI.z+Rk.y*dI.y*Theta.cz); // yy zz

              HessianMatrix.element[n+3][n+3]+=f1_I*(Theta.cz*Rk.z*dI.z+Rk.z*dI.z*Theta.cz); // zz zz


              // the third term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_I*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]+=temp1*Rk.x*dI.x;   // xx xx
              HessianMatrix.element[n][n+1]+=temp1*Rk.x*dI.y;   // xx xy
              HessianMatrix.element[n][n+2]+=temp1*Rk.y*dI.y;   // xx yy
              HessianMatrix.element[n][n+3]+=temp1*Rk.z*dI.z;   // xx zz

              temp1=f2_I*Rk.x*dI.y;
              HessianMatrix.element[n+1][n+1]+=temp1*Rk.x*dI.y; // xy xy
              HessianMatrix.element[n+1][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+1][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]+=temp1*Rk.y*dI.y; // xy yy
              HessianMatrix.element[n+2][n+3]+=temp1*Rk.z*dI.z; // xy zz

              temp1=f2_I*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]+=temp1*Rk.z*dI.z; // xy zz

              // the fifth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // ================================================================
              HessianMatrix.element[n][n]+=f1_I*Rk.x*dI.x;     // xx xx
              HessianMatrix.element[n][n+1]+=f1_I*Rk.x*dI.y;   // xx xy

              HessianMatrix.element[n+1][n+2]+=f1_I*Rk.x*dI.y; // xy yy

              HessianMatrix.element[n+2][n+2]+=f1_I*Rk.y*dI.y; // yy yy

              HessianMatrix.element[n+3][n+3]+=f1_I*Rk.z*dI.z; // zz zz
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

// Hessian: Strain - Strain (part III)
// ===================================
static inline void HessianAtomicCorrectionStrainStrainJ(REAL_MATRIX HessianMatrix,REAL f2_IJ,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,VECTOR Rk)
{
  int n;
  REAL temp1;
  VECTOR dI,dJ;

  n=NumberOfCoordinatesMinimizationVariables;

  dI.x=posA.x-comA.x;
  dI.y=posA.y-comA.y;
  dI.z=posA.z-comA.z;

  dJ.x=posB.x-comB.x;
  dJ.y=posB.y-comB.y;
  dJ.z=posB.z-comB.z;

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      temp1=f2_IJ*dI.x*Rk.x; // xx
      HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
      HessianMatrix.element[n][n]-=2.0*temp1*dJ.y*Rk.y;                   // xx yy
      HessianMatrix.element[n][n]-=2.0*temp1*dJ.z*Rk.z;                   // xx zz
      temp1=f2_IJ*dI.y*Rk.y; // yy
      HessianMatrix.element[n][n]-=temp1*dJ.y*Rk.y;                 // yy yy
      HessianMatrix.element[n][n]-=2.0*temp1*dJ.z*Rk.z;                 // yy zz
      temp1=f2_IJ*dI.z*Rk.z; // zz
      HessianMatrix.element[n][n]-=temp1*dJ.z*Rk.z;                 // zz zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          temp1=f2_IJ*dI.x*Rk.x; // xx
          HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]-=temp1*dJ.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+2]-=temp1*dJ.z*Rk.z;                   // xx zz
          temp1=f2_IJ*dI.y*Rk.y; // yy
          HessianMatrix.element[n+1][n+1]-=temp1*dJ.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+1][n+2]-=temp1*dJ.z*Rk.z;                 // yy zz
          temp1=f2_IJ*dI.z*Rk.z; // zz
          HessianMatrix.element[n+2][n+2]-=temp1*dJ.z*Rk.z;                 // zz zz
          break;
        case REGULAR:
          // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================
          temp1=f2_IJ*dI.x*Rk.x; // xx
          HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
          HessianMatrix.element[n][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y);   // xx xy
          HessianMatrix.element[n][n+2]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z);   // xx xz
          HessianMatrix.element[n][n+3]-=temp1*dJ.y*Rk.y;                   // xx yy
          HessianMatrix.element[n][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z);   // xx yz
          HessianMatrix.element[n][n+5]-=temp1*dJ.z*Rk.z;                   // xx zz

          temp1=0.5*f2_IJ*(dI.x*Rk.y+dI.y*Rk.x); // xy
          HessianMatrix.element[n+1][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y); // xy xy
          HessianMatrix.element[n+1][n+2]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z); // xy xz
          HessianMatrix.element[n+1][n+3]-=temp1*dJ.y*Rk.y;                 // xy yy
          HessianMatrix.element[n+1][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // xy yz
          HessianMatrix.element[n+1][n+5]-=temp1*dJ.z*Rk.z;                 // xy zz

          temp1=0.5*f2_IJ*(dI.x*Rk.z+dI.z*Rk.x); // xz
          HessianMatrix.element[n+2][n+2]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z); // xz xz
          HessianMatrix.element[n+2][n+3]-=temp1*dJ.y*Rk.y;                 // xz yy
          HessianMatrix.element[n+2][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // xz yz
          HessianMatrix.element[n+2][n+5]-=temp1*dJ.z*Rk.z;                 // xz zz

          temp1=f2_IJ*dI.y*Rk.y; // yy
          HessianMatrix.element[n+3][n+3]-=temp1*dJ.y*Rk.y;                 // yy yy
          HessianMatrix.element[n+3][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yy yz
          HessianMatrix.element[n+3][n+5]-=temp1*dJ.z*Rk.z;                 // yy zz

          temp1=0.5*f2_IJ*(dI.y*Rk.z+dI.z*Rk.y); // yz
          HessianMatrix.element[n+4][n+4]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yz yz
          HessianMatrix.element[n+4][n+5]-=temp1*dJ.z*Rk.z;                 // yz zz

          temp1=f2_IJ*dI.z*Rk.z; // zz
          HessianMatrix.element[n+5][n+5]-=temp1*dJ.z*Rk.z;                 // zz zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
          // =================================================================

          temp1=f2_IJ*Rk.x*dI.x; 
          HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
          HessianMatrix.element[n][n+1]-=temp1*Rk.x*dJ.y;   // xx xy
          HessianMatrix.element[n][n+2]-=temp1*Rk.x*dJ.z;   // xx xz
          HessianMatrix.element[n][n+3]-=temp1*Rk.y*dJ.y;   // xx yy
          HessianMatrix.element[n][n+4]-=temp1*Rk.y*dJ.z;   // xx yz
          HessianMatrix.element[n][n+5]-=temp1*Rk.z*dJ.z;   // xx zz

          temp1=f2_IJ*Rk.x*dI.y;
          HessianMatrix.element[n+1][n+1]-=temp1*Rk.x*dJ.y; // xy xy
          HessianMatrix.element[n+1][n+2]-=temp1*Rk.x*dJ.z; // xy xz
          HessianMatrix.element[n+1][n+3]-=temp1*Rk.y*dJ.y; // xy yy
          HessianMatrix.element[n+1][n+4]-=temp1*Rk.y*dJ.z; // xy yz
          HessianMatrix.element[n+1][n+5]-=temp1*Rk.z*dJ.z; // xy zz

          temp1=f2_IJ*Rk.x*dI.z;
          HessianMatrix.element[n+2][n+2]-=temp1*Rk.x*dJ.z; // xz xz
          HessianMatrix.element[n+2][n+3]-=temp1*Rk.y*dJ.y; // xz yy
          HessianMatrix.element[n+2][n+4]-=temp1*Rk.y*dJ.z; // xz yz
          HessianMatrix.element[n+2][n+5]-=temp1*Rk.z*dJ.z; // xz zz

          temp1=f2_IJ*Rk.y*dI.y;
          HessianMatrix.element[n+3][n+3]-=temp1*Rk.y*dJ.y; // yy yy
          HessianMatrix.element[n+3][n+4]-=temp1*Rk.y*dJ.z; // yy yz
          HessianMatrix.element[n+3][n+5]-=temp1*Rk.z*dJ.z; // yy zz

          temp1=f2_IJ*Rk.y*dI.z;
          HessianMatrix.element[n+4][n+4]-=temp1*Rk.y*dJ.z; // yz yz
          HessianMatrix.element[n+4][n+5]-=temp1*Rk.z*dJ.z; // yz zz

          temp1=f2_IJ*Rk.z*dI.z;
          HessianMatrix.element[n+5][n+5]-=temp1*Rk.z*dJ.z; // zz zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_IJ*dI.x*Rk.x; // xx
              HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]-=temp1*dJ.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+2]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z);   // xx yz
              HessianMatrix.element[n][n+3]-=temp1*dJ.z*Rk.z;                   // xx zz

              temp1=f2_IJ*dI.y*Rk.y; // yy
              HessianMatrix.element[n+1][n+1]-=temp1*dJ.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+1][n+2]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yy yz
              HessianMatrix.element[n+1][n+3]-=temp1*dJ.z*Rk.z;                 // yy zz

              temp1=0.5*f2_IJ*(dI.y*Rk.z+dI.z*Rk.y); // yz
              HessianMatrix.element[n+2][n+2]-=0.5*temp1*(dJ.z*Rk.y+dJ.y*Rk.z); // yz yz
              HessianMatrix.element[n+2][n+3]-=temp1*dJ.z*Rk.z;                 // yz zz

              temp1=f2_IJ*dI.z*Rk.z; // zz
              HessianMatrix.element[n+3][n+3]-=temp1*dJ.z*Rk.z;                 // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_IJ*dI.x*Rk.x; // xx
              HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z);   // xx xz
              HessianMatrix.element[n][n+2]-=temp1*dJ.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*dJ.z*Rk.z;                   // xx zz

              temp1=0.5*f2_IJ*(dI.x*Rk.z+dI.z*Rk.x); // xz
              HessianMatrix.element[n+1][n+1]-=0.5*temp1*(dJ.z*Rk.x+dJ.x*Rk.z); // xz xz
              HessianMatrix.element[n+1][n+2]-=temp1*dJ.y*Rk.y;                 // xz yy
              HessianMatrix.element[n+1][n+3]-=temp1*dJ.z*Rk.z;                 // xz zz

              temp1=f2_IJ*dI.y*Rk.y; // yy
              HessianMatrix.element[n+2][n+2]-=temp1*dJ.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*dJ.z*Rk.z;                 // yy zz

              temp1=f2_IJ*dI.z*Rk.z; // zz
              HessianMatrix.element[n+3][n+3]-=temp1*dJ.z*Rk.z;                 // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================
              temp1=f2_IJ*dI.x*Rk.x; // xx
              HessianMatrix.element[n][n]-=temp1*dJ.x*Rk.x;                     // xx xx
              HessianMatrix.element[n][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y);   // xx xy
              HessianMatrix.element[n][n+2]-=temp1*dJ.y*Rk.y;                   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*dJ.z*Rk.z;                   // xx zz

              temp1=0.5*f2_IJ*(dI.x*Rk.y+dI.y*Rk.x); // xy
              HessianMatrix.element[n+1][n+1]-=0.5*temp1*(dJ.y*Rk.x+dJ.x*Rk.y); // xy xy
              HessianMatrix.element[n+1][n+2]-=temp1*dJ.y*Rk.y;                 // xy yy
              HessianMatrix.element[n+1][n+3]-=temp1*dJ.z*Rk.z;                 // xy zz

              temp1=f2_IJ*dI.y*Rk.y; // yy
              HessianMatrix.element[n+2][n+2]-=temp1*dJ.y*Rk.y;                 // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*dJ.z*Rk.z;                 // yy zz

              temp1=f2_IJ*dI.z*Rk.z; // zz
              HessianMatrix.element[n+3][n+3]-=temp1*dJ.z*Rk.z;                 // zz zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================

              temp1=f2_IJ*Rk.x*dI.x; 
              HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
              HessianMatrix.element[n][n+1]-=temp1*Rk.y*dJ.y;   // xx yy
              HessianMatrix.element[n][n+2]-=temp1*Rk.y*dJ.z;   // xx yz
              HessianMatrix.element[n][n+3]-=temp1*Rk.z*dJ.z;   // xx zz

              temp1=f2_IJ*Rk.y*dI.y;
              HessianMatrix.element[n+1][n+1]-=temp1*Rk.y*dJ.y; // yy yy
              HessianMatrix.element[n+1][n+2]-=temp1*Rk.y*dJ.z; // yy yz
              HessianMatrix.element[n+1][n+3]-=temp1*Rk.z*dJ.z; // yy zz

              temp1=f2_IJ*Rk.y*dI.z;
              HessianMatrix.element[n+2][n+2]-=temp1*Rk.y*dJ.z; // yz yz
              HessianMatrix.element[n+2][n+3]-=temp1*Rk.z*dJ.z; // yz zz

              temp1=f2_IJ*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]-=temp1*Rk.z*dJ.z; // zz zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================

              temp1=f2_IJ*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
              HessianMatrix.element[n][n+1]-=temp1*Rk.x*dJ.z;   // xx xz
              HessianMatrix.element[n][n+2]-=temp1*Rk.y*dJ.y;   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*Rk.z*dJ.z;   // xx zz

              temp1=f2_IJ*Rk.x*dI.z;
              HessianMatrix.element[n+1][n+1]-=temp1*Rk.x*dJ.z; // xz xz
              HessianMatrix.element[n+1][n+2]-=temp1*Rk.y*dJ.y; // xz yy
              HessianMatrix.element[n+1][n+3]-=temp1*Rk.z*dJ.z; // xz zz

              temp1=f2_IJ*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]-=temp1*Rk.y*dJ.y; // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*Rk.z*dJ.z; // yy zz

              temp1=f2_IJ*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]-=temp1*Rk.z*dJ.z; // zz zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              // the fourth term of Eq. 46 of Ref. Dubbeldam, Krishna, Snurr, 2009
              // =================================================================

              temp1=f2_IJ*Rk.x*dI.x;
              HessianMatrix.element[n][n  ]-=temp1*Rk.x*dJ.x;   // xx xx
              HessianMatrix.element[n][n+1]-=temp1*Rk.x*dJ.y;   // xx xy
              HessianMatrix.element[n][n+2]-=temp1*Rk.y*dJ.y;   // xx yy
              HessianMatrix.element[n][n+3]-=temp1*Rk.z*dJ.z;   // xx zz

              temp1=f2_IJ*Rk.x*dI.y;
              HessianMatrix.element[n+1][n+1]-=temp1*Rk.x*dJ.y; // xy xy
              HessianMatrix.element[n+1][n+2]-=temp1*Rk.y*dJ.y; // xy yy
              HessianMatrix.element[n+1][n+3]-=temp1*Rk.z*dJ.z; // xy zz

              temp1=f2_IJ*Rk.y*dI.y;
              HessianMatrix.element[n+2][n+2]-=temp1*Rk.y*dJ.y; // yy yy
              HessianMatrix.element[n+2][n+3]-=temp1*Rk.z*dJ.z; // yy zz

              temp1=f2_IJ*Rk.z*dI.z;
              HessianMatrix.element[n+3][n+3]-=temp1*Rk.z*dJ.z; // zz zz
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}


static inline void GradientStrain(REAL fac,REAL *Gradient,REAL_MATRIX3x3 Theta)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]-=fac*(Theta.ax+Theta.by+Theta.cz); // xx+yy+zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]-=fac*Theta.ax; // xx
          Gradient[n+1]-=fac*Theta.by; // yy
          Gradient[n+2]-=fac*Theta.cz; // zz
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]-=fac*Theta.ax; // xx
          Gradient[n+1]-=fac*Theta.bx; // xy
          Gradient[n+2]-=fac*Theta.cx; // xz
          Gradient[n+3]-=fac*Theta.by; // yy
          Gradient[n+4]-=fac*Theta.cy; // yz
          Gradient[n+5]-=fac*Theta.cz; // zz
          break;
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]-=fac*Theta.ax; // xx
              Gradient[n+1]-=fac*Theta.by; // yy
              Gradient[n+2]-=fac*Theta.cy; // yz
              Gradient[n+3]-=fac*Theta.cz; // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]-=fac*Theta.ax; // xx
              Gradient[n+1]-=fac*Theta.cx; // xz
              Gradient[n+2]-=fac*Theta.by; // yy
              Gradient[n+3]-=fac*Theta.cz; // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]-=fac*Theta.ax; // xx
              Gradient[n+1]-=fac*Theta.bx; // xy
              Gradient[n+2]-=fac*Theta.by; // yy
              Gradient[n+3]-=fac*Theta.cz; // zz
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

static inline void GradientStrainI(REAL *Gradient,REAL fac,VECTOR Rk,VECTOR posA,VECTOR comA)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Gradient[n]-=fac*Rk.x*(posA.x-comA.x)+fac*Rk.y*(posA.y-comA.y)+fac*Rk.z*(posA.z-comA.z); // xx+yy+zz
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
          Gradient[n+1]-=fac*Rk.y*(posA.y-comA.y); // yy
          Gradient[n+2]-=fac*Rk.z*(posA.z-comA.z); // zz
          break;
        case REGULAR:
          Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
          Gradient[n+1]-=0.5*fac*(Rk.x*(posA.y-comA.y)+Rk.y*(posA.x-comA.x)); // xy
          Gradient[n+2]-=0.5*fac*(Rk.x*(posA.z-comA.z)+Rk.z*(posA.x-comA.x)); // xz
          Gradient[n+3]-=fac*Rk.y*(posA.y-comA.y);                            // yy
          Gradient[n+4]-=0.5*fac*(Rk.y*(posA.z-comA.z)+Rk.z*(posA.y-comA.y)); // yz
          Gradient[n+5]-=fac*Rk.z*(posA.z-comA.z);                            // zz
          break;
        case REGULAR_UPPER_TRIANGLE:
          Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
          Gradient[n+1]-=fac*Rk.x*(posA.y-comA.y); // xy
          Gradient[n+2]-=fac*Rk.x*(posA.z-comA.z); // xz
          Gradient[n+3]-=fac*Rk.y*(posA.y-comA.y); // yy
          Gradient[n+4]-=fac*Rk.y*(posA.z-comA.z); // yz
          Gradient[n+5]-=fac*Rk.z*(posA.z-comA.z); // zz
          break;
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
              Gradient[n+1]-=fac*Rk.y*(posA.y-comA.y);                            // yy
              Gradient[n+2]-=0.5*fac*(Rk.y*(posA.z-comA.z)+Rk.z*(posA.y-comA.y)); // yz
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z);                            // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
              Gradient[n+1]-=0.5*fac*(Rk.x*(posA.z-comA.z)+Rk.z*(posA.x-comA.x)); // xz
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y);                            // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z);                            // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x);                            // xx
              Gradient[n+1]-=0.5*fac*(Rk.x*(posA.y-comA.y)+Rk.y*(posA.x-comA.x)); // xy
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y);                            // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z);                            // zz
              break;
          }
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
              Gradient[n+1]-=fac*Rk.y*(posA.y-comA.y); // yy
              Gradient[n+2]-=fac*Rk.y*(posA.z-comA.z); // yz
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z); // zz
              break;
            case MONOCLINIC_BETA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
              Gradient[n+1]-=fac*Rk.x*(posA.z-comA.z); // xz
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y); // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z); // zz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              Gradient[n  ]-=fac*Rk.x*(posA.x-comA.x); // xx
              Gradient[n+1]-=fac*Rk.x*(posA.y-comA.y); // xy
              Gradient[n+2]-=fac*Rk.y*(posA.y-comA.y); // yy
              Gradient[n+3]-=fac*Rk.z*(posA.z-comA.z); // zz
              break;
          }
          break;
        default:
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

static inline void HessianAtomicPositionStrainExcluded(REAL_MATRIX HessianMatrix,int index_i,int index_j,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      if(index_i>=0)
      {
        HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
        HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
        HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
      }

      if(index_j>=0)
      {
        HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x+f2*dr.y*dr.y*dr.x+f2*dr.z*dr.z*dr.x;
        HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y+f2*dr.y*dr.y*dr.y+2.0*f1*dr.y+f2*dr.z*dr.z*dr.y;
        HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z+f2*dr.y*dr.y*dr.z+f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
      }
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          if(index_i>=0)
          {
            HessianMatrix.element[index_i][n]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
            HessianMatrix.element[index_i][n+1]+=f2*dr.y*dr.y*dr.x;
            HessianMatrix.element[index_i][n+2]+=f2*dr.z*dr.z*dr.x;

            HessianMatrix.element[index_i+1][n]+=f2*dr.x*dr.x*dr.y;
            HessianMatrix.element[index_i+1][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            HessianMatrix.element[index_i+1][n+2]+=f2*dr.z*dr.z*dr.y;

            HessianMatrix.element[index_i+2][n]+=f2*dr.x*dr.x*dr.z;
            HessianMatrix.element[index_i+2][n+1]+=f2*dr.y*dr.y*dr.z;
            HessianMatrix.element[index_i+2][n+2]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j][n]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x;
            HessianMatrix.element[index_j][n+1]-=f2*dr.y*dr.y*dr.x;
            HessianMatrix.element[index_j][n+2]-=f2*dr.z*dr.z*dr.x;

            HessianMatrix.element[index_j+1][n]-=f2*dr.x*dr.x*dr.y;
            HessianMatrix.element[index_j+1][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y;
            HessianMatrix.element[index_j+1][n+2]-=f2*dr.z*dr.z*dr.y;

            HessianMatrix.element[index_j+2][n]-=f2*dr.x*dr.x*dr.z;
            HessianMatrix.element[index_j+2][n+1]-=f2*dr.y*dr.y*dr.z;
            HessianMatrix.element[index_j+2][n+2]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z;
          }
          break;
        case REGULAR:
        case REGULAR_UPPER_TRIANGLE:
          if(index_i>=0)
          {
            HessianMatrix.element[index_i  ][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
            HessianMatrix.element[index_i  ][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
            HessianMatrix.element[index_i  ][n+2]+=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
            HessianMatrix.element[index_i  ][n+3]+=f2*dr.y*dr.y*dr.x;             // yy x
            HessianMatrix.element[index_i  ][n+4]+=f2*dr.y*dr.z*dr.x;             // yz x
            HessianMatrix.element[index_i  ][n+5]+=f2*dr.z*dr.z*dr.x;             // zz x

            HessianMatrix.element[index_i+1][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
            HessianMatrix.element[index_i+1][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
            HessianMatrix.element[index_i+1][n+2]+=f2*dr.x*dr.z*dr.y;             // xz y
            HessianMatrix.element[index_i+1][n+3]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
            HessianMatrix.element[index_i+1][n+4]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
            HessianMatrix.element[index_i+1][n+5]+=f2*dr.z*dr.z*dr.y;             // zz y

            HessianMatrix.element[index_i+2][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
            HessianMatrix.element[index_i+2][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
            HessianMatrix.element[index_i+2][n+2]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
            HessianMatrix.element[index_i+2][n+3]+=f2*dr.y*dr.y*dr.z;             // yy z
            HessianMatrix.element[index_i+2][n+4]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
            HessianMatrix.element[index_i+2][n+5]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j  ][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
            HessianMatrix.element[index_j  ][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
            HessianMatrix.element[index_j  ][n+2]-=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
            HessianMatrix.element[index_j  ][n+3]-=f2*dr.y*dr.y*dr.x;             // yy x
            HessianMatrix.element[index_j  ][n+4]-=f2*dr.y*dr.z*dr.x;             // yz x
            HessianMatrix.element[index_j  ][n+5]-=f2*dr.z*dr.z*dr.x;             // zz x

            HessianMatrix.element[index_j+1][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
            HessianMatrix.element[index_j+1][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
            HessianMatrix.element[index_j+1][n+2]-=f2*dr.x*dr.z*dr.y;             // xz y
            HessianMatrix.element[index_j+1][n+3]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
            HessianMatrix.element[index_j+1][n+4]-=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
            HessianMatrix.element[index_j+1][n+5]-=f2*dr.z*dr.z*dr.y;             // zz y

            HessianMatrix.element[index_j+2][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
            HessianMatrix.element[index_j+2][n+1]-=f2*dr.x*dr.y*dr.z;             // xy z
            HessianMatrix.element[index_j+2][n+2]-=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
            HessianMatrix.element[index_j+2][n+3]-=f2*dr.y*dr.y*dr.z;             // yy z
            HessianMatrix.element[index_j+2][n+4]-=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
            HessianMatrix.element[index_j+2][n+5]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
          }
          break;                                   
        case MONOCLINIC:
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              if(index_i>=0)
              {
                HessianMatrix.element[index_i  ][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
                HessianMatrix.element[index_i  ][n+1]+=f2*dr.y*dr.y*dr.x;             // yy x
                HessianMatrix.element[index_i  ][n+2]+=f2*dr.y*dr.z*dr.x;             // yz x
                HessianMatrix.element[index_i  ][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

                HessianMatrix.element[index_i+1][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
                HessianMatrix.element[index_i+1][n+2]+=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
                HessianMatrix.element[index_i+1][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2*dr.y*dr.y*dr.z;             // yy z
                HessianMatrix.element[index_i+2][n+2]+=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
                HessianMatrix.element[index_i+2][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j  ][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
                HessianMatrix.element[index_j  ][n+1]-=f2*dr.y*dr.y*dr.x;             // yy x
                HessianMatrix.element[index_j  ][n+2]-=f2*dr.y*dr.z*dr.x;             // yz x
                HessianMatrix.element[index_j  ][n+3]-=f2*dr.z*dr.z*dr.x;             // zz x

                HessianMatrix.element[index_j+1][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
                HessianMatrix.element[index_j+1][n+1]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
                HessianMatrix.element[index_j+1][n+2]-=f2*dr.y*dr.z*dr.y+f1*dr.z;     // yz y
                HessianMatrix.element[index_j+1][n+3]-=f2*dr.z*dr.z*dr.y;             // zz y

                HessianMatrix.element[index_j+2][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
                HessianMatrix.element[index_j+2][n+1]-=f2*dr.y*dr.y*dr.z;             // yy z
                HessianMatrix.element[index_j+2][n+2]-=f2*dr.y*dr.z*dr.z+f1*dr.y;     // yz z
                HessianMatrix.element[index_j+2][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              }
              break;
            case MONOCLINIC_BETA_ANGLE:
              if(index_i>=0)
              {
                HessianMatrix.element[index_i  ][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
                HessianMatrix.element[index_i  ][n+1]+=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
                HessianMatrix.element[index_i  ][n+2]+=f2*dr.y*dr.y*dr.x;             // yy x
                HessianMatrix.element[index_i  ][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

                HessianMatrix.element[index_i+1][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2*dr.x*dr.z*dr.y;             // xz y
                HessianMatrix.element[index_i+1][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
                HessianMatrix.element[index_i+1][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
                HessianMatrix.element[index_i+2][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
                HessianMatrix.element[index_i+2][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j  ][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
                HessianMatrix.element[index_j  ][n+1]-=f2*dr.x*dr.z*dr.x+f1*dr.z;     // xz x
                HessianMatrix.element[index_j  ][n+2]-=f2*dr.y*dr.y*dr.x;             // yy x
                HessianMatrix.element[index_j  ][n+3]-=f2*dr.z*dr.z*dr.x;             // zz x

                HessianMatrix.element[index_j+1][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
                HessianMatrix.element[index_j+1][n+1]-=f2*dr.x*dr.z*dr.y;             // xz y
                HessianMatrix.element[index_j+1][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
                HessianMatrix.element[index_j+1][n+3]-=f2*dr.z*dr.z*dr.y;             // zz y

                HessianMatrix.element[index_j+2][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
                HessianMatrix.element[index_j+2][n+1]-=f2*dr.x*dr.z*dr.z+f1*dr.x;     // xz z
                HessianMatrix.element[index_j+2][n+2]-=f2*dr.y*dr.y*dr.z;             // yy z
                HessianMatrix.element[index_j+2][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              }
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              if(index_i>=0)
              {
                HessianMatrix.element[index_i  ][n  ]+=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
                HessianMatrix.element[index_i  ][n+1]+=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
                HessianMatrix.element[index_i  ][n+2]+=f2*dr.y*dr.y*dr.x;             // yy x
                HessianMatrix.element[index_i  ][n+3]+=f2*dr.z*dr.z*dr.x;             // zz x

                HessianMatrix.element[index_i+1][n  ]+=f2*dr.x*dr.x*dr.y;             // xx y
                HessianMatrix.element[index_i+1][n+1]+=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
                HessianMatrix.element[index_i+1][n+2]+=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
                HessianMatrix.element[index_i+1][n+3]+=f2*dr.z*dr.z*dr.y;             // zz y

                HessianMatrix.element[index_i+2][n  ]+=f2*dr.x*dr.x*dr.z;             // xx z
                HessianMatrix.element[index_i+2][n+1]+=f2*dr.x*dr.y*dr.z;             // xy z
                HessianMatrix.element[index_i+2][n+2]+=f2*dr.y*dr.y*dr.z;             // yy z
                HessianMatrix.element[index_i+2][n+3]+=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j  ][n  ]-=f2*dr.x*dr.x*dr.x+2.0*f1*dr.x; // xx x
                HessianMatrix.element[index_j  ][n+1]-=f2*dr.x*dr.y*dr.x+f1*dr.y;     // xy x
                HessianMatrix.element[index_j  ][n+2]-=f2*dr.y*dr.y*dr.x;             // yy x
                HessianMatrix.element[index_j  ][n+3]-=f2*dr.z*dr.z*dr.x;             // zz x

                HessianMatrix.element[index_j+1][n  ]-=f2*dr.x*dr.x*dr.y;             // xx y
                HessianMatrix.element[index_j+1][n+1]-=f2*dr.x*dr.y*dr.y+f1*dr.x;     // xy y
                HessianMatrix.element[index_j+1][n+2]-=f2*dr.y*dr.y*dr.y+2.0*f1*dr.y; // yy y
                HessianMatrix.element[index_j+1][n+3]-=f2*dr.z*dr.z*dr.y;             // zz y

                HessianMatrix.element[index_j+2][n  ]-=f2*dr.x*dr.x*dr.z;             // xx z
                HessianMatrix.element[index_j+2][n+1]-=f2*dr.x*dr.y*dr.z;             // xy z
                HessianMatrix.element[index_j+2][n+2]-=f2*dr.y*dr.y*dr.z;             // yy z
                HessianMatrix.element[index_j+2][n+3]-=f2*dr.z*dr.z*dr.z+2.0*f1*dr.z; // zz z
              }
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }
}

static inline void HessianAtomicStrainStrainExcluded(REAL_MATRIX HessianMatrix,int index_i,int index_j,REAL f1,REAL f2,VECTOR dr)
{
  int n;

  n=NumberOfCoordinatesMinimizationVariables;
  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x+2.0*f2*dr.x*dr.x*dr.y*dr.y+2.0*f2*dr.x*dr.x*dr.z*dr.z+
                            f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y+2.0*f2*dr.y*dr.y*dr.z*dr.z+f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
        case ANISOTROPIC:
          HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz
          HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
          HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz
          HessianMatrix.element[n+2][n+2]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
          break;
        case REGULAR:
          HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;    // xxxy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;    // xxxz
          HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
          HessianMatrix.element[n][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                 // xxyz
          HessianMatrix.element[n][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.x*dr.z+0.5*f1*dr.y*dr.z;              // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.y*dr.y+f1*dr.x*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dr.y*dr.y*dr.z+0.5*f1*(dr.x*dr.z);            // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dr.y*dr.z*dr.z;                               // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dr.z*dr.x*dr.z+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dr.z*dr.y*dr.y;                               // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dr.z*dr.y*dr.z+0.5*f1*dr.x*dr.y;              // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dr.z*dr.z*dr.z+f1*dr.x*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;                  // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dr.z*dr.y*dr.z+0.5*f1*(dr.y*dr.y+dr.z*dr.z);  // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dr.z*dr.z*dr.z+f1*dr.y*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
          break;                                   
        case REGULAR_UPPER_TRIANGLE:
          HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
          HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;       // xxxy
          HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;       // xxxz
          HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
          HessianMatrix.element[n][n+4]+=f2*dr.x*dr.x*dr.y*dr.z;                    // xxyz
          HessianMatrix.element[n][n+5]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz


          HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+f1*dr.y*dr.y;     // xyxy
          HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.x*dr.z+f1*dr.y*dr.z;     // xyxz
          HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.y*dr.y;                  // xyyy
          HessianMatrix.element[n+1][n+4]+=f2*dr.x*dr.y*dr.y*dr.z;                  // xyyz
          HessianMatrix.element[n+1][n+5]+=f2*dr.x*dr.y*dr.z*dr.z;                  // xyzz

          HessianMatrix.element[n+2][n+2]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
          HessianMatrix.element[n+2][n+3]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
          HessianMatrix.element[n+2][n+4]+=f2*dr.x*dr.z*dr.y*dr.z;                  // xzyz
          HessianMatrix.element[n+2][n+5]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

          HessianMatrix.element[n+3][n+3]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
          HessianMatrix.element[n+3][n+4]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
          HessianMatrix.element[n+3][n+5]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

          HessianMatrix.element[n+4][n+4]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
          HessianMatrix.element[n+4][n+5]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

          HessianMatrix.element[n+5][n+5]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
          break;                                   
        case MONOCLINIC:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                 // xxyz
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;                  // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+0.5*f1*(dr.y*dr.y+dr.z*dr.z);  // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z+f1*dr.y*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;    // xxxz
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+0.5*f1*(dr.x*dr.x+dr.z*dr.z);  // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                               // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z+f1*dr.x*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;  // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;    // xxxy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                 // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                 // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+0.5*f1*(dr.x*dr.x+dr.y*dr.y);  // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.y*dr.y+f1*dr.x*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.z*dr.z;                               // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y;              // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                               // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z;              // zzzz
              break;
          }
          break;  
        case MONOCLINIC_UPPER_TRIANGLE:
          switch(MonoclinicAngleType[CurrentSystem])
          {
            case MONOCLINIC_ALPHA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.z;                    // xxyz
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+1][n+2]+=f2*dr.y*dr.y*dr.y*dr.z+f1*dr.y*dr.z;     // yyyz
              HessianMatrix.element[n+1][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.z*dr.y*dr.z+f1*dr.z*dr.z;     // yzyz
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.z*dr.z*dr.z;                  // yzzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_BETA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.z+f1*dr.x*dr.z;       // xxxz
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz

              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.z*dr.x*dr.z+f1*dr.z*dr.z;     // xzxz
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.z*dr.y*dr.y;                  // xzyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.z*dr.z*dr.z;                  // xzzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
            case MONOCLINIC_GAMMA_ANGLE:
              HessianMatrix.element[n][n  ]+=f2*dr.x*dr.x*dr.x*dr.x+2.0*f1*dr.x*dr.x;   // xxxx
              HessianMatrix.element[n][n+1]+=f2*dr.x*dr.x*dr.x*dr.y+f1*dr.x*dr.y;       // xxxy
              HessianMatrix.element[n][n+2]+=f2*dr.x*dr.x*dr.y*dr.y;                    // xxyy
              HessianMatrix.element[n][n+3]+=f2*dr.x*dr.x*dr.z*dr.z;                    // xxzz


              HessianMatrix.element[n+1][n+1]+=f2*dr.x*dr.y*dr.x*dr.y+f1*dr.y*dr.y;     // xyxy
              HessianMatrix.element[n+1][n+2]+=f2*dr.x*dr.y*dr.y*dr.y;                  // xyyy
              HessianMatrix.element[n+1][n+3]+=f2*dr.x*dr.y*dr.z*dr.z;                  // xyzz

              HessianMatrix.element[n+2][n+2]+=f2*dr.y*dr.y*dr.y*dr.y+2.0*f1*dr.y*dr.y; // yyyy
              HessianMatrix.element[n+2][n+3]+=f2*dr.y*dr.y*dr.z*dr.z;                  // yyzz

              HessianMatrix.element[n+3][n+3]+=f2*dr.z*dr.z*dr.z*dr.z+2.0*f1*dr.z*dr.z; // zzzz
              break;
          }
          break;  
        default:  
          printf("Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
      break;
  }

}


void PrecomputeFixedWaveVectors(void)
{
  int i,j,m,n,f1;
  int I,ig,ia,A,TypeMolA;
  int Mmin,Nmin,Ll,Mm,Nn;
  int A1,A2,index,TotalNumberOfAtoms,Type,TypeMol;
  REAL Cs,Rksq,fac;
  REAL charge,temp,DipoleMagnitudeA,length,PreFactor4,PreFactor8;
  VECTOR Ss,pos,Rk,Rk1,Rk2,posA1,posA2;
  int NKVec;
  REAL USum,USumChargeChargeFramework,USumChargeChargeAdsorbates,USumChargeChargeCations;
  REAL USumChargeBondDipoleFramework,USumChargeBondDipoleAdsorbates,USumChargeBondDipoleCations;
  REAL USumBondDipoleBondDipoleFramework,USumBondDipoleBondDipoleAdsorbates,USumBondDipoleBondDipoleCations;
  REAL USumChargeChargeFrameworkAdsorbates,USumChargeChargeFrameworkCations,USumChargeChargeAdsorbatesCations;
  REAL USumChargeBondDipolesFrameworkAdsorbates,USumChargeBondDipolesFrameworkCations,USumChargeBondDipolesAdsorbatesCations;
  REAL USumBondDipolesBondDipolesFrameworkAdsorbates,USumBondDipolesBondDipolesFrameworkCations,USumBondDipolesBondDipolesAdsorbatesCations;
  REAL USelfSumChargesFramework,USelfSumChargesAdsorbates,USelfSumChargesCations;
  REAL USelfSumBondDipolesFramework,USelfSumBondDipolesAdsorbates,USelfSumBondDipolesCations;
  COMPLEX CksumFramework,CksumChargesFramework,CksumBondDipolesFramework,CksumChargesFrameworkFixed;
  COMPLEX CksumAdsorbates,CksumChargesAdsorbates,CksumBondDipolesAdsorbates,CksumChargesAdsorbatesFixed;
  COMPLEX CksumCations,CksumChargesCations,CksumBondDipolesCations,CksumChargesCationsFixed;
  VECTOR dipoleA;


  if((ChargeMethod!=EWALD)||OmitEwaldFourier) return;

  USum=USumChargeChargeFramework=USumChargeChargeAdsorbates=USumChargeChargeCations=0.0;
  USumChargeBondDipoleFramework=USumChargeBondDipoleAdsorbates=USumChargeBondDipoleCations=0.0;
  USumBondDipoleBondDipoleFramework=USumBondDipoleBondDipoleAdsorbates=USumBondDipoleBondDipoleCations=0.0;

  USumChargeChargeFrameworkAdsorbates=USumChargeChargeFrameworkCations=USumChargeChargeAdsorbatesCations=0.0;
  USumChargeBondDipolesFrameworkAdsorbates=USumChargeBondDipolesFrameworkCations=USumChargeBondDipolesAdsorbatesCations=0.0;
  USumBondDipolesBondDipolesFrameworkAdsorbates=USumBondDipolesBondDipolesFrameworkCations=USumBondDipolesBondDipolesAdsorbatesCations=0.0;

  PreFactor4=COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem]);
  PreFactor8=2.0*PreFactor4;

  ReciprocalCutOffSquared[CurrentSystem]=SQR(1.05*2.0*M_PI*
         MIN2((REAL)kvec[CurrentSystem].x*InverseBoxProperties[CurrentSystem].cx,
         MIN2((REAL)kvec[CurrentSystem].y*InverseBoxProperties[CurrentSystem].cy,
             (REAL)kvec[CurrentSystem].z*InverseBoxProperties[CurrentSystem].cz)));

  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UCationCationChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UHostHostChargeChargeFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostChargeBondDipoleFourier[CurrentSystem]=0.0;

  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=0.0;

  for(i=0;i<NumberOfKVectors[CurrentSystem];i++)
  {
    StoreCkChargeFramework[CurrentSystem][i].re=0.0;
    StoreCkChargeFramework[CurrentSystem][i].im=0.0;
    StoreCkChargeCations[CurrentSystem][i].re=0.0;
    StoreCkChargeCations[CurrentSystem][i].im=0.0;
    StoreCkChargeAdsorbates[CurrentSystem][i].re=0.0;
    StoreCkChargeAdsorbates[CurrentSystem][i].im=0.0;

    StoreCkBondDipolesFramework[CurrentSystem][i].re=0.0;
    StoreCkBondDipolesFramework[CurrentSystem][i].im=0.0;
    StoreCkBondDipolesCations[CurrentSystem][i].re=0.0;
    StoreCkBondDipolesCations[CurrentSystem][i].im=0.0;
    StoreCkBondDipolesAdsorbates[CurrentSystem][i].re=0.0;
    StoreCkBondDipolesAdsorbates[CurrentSystem][i].im=0.0;
  }

  fac=0.0;
  NKVec=0;
  Ss.x=Ss.y=Ss.z=0.0;
  Rk.x=Rk.y=Rk.z=0.0;
  Rk1.x=Rk1.y=Rk1.z=0.0;
  Rk2.x=Rk2.y=Rk2.z=0.0;

  USelfSumChargesFramework=USelfSumChargesAdsorbates=USelfSumChargesCations=0.0;
  USelfSumBondDipolesFramework=USelfSumBondDipolesAdsorbates=USelfSumBondDipolesCations=0.0;
  NetChargeFramework[CurrentSystem]=NetChargeCations[CurrentSystem]=NetChargeAdsorbates[CurrentSystem]=0.0;

  index=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Type=Framework[CurrentSystem].Atoms[f1][i].Type;
      if(PseudoAtoms[Type].HasCharges)
      {
        Charges[index]=PseudoAtoms[Type].Charge;
        USelfSumChargesFramework+=SQR(Charges[index])*Alpha[CurrentSystem]/sqrt(M_PI);
        NetChargeFramework[CurrentSystem]+=Charges[index];
        El[0][index].re=1.0; El[0][index].im=0.0;
        Em[0][index].re=1.0; Em[0][index].im=0.0;
        En[0][index].re=1.0; En[0][index].im=0.0;
        pos=Framework[CurrentSystem].Atoms[f1][i].Position;
        Ss=ConvertFromXYZtoABC(pos);
        Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
        El[1][index].re=cos(Ss.x); El[1][index].im=sin(Ss.x);
        Em[1][index].re=cos(Ss.y); Em[1][index].im=sin(Ss.y);
        En[1][index].re=cos(Ss.z); En[1][index].im=sin(Ss.z);
        index++;
      }
    }
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][i];
      A1=Framework[CurrentSystem].BondDipoles[f1][i].A;
      A2=Framework[CurrentSystem].BondDipoles[f1][i].B;
      posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
      posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      dipoleA=ApplyBoundaryCondition(dipoleA);
      pos.x=posA1.x+0.5*dipoleA.x;
      pos.y=posA1.y+0.5*dipoleA.y;
      pos.z=posA1.z+0.5*dipoleA.z;
      length=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      temp=DipoleMagnitudeA/length;
      BondDipoleMagnitude[index]=DipoleMagnitudeA;
      BondLength[index]=length;
      BondDipoleVector[index].x=temp*dipoleA.x;
      BondDipoleVector[index].y=temp*dipoleA.y;
      BondDipoleVector[index].z=temp*dipoleA.z;
      USelfSumBondDipolesFramework+=2.0*CUBE(Alpha[CurrentSystem])*
           (SQR(BondDipoleVector[index].x)+SQR(BondDipoleVector[index].y)+SQR(BondDipoleVector[index].z))/(3.0*sqrt(M_PI));
      BondDipolePosition[index]=pos;
      Ss=ConvertFromXYZtoABC(pos);
      Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
      El[1][index].re=cos(Ss.x); El[1][index].im=sin(Ss.x);
      Em[1][index].re=cos(Ss.y); Em[1][index].im=sin(Ss.y);
      En[1][index].re=cos(Ss.z); En[1][index].im=sin(Ss.z);
      El[0][index].re=1.0; El[0][index].im=0.0;
      Em[0][index].re=1.0; Em[0][index].im=0.0;
      En[0][index].re=1.0; En[0][index].im=0.0;
      index++;
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeMol=Adsorbates[CurrentSystem][i].Type;
    if(Components[TypeMol].HasCharges)
    {
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
        if(PseudoAtoms[Type].HasCharges)
        {
          Charges[index]=PseudoAtoms[Type].Charge;
          USelfSumChargesAdsorbates+=SQR(Charges[index])*Alpha[CurrentSystem]/sqrt(M_PI);
          NetChargeAdsorbates[CurrentSystem]+=Charges[index];
          El[0][index].re=1.0; El[0][index].im=0.0;
          Em[0][index].re=1.0; Em[0][index].im=0.0;
          En[0][index].re=1.0; En[0][index].im=0.0;
          pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;
          Ss=ConvertFromXYZtoABC(pos);
          Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
          El[1][index].re=cos(Ss.x); El[1][index].im=sin(Ss.x);
          Em[1][index].re=cos(Ss.y); Em[1][index].im=sin(Ss.y);
          En[1][index].re=cos(Ss.z); En[1][index].im=sin(Ss.z);
          index++;
        }
      }
    }
    for(j=0;j<Components[TypeMol].NumberOfBondDipoles;j++)
    {
      DipoleMagnitudeA=Components[TypeMol].BondDipoleMagnitude[j];
      A1=Components[TypeMol].BondDipoles[j].A;
      A2=Components[TypeMol].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      pos.x=posA1.x+0.5*dipoleA.x;
      pos.y=posA1.y+0.5*dipoleA.y;
      pos.z=posA1.z+0.5*dipoleA.z;
      length=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      temp=DipoleMagnitudeA/length;
      BondDipoleMagnitude[index]=DipoleMagnitudeA;
      BondLength[index]=length;
      BondDipoleVector[index].x=temp*dipoleA.x;
      BondDipoleVector[index].y=temp*dipoleA.y;
      BondDipoleVector[index].z=temp*dipoleA.z;
      USelfSumBondDipolesAdsorbates+=2.0*CUBE(Alpha[CurrentSystem])*
           (SQR(BondDipoleVector[index].x)+SQR(BondDipoleVector[index].y)+SQR(BondDipoleVector[index].z))/(3.0*sqrt(M_PI));
      BondDipolePosition[index]=pos;
      Ss=ConvertFromXYZtoABC(pos);
      Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
      El[1][index].re=cos(Ss.x); El[1][index].im=sin(Ss.x);
      Em[1][index].re=cos(Ss.y); Em[1][index].im=sin(Ss.y);
      En[1][index].re=cos(Ss.z); En[1][index].im=sin(Ss.z);
      El[0][index].re=1.0; El[0][index].im=0.0;
      Em[0][index].re=1.0; Em[0][index].im=0.0;
      En[0][index].re=1.0; En[0][index].im=0.0;
      index++;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeMol=Cations[CurrentSystem][i].Type;
    if(Components[TypeMol].HasCharges)
    {
      for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      {
        Type=Cations[CurrentSystem][i].Atoms[j].Type;
        if(PseudoAtoms[Type].HasCharges)
        {
          Charges[index]=PseudoAtoms[Type].Charge;
          USelfSumChargesCations+=SQR(Charges[index])*Alpha[CurrentSystem]/sqrt(M_PI);
          NetChargeCations[CurrentSystem]+=Charges[index];
          El[0][index].re=1.0; El[0][index].im=0.0;
          Em[0][index].re=1.0; Em[0][index].im=0.0;
          En[0][index].re=1.0; En[0][index].im=0.0;
          pos=Cations[CurrentSystem][i].Atoms[j].Position;
          Ss=ConvertFromXYZtoABC(pos);
          Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
          El[1][index].re=cos(Ss.x); El[1][index].im=sin(Ss.x);
          Em[1][index].re=cos(Ss.y); Em[1][index].im=sin(Ss.y);
          En[1][index].re=cos(Ss.z); En[1][index].im=sin(Ss.z);
          index++;
        }
      }
    }
    for(j=0;j<Components[TypeMol].NumberOfBondDipoles;j++)
    {
      DipoleMagnitudeA=Components[TypeMol].BondDipoleMagnitude[j];
      A1=Components[TypeMol].BondDipoles[j].A;
      A2=Components[TypeMol].BondDipoles[j].B;
      posA1=Cations[CurrentSystem][i].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][i].Atoms[A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      pos.x=posA1.x+0.5*dipoleA.x;
      pos.y=posA1.y+0.5*dipoleA.y;
      pos.z=posA1.z+0.5*dipoleA.z;
      length=sqrt(SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z));
      temp=DipoleMagnitudeA/length;
      BondDipoleMagnitude[index]=DipoleMagnitudeA;
      BondLength[index]=length;
      BondDipoleVector[index].x=temp*dipoleA.x;
      BondDipoleVector[index].y=temp*dipoleA.y;
      BondDipoleVector[index].z=temp*dipoleA.z;
      USelfSumBondDipolesCations+=2.0*CUBE(Alpha[CurrentSystem])*
           (SQR(BondDipoleVector[index].x)+SQR(BondDipoleVector[index].y)+SQR(BondDipoleVector[index].z))/(3.0*sqrt(M_PI));
      BondDipolePosition[index]=pos;
      Ss=ConvertFromXYZtoABC(pos);
      Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
      El[1][index].re=cos(Ss.x); El[1][index].im=sin(Ss.x);
      Em[1][index].re=cos(Ss.y); Em[1][index].im=sin(Ss.y);
      En[1][index].re=cos(Ss.z); En[1][index].im=sin(Ss.z);
      El[0][index].re=1.0; El[0][index].im=0.0;
      Em[0][index].re=1.0; Em[0][index].im=0.0;
      En[0][index].re=1.0; En[0][index].im=0.0;
      index++;
    }
  }

  TotalNumberOfAtoms=index;

  for(j=2;j<=kvec[CurrentSystem].y;j++)
    for(i=0;i<TotalNumberOfAtoms;i++)
    {
      Em[j][i].re=Em[j-1][i].re*Em[1][i].re-Em[j-1][i].im*Em[1][i].im;
      Em[j][i].im=Em[j-1][i].im*Em[1][i].re+Em[j-1][i].re*Em[1][i].im;
    }
  for(j=2;j<=kvec[CurrentSystem].z;j++)
    for(i=0;i<TotalNumberOfAtoms;i++)
    {
      En[j][i].re=En[j-1][i].re*En[1][i].re-En[j-1][i].im*En[1][i].im;
      En[j][i].im=En[j-1][i].im*En[1][i].re+En[j-1][i].re*En[1][i].im;
    }



  Mmin=0;
  Nmin=1;
  for(Ll=0;Ll<=kvec[CurrentSystem].x;Ll++)     // Loop Over All K-Vectors
  {
    Rk2.x=2.0*M_PI*Ll*InverseBox[CurrentSystem].ax;
    Rk2.y=2.0*M_PI*Ll*InverseBox[CurrentSystem].bx;
    Rk2.z=2.0*M_PI*Ll*InverseBox[CurrentSystem].cx;
    if(Ll==1)
      for(i=0;i<TotalNumberOfAtoms;i++)
      {
        El[0][i].re=El[1][i].re;
        El[0][i].im=El[1][i].im;
      }
    else if(Ll>1)
      for(i=0;i<TotalNumberOfAtoms;i++)
      {
        Cs=El[0][i].re;
        El[0][i].re=El[1][i].re*Cs-El[1][i].im*El[0][i].im;
        El[0][i].im=El[1][i].re*El[0][i].im+El[1][i].im*Cs;
      }

    for(Mm=Mmin;Mm<=kvec[CurrentSystem].y;Mm++)
    {
      Rk1.x=Rk2.x+2.0*M_PI*Mm*InverseBox[CurrentSystem].ay;
      Rk1.y=Rk2.y+2.0*M_PI*Mm*InverseBox[CurrentSystem].by;
      Rk1.z=Rk2.z+2.0*M_PI*Mm*InverseBox[CurrentSystem].cy;
      m=abs(Mm);
      if(Mm>=0)
        for(i=0;i<TotalNumberOfAtoms;i++)
        {
          lm[i].re=El[0][i].re*Em[m][i].re-El[0][i].im*Em[m][i].im;
          lm[i].im=El[0][i].im*Em[m][i].re+Em[m][i].im*El[0][i].re;
        }
      else
        for(i=0;i<TotalNumberOfAtoms;i++)
        {
          lm[i].re=El[0][i].re*Em[m][i].re+El[0][i].im*Em[m][i].im;
          lm[i].im=El[0][i].im*Em[m][i].re-Em[m][i].im*El[0][i].re;
        }

      for(Nn=Nmin;Nn<=kvec[CurrentSystem].z;Nn++)
      {
        n=abs(Nn);
        Rk.x=Rk1.x+2.0*M_PI*Nn*InverseBox[CurrentSystem].az;
        Rk.y=Rk1.y+2.0*M_PI*Nn*InverseBox[CurrentSystem].bz;
        Rk.z=Rk1.z+2.0*M_PI*Nn*InverseBox[CurrentSystem].cz;

        Rksq=SQR(Rk.x)+SQR(Rk.y)+SQR(Rk.z);

        //if(Rksq<ReciprocalCutOffSquared[CurrentSystem])
        {
          CksumFramework.re=CksumFramework.im=0.0;
          CksumChargesFramework.re=CksumChargesFramework.im=0.0;
          CksumBondDipolesFramework.re=CksumBondDipolesFramework.im=0.0;

          CksumChargesFrameworkFixed.re=CksumChargesFrameworkFixed.im=0.0;
          CksumChargesAdsorbatesFixed.re=CksumChargesAdsorbatesFixed.im=0.0;
          CksumChargesCationsFixed.re=CksumChargesCationsFixed.im=0.0;

          CksumAdsorbates.re=CksumAdsorbates.im=0.0;
          CksumChargesAdsorbates.re=CksumChargesAdsorbates.im=0.0;
          CksumBondDipolesAdsorbates.re=CksumBondDipolesAdsorbates.im=0.0;

          CksumCations.re=CksumCations.im=0.0;
          CksumChargesCations.re=CksumChargesCations.im=0.0;
          CksumBondDipolesCations.re=CksumBondDipolesCations.im=0.0;

          if(Nn>=0)
          {
            index=0;
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
              {
                Type=Framework[CurrentSystem].Atoms[f1][i].Type;
                if(PseudoAtoms[Type].HasCharges)
                {
                  charge=PseudoAtoms[Type].Charge;
                  Ck[index].re=lm[index].re*En[n][index].re-lm[index].im*En[n][index].im;
                  Ck[index].im=lm[index].im*En[n][index].re+lm[index].re*En[n][index].im;
                  CksumChargesFramework.re+=charge*Ck[index].re;
                  CksumChargesFramework.im+=charge*Ck[index].im;

                  if(Framework[CurrentSystem].Atoms[f1][i].Fixed)
                  {
                    CksumChargesFrameworkFixed.re+=charge*Ck[index].re;
                    CksumChargesFrameworkFixed.im+=charge*Ck[index].im;
                  }
                  index++;
                }
              }
            }

            for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
            {
              TypeMolA=Adsorbates[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if((Adsorbates[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Adsorbates[CurrentSystem][I].Groups[ig].FixedOrientation))
                      {
                        CksumChargesAdsorbatesFixed.re+=Charges[index]*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        CksumChargesAdsorbatesFixed.im+=Charges[index]*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                    charge=PseudoAtoms[Type].Charge;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if(Adsorbates[CurrentSystem][I].Atoms[A].Fixed)
                      {
                        CksumChargesAdsorbatesFixed.re+=charge*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        CksumChargesAdsorbatesFixed.im+=charge*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
              }
            }

            for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
            {
              TypeMolA=Cations[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Cations[CurrentSystem][I].Atoms[A].Type;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if((Cations[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Cations[CurrentSystem][I].Groups[ig].FixedOrientation))
                      {
                        CksumChargesCationsFixed.re+=Charges[index]*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        CksumChargesCationsFixed.im+=Charges[index]*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Cations[CurrentSystem][I].Atoms[A].Type;
                    charge=PseudoAtoms[Type].Charge;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if(Cations[CurrentSystem][I].Atoms[A].Fixed)
                      {
                        CksumChargesCationsFixed.re+=charge*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        CksumChargesCationsFixed.im+=charge*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
              }
            }


          }
          else
          {
            index=0;
            for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
            {
              for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
              {
                Type=Framework[CurrentSystem].Atoms[f1][i].Type;
                if(PseudoAtoms[Type].HasCharges)
                {
                  charge=PseudoAtoms[Type].Charge;
                  Ck[index].re=lm[index].re*En[n][index].re+lm[index].im*En[n][index].im;
                  Ck[index].im=lm[index].im*En[n][index].re-lm[index].re*En[n][index].im;
                  CksumChargesFramework.re+=charge*Ck[index].re;
                  CksumChargesFramework.im+=charge*Ck[index].im;

                  if(Framework[CurrentSystem].Atoms[f1][i].Fixed)
                  {
                    CksumChargesFrameworkFixed.re+=charge*Ck[index].re;
                    CksumChargesFrameworkFixed.im+=charge*Ck[index].im;
                  }
                  index++;
                }
              }
            }

            for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
            {
              TypeMolA=Adsorbates[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if((Adsorbates[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Adsorbates[CurrentSystem][I].Groups[ig].FixedOrientation))
                      {
                        CksumChargesAdsorbatesFixed.re+=Charges[index]*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                        CksumChargesAdsorbatesFixed.im+=Charges[index]*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                    charge=PseudoAtoms[Type].Charge;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if(Adsorbates[CurrentSystem][I].Atoms[A].Fixed)
                      {
                        CksumChargesAdsorbatesFixed.re+=charge*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                        CksumChargesAdsorbatesFixed.im+=charge*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
              }
            }

            for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
            {
              TypeMolA=Cations[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Cations[CurrentSystem][I].Atoms[A].Type;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if((Cations[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Cations[CurrentSystem][I].Groups[ig].FixedOrientation))
                      {
                        CksumChargesCationsFixed.re+=Charges[index]*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                        CksumChargesCationsFixed.im+=Charges[index]*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    Type=Cations[CurrentSystem][I].Atoms[A].Type;
                    charge=PseudoAtoms[Type].Charge;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      if(Cations[CurrentSystem][I].Atoms[A].Fixed)
                      {
                        CksumChargesCationsFixed.re+=charge*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                        CksumChargesCationsFixed.im+=charge*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                      }
                      index++;
                    }
                  }
                }
              }
            }



          }

          // Store the Ck sums 
          StoreCkChargeFramework[CurrentSystem][NKVec]=CksumChargesFrameworkFixed;
          StoreCkChargeAdsorbates[CurrentSystem][NKVec]=CksumChargesAdsorbatesFixed;
          StoreCkChargeCations[CurrentSystem][NKVec]=CksumChargesCationsFixed;

          NKVec++;
        }
      }
      Nmin=-kvec[CurrentSystem].z;
    }
    Mmin=-kvec[CurrentSystem].y;
  }
}


// Remark 10 February 2011: Check and fix atoms without charge
int CalculateEwaldFourierDerivatives(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian)
{
  int i,j,l,m,n,f1,f2;
  int Mmin,Nmin,Ll,Mm,Nn;
  int TotalNumberOfAtoms,Type,typeA,typeB;
  REAL Cs,Rksq,fac;
  REAL charge,exp_term,energy_term;
  VECTOR Ss,pos,Rk,Rk1,Rk2;
  int NKVec;
  REAL ChargeA,ChargeB;
  int A,B,NumberOfExcludedPairs,TypeA,TypeB;
  REAL r,rr;
  POINT posA,posB,comA,comB;
  VECTOR dr;
  int index,index_i,index_j,index_i2,index_j2;
  REAL_MATRIX3x3 Hessian,Stress;
  REAL USum,Uself_sum,DF,DDF;
  COMPLEX Cksum,CksumFramework;
  REAL InverseLamdaSquared;
  REAL Usum_framework;
  COMPLEX ctemp;
  int TypeMolA,TypeMolB;
  int I,J,ia,ja,ig,jg;
  int index1,index2;
  int index1_rigid,index2_rigid;
  REAL temp1,temp2,temp3;
  REAL_MATRIX3x3 S,Theta;
  int grpA,grpB;
  REAL f,f1_I,f2_I,f2_IJ;
  VECTOR dot_product_i,dot_product_j;
  REAL dot_product_AX,dot_product_BY,dot_product_CZ;
  REAL dot_product_AY,dot_product_AZ,dot_product_BZ;

  if(ChargeMethod==NONE) return 0;

  ReciprocalCutOffSquared[CurrentSystem]=SQR(1.05*2.0*M_PI*
     MIN2((REAL)kvec[CurrentSystem].x*InverseBoxProperties[CurrentSystem].cx,
     MIN2((REAL)kvec[CurrentSystem].y*InverseBoxProperties[CurrentSystem].cy,
         (REAL)kvec[CurrentSystem].z*InverseBoxProperties[CurrentSystem].cz)));

  fac=0.0;
  NKVec=0;
  USum=0.0;
  Uself_sum=0.0;
  Usum_framework=0.0;
  Ss.x=Ss.y=Ss.z=0.0;
  Rk.x=Rk.y=Rk.z=0.0;
  Rk1.x=Rk1.y=Rk1.z=0.0;
  Rk2.x=Rk2.y=Rk2.z=0.0;

  Stress.ax=Stress.bx=Stress.cx=0.0;
  Stress.ay=Stress.by=Stress.cy=0.0;
  Stress.az=Stress.bz=Stress.cz=0.0;

  index=0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Type=Framework[CurrentSystem].Atoms[f1][i].Type;
        if((!Framework[CurrentSystem].Atoms[f1][i].Fixed)&&(PseudoAtoms[Type].HasCharges))
        {
          Em[0][index].re=1.0;
          Em[0][index].im=0.0;
          El[0][index].re=1.0;
          El[0][index].im=0.0;
          En[0][index].re=1.0;
          En[0][index].im=0.0;
          pos=Framework[CurrentSystem].Atoms[f1][i].Position;
          Charges[index]=PseudoAtoms[Type].Charge;
          Ss=ConvertFromXYZtoABC(pos);
          Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
          El[1][index].re=cos(Ss.x);
          El[1][index].im=sin(Ss.x);
          Em[1][index].re=cos(Ss.y);
          Em[1][index].im=sin(Ss.y);
          En[1][index].re=cos(Ss.z);
          En[1][index].im=sin(Ss.z);
          index++;
          Uself_sum+=SQR(PseudoAtoms[Type].Charge);
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        if(!((Adsorbates[CurrentSystem][m].Groups[l].FixedCenterOfMass)&&(Adsorbates[CurrentSystem][m].Groups[l].FixedOrientation)))
        {
          for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
          {
            A=Components[Type].Groups[l].Atoms[i];
            if(PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].HasCharges)
            {
              pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;
              Em[0][index].re=1.0;
              Em[0][index].im=0.0;
              El[0][index].re=1.0;
              El[0][index].im=0.0;
              En[0][index].re=1.0;
              En[0][index].im=0.0;
              Charges[index]=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Charge;
              Ss=ConvertFromXYZtoABC(pos);
              Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
              El[1][index].re=cos(Ss.x);
              El[1][index].im=sin(Ss.x);
              Em[1][index].re=cos(Ss.y);
              Em[1][index].im=sin(Ss.y);
              En[1][index].re=cos(Ss.z);
              En[1][index].im=sin(Ss.z);
              index++;
            }
          }
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          if((PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].HasCharges)&&(!Adsorbates[CurrentSystem][m].Atoms[A].Fixed))
          {
            Em[0][index].re=1.0;
            Em[0][index].im=0.0;
            El[0][index].re=1.0;
            El[0][index].im=0.0;
            En[0][index].re=1.0;
            En[0][index].im=0.0;
            Charges[index]=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Charge;
            pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;
            Ss=ConvertFromXYZtoABC(pos);
            Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
            El[1][index].re=cos(Ss.x);
            El[1][index].im=sin(Ss.x);
            Em[1][index].re=cos(Ss.y);
            Em[1][index].im=sin(Ss.y);
            En[1][index].re=cos(Ss.z);
            En[1][index].im=sin(Ss.z);
            index++;
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        if(!((Cations[CurrentSystem][m].Groups[l].FixedCenterOfMass)&&(Cations[CurrentSystem][m].Groups[l].FixedOrientation)))
        {
          for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
          {
            A=Components[Type].Groups[l].Atoms[i];
            if(PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].HasCharges)
            {
              pos=Cations[CurrentSystem][m].Atoms[A].Position;
              Em[0][index].re=1.0;
              Em[0][index].im=0.0;
              El[0][index].re=1.0;
              El[0][index].im=0.0;
              En[0][index].re=1.0;
              En[0][index].im=0.0;
              Charges[index]=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Charge;
              Ss=ConvertFromXYZtoABC(pos);
              Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
              El[1][index].re=cos(Ss.x);
              El[1][index].im=sin(Ss.x);
              Em[1][index].re=cos(Ss.y);
              Em[1][index].im=sin(Ss.y);
              En[1][index].re=cos(Ss.z);
              En[1][index].im=sin(Ss.z);
              index++;
            }
          }
        }
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          if((PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].HasCharges)&&(!Cations[CurrentSystem][m].Atoms[A].Fixed))
          {
            Em[0][index].re=1.0;
            Em[0][index].im=0.0;
            El[0][index].re=1.0;
            El[0][index].im=0.0;
            En[0][index].re=1.0;
            En[0][index].im=0.0;
            Charges[index]=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Charge;
            pos=Cations[CurrentSystem][m].Atoms[A].Position;
            Ss=ConvertFromXYZtoABC(pos);
            Ss.x*=2.0*M_PI; Ss.y*=2.0*M_PI; Ss.z*=2.0*M_PI;
            El[1][index].re=cos(Ss.x);
            El[1][index].im=sin(Ss.x);
            Em[1][index].re=cos(Ss.y);
            Em[1][index].im=sin(Ss.y);
            En[1][index].re=cos(Ss.z);
            En[1][index].im=sin(Ss.z);
            index++;
          }
        }
      }
    }
  }

  TotalNumberOfAtoms=index;

  for(j=2;j<=kvec[CurrentSystem].y;j++)
    for(i=0;i<TotalNumberOfAtoms;i++)
    {
      Em[j][i].re=Em[j-1][i].re*Em[1][i].re-Em[j-1][i].im*Em[1][i].im;
      Em[j][i].im=Em[j-1][i].im*Em[1][i].re+Em[j-1][i].re*Em[1][i].im;
    }
  for(j=2;j<=kvec[CurrentSystem].z;j++)
    for(i=0;i<TotalNumberOfAtoms;i++)
    {
      En[j][i].re=En[j-1][i].re*En[1][i].re-En[j-1][i].im*En[1][i].im;
      En[j][i].im=En[j-1][i].im*En[1][i].re+En[j-1][i].re*En[1][i].im;
    }


  Mmin=0;
  Nmin=1;
  for(Ll=0;Ll<=kvec[CurrentSystem].x;Ll++)     // Loop Over All K-Vectors
  {
    Rk2.x=2.0*M_PI*Ll*InverseBox[CurrentSystem].ax;
    Rk2.y=2.0*M_PI*Ll*InverseBox[CurrentSystem].bx;
    Rk2.z=2.0*M_PI*Ll*InverseBox[CurrentSystem].cx;
    if(Ll==1)
      for(i=0;i<TotalNumberOfAtoms;i++)
      {
        El[0][i].re=El[1][i].re;
        El[0][i].im=El[1][i].im;
      }
    else if(Ll>1)
      for(i=0;i<TotalNumberOfAtoms;i++)
      {
        Cs=El[0][i].re;
        El[0][i].re=El[1][i].re*Cs-El[1][i].im*El[0][i].im;
        El[0][i].im=El[1][i].re*El[0][i].im+El[1][i].im*Cs;
      }

    for(Mm=Mmin;Mm<=kvec[CurrentSystem].y;Mm++)
    {
      Rk1.x=Rk2.x+2.0*M_PI*Mm*InverseBox[CurrentSystem].ay;
      Rk1.y=Rk2.y+2.0*M_PI*Mm*InverseBox[CurrentSystem].by;
      Rk1.z=Rk2.z+2.0*M_PI*Mm*InverseBox[CurrentSystem].cy;
      m=abs(Mm);
      if(Mm>=0)
      {
        for(i=0;i<TotalNumberOfAtoms;i++)
        {
          lm[i].re=El[0][i].re*Em[m][i].re-El[0][i].im*Em[m][i].im;
          lm[i].im=El[0][i].im*Em[m][i].re+Em[m][i].im*El[0][i].re;
        }
      }
      else
      {
        for(i=0;i<TotalNumberOfAtoms;i++)
        {
          lm[i].re=El[0][i].re*Em[m][i].re+El[0][i].im*Em[m][i].im;
          lm[i].im=El[0][i].im*Em[m][i].re-Em[m][i].im*El[0][i].re;
        }
      }

      for(Nn=Nmin;Nn<=kvec[CurrentSystem].z;Nn++)
      {
        n=abs(Nn);
        Rk.x=Rk1.x+2.0*M_PI*Nn*InverseBox[CurrentSystem].az;
        Rk.y=Rk1.y+2.0*M_PI*Nn*InverseBox[CurrentSystem].bz;
        Rk.z=Rk1.z+2.0*M_PI*Nn*InverseBox[CurrentSystem].cz;
        Rksq=SQR(Rk.x)+SQR(Rk.y)+SQR(Rk.z);

        //if(Rksq<ReciprocalCutOffSquared[CurrentSystem])
        {
          Cksum.re=0.0;
          Cksum.im=0.0;

          CksumFramework.re=0.0;
          CksumFramework.im=0.0;

          if(Nn>=0)
          {
            index=0;
            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
                {
                  if(!Framework[CurrentSystem].Atoms[f1][i].Fixed)
                  {
                    Type=Framework[CurrentSystem].Atoms[f1][i].Type;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      charge=PseudoAtoms[Type].Charge;
                      ctemp.re=charge*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                      ctemp.im=charge*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                      Ck[index].re=ctemp.re;
                      Ck[index].im=ctemp.im;
                      Cksum.re+=ctemp.re;
                      Cksum.im+=ctemp.im;
                      index++;
                    }
                  }
                }
              }
            }

            for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
            {
              TypeMolA=Adsorbates[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  if(!((Adsorbates[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Adsorbates[CurrentSystem][I].Groups[ig].FixedOrientation)))
                  {
                    for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                    {
                      A=Components[TypeMolA].Groups[ig].Atoms[ia];
                      Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                        ctemp.re=Charges[index]*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        ctemp.im=Charges[index]*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                        Ck[index].re=ctemp.re;
                        Ck[index].im=ctemp.im;
                        Cksum.re+=ctemp.re;
                        Cksum.im+=ctemp.im;
                        index++;
                      }
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    if(!Adsorbates[CurrentSystem][I].Atoms[A].Fixed)
                    {
                      Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                        charge=PseudoAtoms[Type].Charge;
                        ctemp.re=charge*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        ctemp.im=charge*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                        Ck[index].re=ctemp.re;
                        Ck[index].im=ctemp.im;
                        Cksum.re+=ctemp.re;
                        Cksum.im+=ctemp.im;
                        index++;
                      }
                    }
                  }
                }
              }
            }

            for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
            {
              TypeMolA=Cations[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  if(!((Cations[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Cations[CurrentSystem][I].Groups[ig].FixedOrientation)))
                  {
                    for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                    {
                      A=Components[TypeMolA].Groups[ig].Atoms[ia];
                      Type=Cations[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                        ctemp.re=Charges[index]*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        ctemp.im=Charges[index]*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                        Ck[index].re=ctemp.re;
                        Ck[index].im=ctemp.im;
                        Cksum.re+=ctemp.re;
                        Cksum.im+=ctemp.im;
                        index++;
                      }
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    if(!Cations[CurrentSystem][I].Atoms[A].Fixed)
                    {
                      Type=Cations[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                        charge=PseudoAtoms[Type].Charge;
                        ctemp.re=charge*(lm[index].re*En[n][index].re-lm[index].im*En[n][index].im);
                        ctemp.im=charge*(lm[index].im*En[n][index].re+lm[index].re*En[n][index].im);
                        Ck[index].re=ctemp.re;
                        Ck[index].im=ctemp.im;
                        Cksum.re+=ctemp.re;
                        Cksum.im+=ctemp.im;
                        index++;
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {
            index=0;
            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
                {
                  if(!Framework[CurrentSystem].Atoms[f1][i].Fixed)
                  {
                    Type=Framework[CurrentSystem].Atoms[f1][i].Type;
                    if(PseudoAtoms[Type].HasCharges)
                    {
                      charge=PseudoAtoms[Type].Charge;
                      ctemp.re=charge*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                      ctemp.im=charge*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                      Ck[index].re=ctemp.re;
                      Ck[index].im=ctemp.im;
                      Cksum.re+=ctemp.re;
                      Cksum.im+=ctemp.im;
                      index++;
                    }
                  }
                }
              }
            }

            for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
            {
              TypeMolA=Adsorbates[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  if(!((Adsorbates[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Adsorbates[CurrentSystem][I].Groups[ig].FixedOrientation)))
                  {
                    for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                    {
                      A=Components[TypeMolA].Groups[ig].Atoms[ia];
                      Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                        ctemp.re=Charges[index]*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                        ctemp.im=Charges[index]*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                        Ck[index].re=ctemp.re;
                        Ck[index].im=ctemp.im;
                        Cksum.re+=ctemp.re;
                        Cksum.im+=ctemp.im;
                        index++;
                      }
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    if(!Adsorbates[CurrentSystem][I].Atoms[A].Fixed)
                    {
                      Type=Adsorbates[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                      charge=PseudoAtoms[Type].Charge;
                      ctemp.re=charge*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                      ctemp.im=charge*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                      Ck[index].re=ctemp.re;
                      Ck[index].im=ctemp.im;
                      Cksum.re+=ctemp.re;
                      Cksum.im+=ctemp.im;
                      index++;

                      }
                    }
                  }
                }
              }
            }

            for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
            {
              TypeMolA=Cations[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                {
                  if(!((Cations[CurrentSystem][I].Groups[ig].FixedCenterOfMass)&&(Cations[CurrentSystem][I].Groups[ig].FixedOrientation)))
                  {
                    for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                    {
                      A=Components[TypeMolA].Groups[ig].Atoms[ia];
                      Type=Cations[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                        ctemp.re=Charges[index]*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                        ctemp.im=Charges[index]*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                        Ck[index].re=ctemp.re;
                        Ck[index].im=ctemp.im;
                        Cksum.re+=ctemp.re;
                        Cksum.im+=ctemp.im;
                        index++;
                      }
                    }
                  }
                }
                else
                {
                  for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                  {
                    A=Components[TypeMolA].Groups[ig].Atoms[ia];
                    if(!Cations[CurrentSystem][I].Atoms[A].Fixed)
                    {
                      Type=Cations[CurrentSystem][I].Atoms[A].Type;
                      if(PseudoAtoms[Type].HasCharges)
                      {
                        charge=PseudoAtoms[Type].Charge;
                        ctemp.re=charge*(lm[index].re*En[n][index].re+lm[index].im*En[n][index].im);
                        ctemp.im=charge*(lm[index].im*En[n][index].re-lm[index].re*En[n][index].im);
                        Ck[index].re=ctemp.re;
                        Ck[index].im=ctemp.im;
                        Cksum.re+=ctemp.re;
                        Cksum.im+=ctemp.im;
                        index++;
                      }
                    }
                  }
                }
              }
            }
          }

          Cksum.re+=StoreCkChargeFramework[CurrentSystem][NKVec].re;
          Cksum.im+=StoreCkChargeFramework[CurrentSystem][NKVec].im;

          Cksum.re+=StoreCkChargeAdsorbates[CurrentSystem][NKVec].re;
          Cksum.im+=StoreCkChargeAdsorbates[CurrentSystem][NKVec].im;

          Cksum.re+=StoreCkChargeCations[CurrentSystem][NKVec].re;
          Cksum.im+=StoreCkChargeCations[CurrentSystem][NKVec].im;

          exp_term=exp((-0.25/SQR(Alpha[CurrentSystem]))*Rksq)/Rksq;
          energy_term=(SQR(Cksum.re)+SQR(Cksum.im))*exp_term; 

          f=COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*energy_term;

          // energy sums
          USum+=f;

          // HERE
          if(Framework[CurrentSystem].FrameworkModel!=FLEXIBLE)
          {
            USum-=COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                  (SQR(StoreCkChargeFramework[CurrentSystem][NKVec].re)+SQR(StoreCkChargeFramework[CurrentSystem][NKVec].im));
            f-=COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                  (SQR(StoreCkChargeFramework[CurrentSystem][NKVec].re)+SQR(StoreCkChargeFramework[CurrentSystem][NKVec].im));
          }

          InverseLamdaSquared=0.25/(SQR(Alpha[CurrentSystem]))+1.0/Rksq;

          Theta.ax=1.0-2.0*Rk.x*Rk.x*InverseLamdaSquared;
          Theta.ay=-2.0*Rk.x*Rk.y*InverseLamdaSquared;
          Theta.az=-2.0*Rk.x*Rk.z*InverseLamdaSquared;

          Theta.bx=-2.0*Rk.y*Rk.x*InverseLamdaSquared;
          Theta.by=1.0-2.0*Rk.y*Rk.y*InverseLamdaSquared;
          Theta.bz=-2.0*Rk.y*Rk.z*InverseLamdaSquared;

          Theta.cx=-2.0*Rk.z*Rk.x*InverseLamdaSquared;
          Theta.cy=-2.0*Rk.z*Rk.y*InverseLamdaSquared;
          Theta.cz=1.0-2.0*Rk.z*Rk.z*InverseLamdaSquared;

          StrainFirstDerivative->ax-=f*Theta.ax;
          StrainFirstDerivative->bx-=f*Theta.bx;
          StrainFirstDerivative->cx-=f*Theta.cx;

          StrainFirstDerivative->ay-=f*Theta.ay;
          StrainFirstDerivative->by-=f*Theta.by;
          StrainFirstDerivative->cy-=f*Theta.cy;

          StrainFirstDerivative->az-=f*Theta.az;
          StrainFirstDerivative->bz-=f*Theta.bz;
          StrainFirstDerivative->cz-=f*Theta.cz;

          // Stress tensor
          index1=0;
          if(ComputeGradient)
          {
            GradientStrain(f,Gradient,Theta);

            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
                {
                  Type=Framework[CurrentSystem].Atoms[f1][i].Type;
                  if(PseudoAtoms[Type].HasCharges)
                  { 
                    index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

                    if(index_i>=0)
                    {
                      f1_I=2.0*COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                           (-Ck[index1].im*Cksum.re+Ck[index1].re*Cksum.im);
                      Gradient[index_i]+=f1_I*Rk.x;
                      Gradient[index_i+1]+=f1_I*Rk.y;
                      Gradient[index_i+2]+=f1_I*Rk.z;
                      index1++;
                    }
                  }
                }
              }
            }
          }

          for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
          {
            TypeMolA=Adsorbates[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              comA.x=comA.y=comA.z=0.0;
              if(Components[TypeMolA].Groups[ig].Rigid)
                comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

              for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
              {
                i=Components[TypeMolA].Groups[ig].Atoms[ia];

                if(PseudoAtoms[Adsorbates[CurrentSystem][I].Atoms[i].Type].HasCharges)
                {
                  f1_I=2.0*COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                       (-Ck[index1].im*Cksum.re+Ck[index1].re*Cksum.im);

                  posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
                    index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;

                    index1_rigid=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;
                  }
                  else
                  {
                    index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
                    index_i2=-1;
                    index1_rigid=-1;
                  }

                  if(ComputeGradient)
                  {
                    if(index_i>=0)
                    {
                      Gradient[index_i]+=f1_I*Rk.x;
                      Gradient[index_i+1]+=f1_I*Rk.y;
                      Gradient[index_i+2]+=f1_I*Rk.z;
                    }
                  }

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    if(ComputeGradient)
                    {
                      if(index_i2>=0)
                      {
                        Gradient[index_i2]+=f1_I*(Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z);
                        Gradient[index_i2+1]+=f1_I*(Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z);
                        Gradient[index_i2+2]+=f1_I*(Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z);
                      }

                      GradientStrainI(Gradient,f1_I,Rk,posA,comA);
                    }

                    pos=Components[TypeMolA].Positions[i];
                    temp1=0.5*f1_I*((posA.y-comA.y)*Rk.x+(posA.x-comA.x)*Rk.y);
                    temp2=0.5*f1_I*((posA.z-comA.z)*Rk.x+(posA.x-comA.x)*Rk.z);
                    temp3=0.5*f1_I*((posA.z-comA.z)*Rk.y+(posA.y-comA.y)*Rk.z);

                    StrainFirstDerivative->ax-=f1_I*(posA.x-comA.x)*Rk.x;
                    StrainFirstDerivative->bx-=temp1;
                    StrainFirstDerivative->cx-=temp2;
                    StrainFirstDerivative->ay-=temp1;
                    StrainFirstDerivative->by-=f1_I*(posA.y-comA.y)*Rk.y;
                    StrainFirstDerivative->cy-=temp3;
                    StrainFirstDerivative->az-=temp2;
                    StrainFirstDerivative->bz-=temp3;
                    StrainFirstDerivative->cz-=f1_I*(posA.z-comA.z)*Rk.z;
                  }

                  if((index_i>=0)||(index_i2>=0))
                    index1++;
                }
              }
            }
          }

          for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
          {
            TypeMolA=Cations[CurrentSystem][I].Type;
            for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
            {
              comA.x=comA.y=comA.z=0.0;
              if(Components[TypeMolA].Groups[ig].Rigid)
                comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

              for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
              {
                i=Components[TypeMolA].Groups[ig].Atoms[ia];

                if(PseudoAtoms[Cations[CurrentSystem][I].Atoms[i].Type].HasCharges)
                {
                  f1_I=2.0*COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                       (-Ck[index1].im*Cksum.re+Ck[index1].re*Cksum.im);

                  posA=Cations[CurrentSystem][I].Atoms[i].Position;

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
                    index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;

                    index1_rigid=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;
                  }
                  else
                  {
                    index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
                    index_i2=-1;
                    index1_rigid=-1;
                  }

                  if(ComputeGradient)
                  {
                    if(index_i>=0)
                    {
                      Gradient[index_i]+=f1_I*Rk.x;
                      Gradient[index_i+1]+=f1_I*Rk.y;
                      Gradient[index_i+2]+=f1_I*Rk.z;
                    }
                  }

                  if(Components[TypeMolA].Groups[ig].Rigid)
                  {
                    if(ComputeGradient)
                    {
                      if(index_i2>=0)
                      {
                        Gradient[index_i2]+=f1_I*(Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z);
                        Gradient[index_i2+1]+=f1_I*(Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z);
                        Gradient[index_i2+2]+=f1_I*(Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z);
                      }

                      GradientStrainI(Gradient,f1_I,Rk,posA,comA);
                    }

                    pos=Components[TypeMolA].Positions[i];
                    temp1=0.5*f1_I*((posA.y-comA.y)*Rk.x+(posA.x-comA.x)*Rk.y);
                    temp2=0.5*f1_I*((posA.z-comA.z)*Rk.x+(posA.x-comA.x)*Rk.z);
                    temp3=0.5*f1_I*((posA.z-comA.z)*Rk.y+(posA.y-comA.y)*Rk.z);

                    StrainFirstDerivative->ax-=f1_I*(posA.x-comA.x)*Rk.x;
                    StrainFirstDerivative->bx-=temp1;
                    StrainFirstDerivative->cx-=temp2;
                    StrainFirstDerivative->ay-=temp1;
                    StrainFirstDerivative->by-=f1_I*(posA.y-comA.y)*Rk.y;
                    StrainFirstDerivative->cy-=temp3;
                    StrainFirstDerivative->az-=temp2;
                    StrainFirstDerivative->bz-=temp3;
                    StrainFirstDerivative->cz-=f1_I*(posA.z-comA.z)*Rk.z;
                  }

                  if((index_i>=0)||(index_i2>=0))
                    index1++;
                }
              }
            }
          }


          if(ComputeHessian)
          {
            HessianAtomicStrainStrain(HessianMatrix,f,InverseLamdaSquared,Rk,Rksq,Theta);

            index1=0;
            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
                {
                  if(PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].HasCharges)
                  {
                    index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;
                    index_i2=-1;

                    f1_I=2.0*COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                         (-Ck[index1].im*Cksum.re+Ck[index1].re*Cksum.im);
                    f2_I=2.0*COULOMBIC_CONVERSION_FACTOR*exp_term*(4.0*M_PI/Volume[CurrentSystem])*
                         (-Ck[index1].re*Cksum.re-Ck[index1].im*Cksum.im);

                    if(index_i>=0)
                    {
                      HessianMatrix.element[index_i][index_i]+=f2_I*Rk.x*Rk.x;
                      HessianMatrix.element[index_i][index_i+1]+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i][index_i+2]+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+1][index_i+1]+=f2_I*Rk.y*Rk.y;
                      HessianMatrix.element[index_i+1][index_i+2]+=f2_I*Rk.y*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+2]+=f2_I*Rk.z*Rk.z;

                      HessianMatrix.element[index_i+1][index_i]+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i+2][index_i]+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+1]+=f2_I*Rk.y*Rk.z;

                      HessianAtomicPositionStrain(HessianMatrix,index_i,f1_I,Theta,Rk);
 
                      index2=0;
                      for(f2=0;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
                      {
                        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
                        {
                          if(PseudoAtoms[Framework[CurrentSystem].Atoms[f2][j].Type].HasCharges)
                          {
                            index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;
                            index_j2=-1;
                            if((index_i>=0)&&(index_j>=0))
                            {
                              if(index_i<=index_j)
                              {
                                f2_IJ=2.0*COULOMBIC_CONVERSION_FACTOR*exp_term*(4.0*M_PI/Volume[CurrentSystem])*
                                      (-Ck[index2].re*Ck[index1].re-Ck[index2].im*Ck[index1].im);
                                HessianMatrix.element[index_i][index_j]-=f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1]-=f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2]-=f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j]-=f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1]-=f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2]-=f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j]-=f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1]-=f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2]-=f2_IJ*Rk.z*Rk.z;
                              }
                            }
                            if((index_j>=0)||(index_j2>=0))
                              index2++;
                          }
                        }
                      }

                      for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
                      {
                        TypeMolB=Adsorbates[CurrentSystem][J].Type;
                        for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                        {
                          if(Components[TypeMolB].Groups[jg].Rigid)
                            comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                          for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                          {
                            j=Components[TypeMolB].Groups[jg].Atoms[ja];
                            if(PseudoAtoms[Adsorbates[CurrentSystem][J].Atoms[j].Type].HasCharges)
                            {
                              f2_IJ=8.0*COULOMBIC_CONVERSION_FACTOR*exp_term*(M_PI/Volume[CurrentSystem])*
                                    (-Ck[index2].re*Ck[index1].re-Ck[index2].im*Ck[index1].im);

                              if(Components[TypeMolB].Groups[jg].Rigid)
                                pos=Components[TypeMolB].Positions[j];

                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                                index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                                index2_rigid=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                                dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                                dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                                dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                              }
                              else
                              {
                                index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                                index_j2=-1;
                                index2_rigid=-1;

                                dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                              }

                              posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;
                              comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                              if(Components[TypeMolB].Groups[jg].Rigid)
                                HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                              if(index_i<=index_j)
                              {
                                if(index_j>=0)
                                {
                                  HessianMatrix.element[index_i][index_j]-=f2_IJ*Rk.x*Rk.x;
                                  HessianMatrix.element[index_i][index_j+1]-=f2_IJ*Rk.x*Rk.y;
                                  HessianMatrix.element[index_i][index_j+2]-=f2_IJ*Rk.x*Rk.z;
                                  HessianMatrix.element[index_i+1][index_j]-=f2_IJ*Rk.y*Rk.x;
                                  HessianMatrix.element[index_i+1][index_j+1]-=f2_IJ*Rk.y*Rk.y;
                                  HessianMatrix.element[index_i+1][index_j+2]-=f2_IJ*Rk.y*Rk.z;
                                  HessianMatrix.element[index_i+2][index_j]-=f2_IJ*Rk.z*Rk.x;
                                  HessianMatrix.element[index_i+2][index_j+1]-=f2_IJ*Rk.z*Rk.y;
                                  HessianMatrix.element[index_i+2][index_j+2]-=f2_IJ*Rk.z*Rk.z;
                                }

                                // com of I with orientation of J
                                if(Components[TypeMolB].Groups[jg].Rigid)
                                {
                                  if((index_i>=0)&&(index_j2>=0))
                                  {
                                    HessianMatrix.element[index_i][index_j2]-=f2_IJ*Rk.x*dot_product_j.x;
                                    HessianMatrix.element[index_i][index_j2+1]-=f2_IJ*Rk.x*dot_product_j.y;
                                    HessianMatrix.element[index_i][index_j2+2]-=f2_IJ*Rk.x*dot_product_j.z;
                                    HessianMatrix.element[index_i+1][index_j2]-=f2_IJ*Rk.y*dot_product_j.x;
                                    HessianMatrix.element[index_i+1][index_j2+1]-=f2_IJ*Rk.y*dot_product_j.y;
                                    HessianMatrix.element[index_i+1][index_j2+2]-=f2_IJ*Rk.y*dot_product_j.z;
                                    HessianMatrix.element[index_i+2][index_j2]-=f2_IJ*Rk.z*dot_product_j.x;
                                    HessianMatrix.element[index_i+2][index_j2+1]-=f2_IJ*Rk.z*dot_product_j.y;
                                    HessianMatrix.element[index_i+2][index_j2+2]-=f2_IJ*Rk.z*dot_product_j.z;
                                  }
                                }
                              }
                              if((index_j>=0)||(index_j2>=0))
                                index2++;
                            }
                          }
                        }
                      }

                      for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
                      {
                        TypeMolB=Cations[CurrentSystem][J].Type;
                        for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                        {
                          if(Components[TypeMolB].Groups[jg].Rigid)
                            comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                          for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                          {
                            j=Components[TypeMolB].Groups[jg].Atoms[ja];
                            if(PseudoAtoms[Cations[CurrentSystem][J].Atoms[j].Type].HasCharges)
                            {
                              f2_IJ=8.0*COULOMBIC_CONVERSION_FACTOR*exp_term*(M_PI/Volume[CurrentSystem])*
                                    (-Ck[index2].re*Ck[index1].re-Ck[index2].im*Ck[index1].im);

                              if(Components[TypeMolB].Groups[jg].Rigid)
                                pos=Components[TypeMolB].Positions[j];

                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                                index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;

                                index2_rigid=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                                dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                                dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                                dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                              }
                              else
                              {
                                index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                                index_j2=-1;
                                index2_rigid=-1;

                                dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                              }

                              posB=Cations[CurrentSystem][J].Atoms[j].Position;
                              comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                              // PLOP
                              if(Components[TypeMolB].Groups[jg].Rigid)
                                HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                              if(index_i<=index_j)
                              {
                                if((index_i>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i][index_j]-=f2_IJ*Rk.x*Rk.x;
                                  HessianMatrix.element[index_i][index_j+1]-=f2_IJ*Rk.x*Rk.y;
                                  HessianMatrix.element[index_i][index_j+2]-=f2_IJ*Rk.x*Rk.z;
                                  HessianMatrix.element[index_i+1][index_j]-=f2_IJ*Rk.y*Rk.x;
                                  HessianMatrix.element[index_i+1][index_j+1]-=f2_IJ*Rk.y*Rk.y;
                                  HessianMatrix.element[index_i+1][index_j+2]-=f2_IJ*Rk.y*Rk.z;
                                  HessianMatrix.element[index_i+2][index_j]-=f2_IJ*Rk.z*Rk.x;
                                  HessianMatrix.element[index_i+2][index_j+1]-=f2_IJ*Rk.z*Rk.y;
                                  HessianMatrix.element[index_i+2][index_j+2]-=f2_IJ*Rk.z*Rk.z;
                                }

                                // com of I with orientation of J
                                if(Components[TypeMolB].Groups[jg].Rigid)
                                {
                                  if((index_i>=0)&&(index_j2>=0))
                                  {
                                    HessianMatrix.element[index_i][index_j2]-=f2_IJ*Rk.x*dot_product_j.x;
                                    HessianMatrix.element[index_i][index_j2+1]-=f2_IJ*Rk.x*dot_product_j.y;
                                    HessianMatrix.element[index_i][index_j2+2]-=f2_IJ*Rk.x*dot_product_j.z;
                                    HessianMatrix.element[index_i+1][index_j2]-=f2_IJ*Rk.y*dot_product_j.x;
                                    HessianMatrix.element[index_i+1][index_j2+1]-=f2_IJ*Rk.y*dot_product_j.y;
                                    HessianMatrix.element[index_i+1][index_j2+2]-=f2_IJ*Rk.y*dot_product_j.z;
                                    HessianMatrix.element[index_i+2][index_j2]-=f2_IJ*Rk.z*dot_product_j.x;
                                    HessianMatrix.element[index_i+2][index_j2+1]-=f2_IJ*Rk.z*dot_product_j.y;
                                    HessianMatrix.element[index_i+2][index_j2+2]-=f2_IJ*Rk.z*dot_product_j.z;
                                  }
                                }
                              }
                              if((index_j>=0)||(index_j2>=0))
                                index2++;
                            }
                          }
                        }
                      }
                      if((index_i>=0)||(index_i2>=0))
                        index1++;
                    }
                  }//added
                }
              }
            }

            // AF, AA and AC
            for(I=0;I<NumberOfAdsorbateMolecules[CurrentSystem];I++)
            {
              TypeMolA=Adsorbates[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                  comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                {
                  i=Components[TypeMolA].Groups[ig].Atoms[ia];
                  if(PseudoAtoms[Adsorbates[CurrentSystem][I].Atoms[i].Type].HasCharges)
                  {
                    f1_I=2.0*COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                        (-Ck[index1].im*Cksum.re+Ck[index1].re*Cksum.im);
                    f2_I=(8.0*COULOMBIC_CONVERSION_FACTOR*exp_term*M_PI/Volume[CurrentSystem])*
                         (-Ck[index1].re*Cksum.re-Ck[index1].im*Cksum.im);

                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      index_i=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndex;
                      index_i2=Adsorbates[CurrentSystem][I].Groups[ig].HessianIndexOrientation;

                      index1_rigid=Adsorbates[CurrentSystem][I].Atoms[i].HessianAtomIndex;

                      dot_product_i.x=Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z;
                      dot_product_i.y=Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z;
                      dot_product_i.z=Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z;

                      dot_product_AX=Rk.x*DDVecAX[index1_rigid].x+Rk.y*DDVecAX[index1_rigid].y+Rk.z*DDVecAX[index1_rigid].z;
                      dot_product_BY=Rk.x*DDVecBY[index1_rigid].x+Rk.y*DDVecBY[index1_rigid].y+Rk.z*DDVecBY[index1_rigid].z;
                      dot_product_CZ=Rk.x*DDVecCZ[index1_rigid].x+Rk.y*DDVecCZ[index1_rigid].y+Rk.z*DDVecCZ[index1_rigid].z;
                      dot_product_AY=Rk.x*DDVecAY[index1_rigid].x+Rk.y*DDVecAY[index1_rigid].y+Rk.z*DDVecAY[index1_rigid].z;
                      dot_product_AZ=Rk.x*DDVecAZ[index1_rigid].x+Rk.y*DDVecAZ[index1_rigid].y+Rk.z*DDVecAZ[index1_rigid].z;
                      dot_product_BZ=Rk.x*DDVecBZ[index1_rigid].x+Rk.y*DDVecBZ[index1_rigid].y+Rk.z*DDVecBZ[index1_rigid].z;
                    }
                    else
                    {
                      index_i=Adsorbates[CurrentSystem][I].Atoms[i].HessianIndex;
                      index_i2=-1;
                      index1_rigid=-1;

                      dot_product_i.x=dot_product_i.y=dot_product_i.z=0.0;
                      dot_product_AX=dot_product_BY=dot_product_CZ=0.0;
                      dot_product_AY=dot_product_AZ=dot_product_BZ=0.0;
                    }

                    posA=Adsorbates[CurrentSystem][I].Atoms[i].Position;
                    comA=Adsorbates[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                    if(index_i>=0)
                    {
                      HessianMatrix.element[index_i][index_i]+=f2_I*Rk.x*Rk.x;
                      HessianMatrix.element[index_i][index_i+1]+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i][index_i+2]+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+1][index_i+1]+=f2_I*Rk.y*Rk.y;
                      HessianMatrix.element[index_i+1][index_i+2]+=f2_I*Rk.y*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+2]+=f2_I*Rk.z*Rk.z;

                      HessianMatrix.element[index_i+1][index_i]+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i+2][index_i]+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+1]+=f2_I*Rk.y*Rk.z;
                    }


                    // Crossterm: derivative of the energy with respect to strain and position
                    HessianAtomicPositionStrain(HessianMatrix,index_i,f1_I,Theta,Rk);


                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      pos=Components[TypeMolA].Positions[i];

                      HessianCenterOfMassStrainI(HessianMatrix,index_i,f2_I,posA,comA,Rk);

                      if(index_i2>=0)
                      {
                        HessianMatrix.element[index_i2][index_i2]+=f2_I*SQR(dot_product_i.x);
                        HessianMatrix.element[index_i2+1][index_i2+1]+=f2_I*SQR(dot_product_i.y);
                        HessianMatrix.element[index_i2+2][index_i2+2]+=f2_I*SQR(dot_product_i.z);
                        HessianMatrix.element[index_i2][index_i2+1]+=f2_I*dot_product_i.x*dot_product_i.y;
                        HessianMatrix.element[index_i2][index_i2+2]+=f2_I*dot_product_i.x*dot_product_i.z;
                        HessianMatrix.element[index_i2+1][index_i2+2]+=f2_I*dot_product_i.y*dot_product_i.z;
                      }

                      if((index_i>=0)&&(index_i2>=0))
                      {
                        HessianMatrix.element[index_i][index_i2]+=f2_I*Rk.x*dot_product_i.x;
                        HessianMatrix.element[index_i][index_i2+1]+=f2_I*Rk.x*dot_product_i.y;
                        HessianMatrix.element[index_i][index_i2+2]+=f2_I*Rk.x*dot_product_i.z;
                        HessianMatrix.element[index_i+1][index_i2]+=f2_I*Rk.y*dot_product_i.x;
                        HessianMatrix.element[index_i+1][index_i2+1]+=f2_I*Rk.y*dot_product_i.y;
                        HessianMatrix.element[index_i+1][index_i2+2]+=f2_I*Rk.y*dot_product_i.z;
                        HessianMatrix.element[index_i+2][index_i2]+=f2_I*Rk.z*dot_product_i.x;
                        HessianMatrix.element[index_i+2][index_i2+1]+=f2_I*Rk.z*dot_product_i.y;
                        HessianMatrix.element[index_i+2][index_i2+2]+=f2_I*Rk.z*dot_product_i.z;
                      }

                      if(index_i2>=0)
                      {
                        HessianMatrix.element[index_i2][index_i2]+=f1_I*dot_product_AX;
                        HessianMatrix.element[index_i2+1][index_i2+1]+=f1_I*dot_product_BY;
                        HessianMatrix.element[index_i2+2][index_i2+2]+=f1_I*dot_product_CZ;
                        HessianMatrix.element[index_i2][index_i2+1]+=f1_I*dot_product_AY;
                        HessianMatrix.element[index_i2][index_i2+2]+=f1_I*dot_product_AZ;
                        HessianMatrix.element[index_i2+1][index_i2+2]+=f1_I*dot_product_BZ;
                      }

                      HessianOrientationStrainI(HessianMatrix,index_i2,index1_rigid,f1_I,f2_I,posA,comA,Rk,Theta);

                      // correction term for first Born term
                      HessianAtomicCorrectionStrainStrainI(HessianMatrix,f1_I,f2_I,posA,comA,Rk,Theta);
                    }

                    index2=0;
                    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
                         for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
                         {
                           if(PseudoAtoms[Framework[CurrentSystem].Atoms[f1][j].Type].HasCharges)
                           {
                             index_j=Framework[CurrentSystem].Atoms[f1][j].HessianIndex;
                             index_j2=-1;

                             if((index_j>=0)||(index_j2>=0))
                               index2++;
                           }
                         }

                    for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Adsorbates[CurrentSystem][J].Type;
 
                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Adsorbates[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          if(PseudoAtoms[Adsorbates[CurrentSystem][J].Atoms[j].Type].HasCharges)
                          {
                            f2_IJ=(8.0*COULOMBIC_CONVERSION_FACTOR*exp_term*M_PI/Volume[CurrentSystem])*
                                  (-Ck[index2].re*Ck[index1].re-Ck[index2].im*Ck[index1].im);
                            if(Components[TypeMolB].Groups[jg].Rigid)
                              pos=Components[TypeMolB].Positions[j];


                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                              index2_rigid=Adsorbates[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=-1;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Adsorbates[CurrentSystem][J].Atoms[j].Position;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_i2,index1_rigid,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posA,comA,posB,comB,Rk);

                            if(index_i<=index_j)
                            {
                              if((index_i>=0)&&(index_j>=0))
                              {
                                HessianMatrix.element[index_i][index_j]-=f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1]-=f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2]-=f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j]-=f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1]-=f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2]-=f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j]-=f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1]-=f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2]-=f2_IJ*Rk.z*Rk.z;
                              }

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i][index_j2]-=f2_IJ*Rk.x*dot_product_j.x;
                                  HessianMatrix.element[index_i][index_j2+1]-=f2_IJ*Rk.x*dot_product_j.y;
                                  HessianMatrix.element[index_i][index_j2+2]-=f2_IJ*Rk.x*dot_product_j.z;
                                  HessianMatrix.element[index_i+1][index_j2]-=f2_IJ*Rk.y*dot_product_j.x;
                                  HessianMatrix.element[index_i+1][index_j2+1]-=f2_IJ*Rk.y*dot_product_j.y;
                                  HessianMatrix.element[index_i+1][index_j2+2]-=f2_IJ*Rk.y*dot_product_j.z;
                                  HessianMatrix.element[index_i+2][index_j2]-=f2_IJ*Rk.z*dot_product_j.x;
                                  HessianMatrix.element[index_i+2][index_j2+1]-=f2_IJ*Rk.z*dot_product_j.y;
                                  HessianMatrix.element[index_i+2][index_j2+2]-=f2_IJ*Rk.z*dot_product_j.z;
                                }
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j]-=f2_IJ*Rk.x*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j]-=f2_IJ*Rk.x*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j]-=f2_IJ*Rk.x*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+1]-=f2_IJ*Rk.y*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+1]-=f2_IJ*Rk.y*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+1]-=f2_IJ*Rk.y*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+2]-=f2_IJ*Rk.z*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+2]-=f2_IJ*Rk.z*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+2]-=f2_IJ*Rk.z*dot_product_i.z;
                                }
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j2]-=f2_IJ*dot_product_i.x*dot_product_j.x;
                                  HessianMatrix.element[index_i2][index_j2+1]-=f2_IJ*dot_product_i.x*dot_product_j.y;
                                  HessianMatrix.element[index_i2][index_j2+2]-=f2_IJ*dot_product_i.x*dot_product_j.z;
                                  HessianMatrix.element[index_i2+1][index_j2]-=f2_IJ*dot_product_i.y*dot_product_j.x;
                                  HessianMatrix.element[index_i2+1][index_j2+1]-=f2_IJ*dot_product_i.y*dot_product_j.y;
                                  HessianMatrix.element[index_i2+1][index_j2+2]-=f2_IJ*dot_product_i.y*dot_product_j.z;
                                  HessianMatrix.element[index_i2+2][index_j2]-=f2_IJ*dot_product_i.z*dot_product_j.x;
                                  HessianMatrix.element[index_i2+2][index_j2+1]-=f2_IJ*dot_product_i.z*dot_product_j.y;
                                  HessianMatrix.element[index_i2+2][index_j2+2]-=f2_IJ*dot_product_i.z*dot_product_j.z;
                                }
                              }
                            }
                            if((index_j>=0)||(index_j2>=0))
                              index2++;
                          }
                        }
                      }
                    }

                    for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Cations[CurrentSystem][J].Type;

                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          if(PseudoAtoms[Cations[CurrentSystem][J].Atoms[j].Type].HasCharges)
                          {
                            f2_IJ=(8.0*COULOMBIC_CONVERSION_FACTOR*exp_term*M_PI/Volume[CurrentSystem])*
                                  (-Ck[index2].re*Ck[index1].re-Ck[index2].im*Ck[index1].im);

                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                              index2_rigid=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;
                              pos=Components[TypeMolB].Positions[j];

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=-1;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Cations[CurrentSystem][J].Atoms[j].Position;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_i2,index1_rigid,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posA,comA,posB,comB,Rk);

                            // TESTING: CATION-ADSORBATE
                            if(Components[TypeMolA].Groups[ig].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_j,f2_IJ,posA,comA,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_j2,index2_rigid,f2_IJ,posA,comA,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posB,comB,posA,comA,Rk);

                            if(index_i<=index_j)
                            {
                              if((index_i>=0)&&(index_j>=0))
                              {
                                HessianMatrix.element[index_i][index_j]-=f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1]-=f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2]-=f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j]-=f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1]-=f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2]-=f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j]-=f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1]-=f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2]-=f2_IJ*Rk.z*Rk.z;
                              }

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i][index_j2]-=f2_IJ*Rk.x*dot_product_j.x;
                                  HessianMatrix.element[index_i][index_j2+1]-=f2_IJ*Rk.x*dot_product_j.y;
                                  HessianMatrix.element[index_i][index_j2+2]-=f2_IJ*Rk.x*dot_product_j.z;
                                  HessianMatrix.element[index_i+1][index_j2]-=f2_IJ*Rk.y*dot_product_j.x;
                                  HessianMatrix.element[index_i+1][index_j2+1]-=f2_IJ*Rk.y*dot_product_j.y;
                                  HessianMatrix.element[index_i+1][index_j2+2]-=f2_IJ*Rk.y*dot_product_j.z;
                                  HessianMatrix.element[index_i+2][index_j2]-=f2_IJ*Rk.z*dot_product_j.x;
                                  HessianMatrix.element[index_i+2][index_j2+1]-=f2_IJ*Rk.z*dot_product_j.y;
                                  HessianMatrix.element[index_i+2][index_j2+2]-=f2_IJ*Rk.z*dot_product_j.z;
                                }
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j]-=f2_IJ*Rk.x*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j]-=f2_IJ*Rk.x*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j]-=f2_IJ*Rk.x*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+1]-=f2_IJ*Rk.y*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+1]-=f2_IJ*Rk.y*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+1]-=f2_IJ*Rk.y*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+2]-=f2_IJ*Rk.z*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+2]-=f2_IJ*Rk.z*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+2]-=f2_IJ*Rk.z*dot_product_i.z;
                                }
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j2]-=f2_IJ*dot_product_i.x*dot_product_j.x;
                                  HessianMatrix.element[index_i2][index_j2+1]-=f2_IJ*dot_product_i.x*dot_product_j.y;
                                  HessianMatrix.element[index_i2][index_j2+2]-=f2_IJ*dot_product_i.x*dot_product_j.z;
                                  HessianMatrix.element[index_i2+1][index_j2]-=f2_IJ*dot_product_i.y*dot_product_j.x;
                                  HessianMatrix.element[index_i2+1][index_j2+1]-=f2_IJ*dot_product_i.y*dot_product_j.y;
                                  HessianMatrix.element[index_i2+1][index_j2+2]-=f2_IJ*dot_product_i.y*dot_product_j.z;
                                  HessianMatrix.element[index_i2+2][index_j2]-=f2_IJ*dot_product_i.z*dot_product_j.x;
                                  HessianMatrix.element[index_i2+2][index_j2+1]-=f2_IJ*dot_product_i.z*dot_product_j.y;
                                  HessianMatrix.element[index_i2+2][index_j2+2]-=f2_IJ*dot_product_i.z*dot_product_j.z;
                                }
                              }
                            }
                            if((index_j>=0)||(index_j2>=0))
                              index2++;
                          }
                        }
                      }
                    }

                    if((index_i>=0)||(index_i2>=0))
                      index1++;
                  }
                }
              }
            }

            // CF, CA and CC
            for(I=0;I<NumberOfCationMolecules[CurrentSystem];I++)
            {
              TypeMolA=Cations[CurrentSystem][I].Type;
              for(ig=0;ig<Components[TypeMolA].NumberOfGroups;ig++)
              {
                if(Components[TypeMolA].Groups[ig].Rigid)
                  comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                for(ia=0;ia<Components[TypeMolA].Groups[ig].NumberOfGroupAtoms;ia++)
                {
                  i=Components[TypeMolA].Groups[ig].Atoms[ia];
                  if(PseudoAtoms[Cations[CurrentSystem][I].Atoms[i].Type].HasCharges)
                  {
                    f1_I=2.0*COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                        (-Ck[index1].im*Cksum.re+Ck[index1].re*Cksum.im);
                    f2_I=(8.0*COULOMBIC_CONVERSION_FACTOR*exp_term*M_PI/Volume[CurrentSystem])*
                         (-Ck[index1].re*Cksum.re-Ck[index1].im*Cksum.im);
                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      index_i=Cations[CurrentSystem][I].Groups[ig].HessianIndex;
                      index_i2=Cations[CurrentSystem][I].Groups[ig].HessianIndexOrientation;
                      index1_rigid=Cations[CurrentSystem][I].Atoms[i].HessianAtomIndex;

                      dot_product_i.x=Rk.x*DVecX[index1_rigid].x+Rk.y*DVecX[index1_rigid].y+Rk.z*DVecX[index1_rigid].z;
                      dot_product_i.y=Rk.x*DVecY[index1_rigid].x+Rk.y*DVecY[index1_rigid].y+Rk.z*DVecY[index1_rigid].z;
                      dot_product_i.z=Rk.x*DVecZ[index1_rigid].x+Rk.y*DVecZ[index1_rigid].y+Rk.z*DVecZ[index1_rigid].z;

                      dot_product_AX=Rk.x*DDVecAX[index1_rigid].x+Rk.y*DDVecAX[index1_rigid].y+Rk.z*DDVecAX[index1_rigid].z;
                      dot_product_BY=Rk.x*DDVecBY[index1_rigid].x+Rk.y*DDVecBY[index1_rigid].y+Rk.z*DDVecBY[index1_rigid].z;
                      dot_product_CZ=Rk.x*DDVecCZ[index1_rigid].x+Rk.y*DDVecCZ[index1_rigid].y+Rk.z*DDVecCZ[index1_rigid].z;
                      dot_product_AY=Rk.x*DDVecAY[index1_rigid].x+Rk.y*DDVecAY[index1_rigid].y+Rk.z*DDVecAY[index1_rigid].z;
                      dot_product_AZ=Rk.x*DDVecAZ[index1_rigid].x+Rk.y*DDVecAZ[index1_rigid].y+Rk.z*DDVecAZ[index1_rigid].z;
                      dot_product_BZ=Rk.x*DDVecBZ[index1_rigid].x+Rk.y*DDVecBZ[index1_rigid].y+Rk.z*DDVecBZ[index1_rigid].z;
                    }
                    else
                    {
                      index_i=Cations[CurrentSystem][I].Atoms[i].HessianIndex;
                      index_i2=-1;
                      index1_rigid=-1;

                      dot_product_i.x=dot_product_i.y=dot_product_i.z=0.0;
                      dot_product_AX=dot_product_BY=dot_product_CZ=0.0;
                      dot_product_AY=dot_product_AZ=dot_product_BZ=0.0;
                    }

                    posA=Cations[CurrentSystem][I].Atoms[i].Position;
                    comA=Cations[CurrentSystem][I].Groups[ig].CenterOfMassPosition;

                    if(index_i>=0)
                    {
                      HessianMatrix.element[index_i][index_i]+=f2_I*Rk.x*Rk.x;
                      HessianMatrix.element[index_i][index_i+1]+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i][index_i+2]+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+1][index_i+1]+=f2_I*Rk.y*Rk.y;
                      HessianMatrix.element[index_i+1][index_i+2]+=f2_I*Rk.y*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+2]+=f2_I*Rk.z*Rk.z;

                      HessianMatrix.element[index_i+1][index_i]+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i+2][index_i]+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+1]+=f2_I*Rk.y*Rk.z;
                    }


                    // Crossterm: derivative of the energy with respect to strain and position
                    HessianAtomicPositionStrain(HessianMatrix,index_i,f1_I,Theta,Rk);

                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      pos=Components[TypeMolA].Positions[i];

                      HessianCenterOfMassStrainI(HessianMatrix,index_i,f2_I,posA,comA,Rk);

                      if(index_i2>=0)
                      {
                        HessianMatrix.element[index_i2][index_i2]+=f2_I*dot_product_i.x*dot_product_i.x;
                        HessianMatrix.element[index_i2+1][index_i2+1]+=f2_I*dot_product_i.y*dot_product_i.y;
                        HessianMatrix.element[index_i2+2][index_i2+2]+=f2_I*dot_product_i.z*dot_product_i.z;
                        HessianMatrix.element[index_i2][index_i2+1]+=f2_I*dot_product_i.x*dot_product_i.y;
                        HessianMatrix.element[index_i2][index_i2+2]+=f2_I*dot_product_i.x*dot_product_i.z;
                        HessianMatrix.element[index_i2+1][index_i2+2]+=f2_I*dot_product_i.y*dot_product_i.z;
                      }

                      if((index_i>=0)&&(index_i2>=0))
                      { 
                        HessianMatrix.element[index_i][index_i2]+=f2_I*Rk.x*dot_product_i.x;
                        HessianMatrix.element[index_i][index_i2+1]+=f2_I*Rk.x*dot_product_i.y;
                        HessianMatrix.element[index_i][index_i2+2]+=f2_I*Rk.x*dot_product_i.z;
                        HessianMatrix.element[index_i+1][index_i2]+=f2_I*Rk.y*dot_product_i.x;
                        HessianMatrix.element[index_i+1][index_i2+1]+=f2_I*Rk.y*dot_product_i.y;
                        HessianMatrix.element[index_i+1][index_i2+2]+=f2_I*Rk.y*dot_product_i.z;
                        HessianMatrix.element[index_i+2][index_i2]+=f2_I*Rk.z*dot_product_i.x;
                        HessianMatrix.element[index_i+2][index_i2+1]+=f2_I*Rk.z*dot_product_i.y;
                        HessianMatrix.element[index_i+2][index_i2+2]+=f2_I*Rk.z*dot_product_i.z;
                      }

                      if(index_i2>=0)
                      {
                        pos=Components[TypeMolA].Positions[i];
                        HessianMatrix.element[index_i2][index_i2]+=f1_I*dot_product_AX;
                        HessianMatrix.element[index_i2+1][index_i2+1]+=f1_I*dot_product_BY;
                        HessianMatrix.element[index_i2+2][index_i2+2]+=f1_I*dot_product_CZ;
                        HessianMatrix.element[index_i2][index_i2+1]+=f1_I*dot_product_AY;
                        HessianMatrix.element[index_i2][index_i2+2]+=f1_I*dot_product_AZ;
                        HessianMatrix.element[index_i2+1][index_i2+2]+=f1_I*dot_product_BZ;
                      }

                      // derivative of stress with respect to the orientation 
                      HessianOrientationStrainI(HessianMatrix,index_i2,index1_rigid,f1_I,f2_I,posA,comA,Rk,Theta);

                      // correction term for first Born term
                      HessianAtomicCorrectionStrainStrainI(HessianMatrix,f1_I,f2_I,posA,comA,Rk,Theta);
                    }

                    index2=0;
                    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
                      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
                         for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
                         {
                           index_j=Framework[CurrentSystem].Atoms[f1][j].HessianIndex;
                           index_j2=-1;

                           if((index_j>=0)||(index_j2>=0))
                             index2++;
                         }

                    for(J=0;J<NumberOfAdsorbateMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Adsorbates[CurrentSystem][J].Type;
 
                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          if(PseudoAtoms[Adsorbates[CurrentSystem][J].Atoms[j].Type].HasCharges)
                          {
                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Adsorbates[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                            }
                            else
                            {
                              index_j=Adsorbates[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=-1;
                            }

                            // SKIPPING

                            if((index_j>=0)||(index_j2>=0))
                              index2++;
                          }
                        }
                      }
                    }

                    for(J=0;J<NumberOfCationMolecules[CurrentSystem];J++)
                    {
                      TypeMolB=Cations[CurrentSystem][J].Type;

                      for(jg=0;jg<Components[TypeMolB].NumberOfGroups;jg++)
                      {
                        if(Components[TypeMolB].Groups[jg].Rigid)
                          comB=Cations[CurrentSystem][J].Groups[jg].CenterOfMassPosition;

                        for(ja=0;ja<Components[TypeMolB].Groups[jg].NumberOfGroupAtoms;ja++)
                        {
                          j=Components[TypeMolB].Groups[jg].Atoms[ja];
                          if(PseudoAtoms[Cations[CurrentSystem][J].Atoms[j].Type].HasCharges)
                          {
                            f2_IJ=8.0*COULOMBIC_CONVERSION_FACTOR*exp_term*(M_PI/Volume[CurrentSystem])*
                                  (-Ck[index2].re*Ck[index1].re-Ck[index2].im*Ck[index1].im);
                            if(Components[TypeMolB].Groups[jg].Rigid)
                              pos=Components[TypeMolB].Positions[j];


                            if(Components[TypeMolB].Groups[jg].Rigid)
                            {
                              index_j=Cations[CurrentSystem][J].Groups[jg].HessianIndex;
                              index_j2=Cations[CurrentSystem][J].Groups[jg].HessianIndexOrientation;
                              index2_rigid=Cations[CurrentSystem][J].Atoms[j].HessianAtomIndex;

                              dot_product_j.x=Rk.x*DVecX[index2_rigid].x+Rk.y*DVecX[index2_rigid].y+Rk.z*DVecX[index2_rigid].z;
                              dot_product_j.y=Rk.x*DVecY[index2_rigid].x+Rk.y*DVecY[index2_rigid].y+Rk.z*DVecY[index2_rigid].z;
                              dot_product_j.z=Rk.x*DVecZ[index2_rigid].x+Rk.y*DVecZ[index2_rigid].y+Rk.z*DVecZ[index2_rigid].z;
                            }
                            else
                            {
                              index_j=Cations[CurrentSystem][J].Atoms[j].HessianIndex;
                              index_j2=-1;
                              index2_rigid=-1;

                              dot_product_j.x=dot_product_j.y=dot_product_j.z=0.0;
                            }

                            posB=Cations[CurrentSystem][J].Atoms[j].Position;

                            if(Components[TypeMolB].Groups[jg].Rigid)
                              HessianCenterOfMassStrainJ(HessianMatrix,index_i,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianOrientationStrainJ(HessianMatrix,index_i2,index1_rigid,f2_IJ,posB,comB,Rk);

                            if(Components[TypeMolA].Groups[ig].Rigid&&Components[TypeMolB].Groups[jg].Rigid)
                              HessianAtomicCorrectionStrainStrainJ(HessianMatrix,f2_IJ,posA,comA,posB,comB,Rk);

                            if(index_i<=index_j)
                            {
                              if((index_i>=0)&(index_j>=0))
                              { 
                                HessianMatrix.element[index_i][index_j]-=f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1]-=f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2]-=f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j]-=f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1]-=f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2]-=f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j]-=f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1]-=f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2]-=f2_IJ*Rk.z*Rk.z;
                              }

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i][index_j2]-=f2_IJ*Rk.x*dot_product_j.x;
                                  HessianMatrix.element[index_i][index_j2+1]-=f2_IJ*Rk.x*dot_product_j.y;
                                  HessianMatrix.element[index_i][index_j2+2]-=f2_IJ*Rk.x*dot_product_j.z;
                                  HessianMatrix.element[index_i+1][index_j2]-=f2_IJ*Rk.y*dot_product_j.x;
                                  HessianMatrix.element[index_i+1][index_j2+1]-=f2_IJ*Rk.y*dot_product_j.y;
                                  HessianMatrix.element[index_i+1][index_j2+2]-=f2_IJ*Rk.y*dot_product_j.z;
                                  HessianMatrix.element[index_i+2][index_j2]-=f2_IJ*Rk.z*dot_product_j.x;
                                  HessianMatrix.element[index_i+2][index_j2+1]-=f2_IJ*Rk.z*dot_product_j.y;
                                  HessianMatrix.element[index_i+2][index_j2+2]-=f2_IJ*Rk.z*dot_product_j.z;
                                }
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j]-=f2_IJ*Rk.x*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j]-=f2_IJ*Rk.x*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j]-=f2_IJ*Rk.x*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+1]-=f2_IJ*Rk.y*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+1]-=f2_IJ*Rk.y*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+1]-=f2_IJ*Rk.y*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+2]-=f2_IJ*Rk.z*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+2]-=f2_IJ*Rk.z*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+2]-=f2_IJ*Rk.z*dot_product_i.z;
                                }
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j2]-=f2_IJ*dot_product_i.x*dot_product_j.x;
                                  HessianMatrix.element[index_i2][index_j2+1]-=f2_IJ*dot_product_i.x*dot_product_j.y;
                                  HessianMatrix.element[index_i2][index_j2+2]-=f2_IJ*dot_product_i.x*dot_product_j.z;
                                  HessianMatrix.element[index_i2+1][index_j2]-=f2_IJ*dot_product_i.y*dot_product_j.x;
                                  HessianMatrix.element[index_i2+1][index_j2+1]-=f2_IJ*dot_product_i.y*dot_product_j.y;
                                  HessianMatrix.element[index_i2+1][index_j2+2]-=f2_IJ*dot_product_i.y*dot_product_j.z;
                                  HessianMatrix.element[index_i2+2][index_j2]-=f2_IJ*dot_product_i.z*dot_product_j.x;
                                  HessianMatrix.element[index_i2+2][index_j2+1]-=f2_IJ*dot_product_i.z*dot_product_j.y;
                                  HessianMatrix.element[index_i2+2][index_j2+2]-=f2_IJ*dot_product_i.z*dot_product_j.z;
                                }
                              }
                            }
                            if((index_j>=0)||(index_j2>=0))
                              index2++;
                          }
                        }
                      }
                    }
                    if((index_i>=0)||(index_i2>=0))
                      index1++;
                  }
                }
              }
            }

          }


          NKVec++;
        }
      }
      Nmin=-kvec[CurrentSystem].z;
    }
    Mmin=-kvec[CurrentSystem].y;
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {

    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];i++)
      {
        A=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i].B;

        index_i=Framework[CurrentSystem].Atoms[f1][A].HessianIndex;
        index_j=Framework[CurrentSystem].Atoms[f1][B].HessianIndex;

        typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
        typeB=Framework[CurrentSystem].Atoms[f1][B].Type;

        ChargeA=PseudoAtoms[typeA].Charge;
        ChargeB=PseudoAtoms[typeB].Charge;

        posA=Framework[CurrentSystem].Atoms[f1][A].Position;
        posB=Framework[CurrentSystem].Atoms[f1][B].Position;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        // Coulomb subtraction has a problem at small distances, leading to large negative values in the Hessian
        // i.e. minimization of Coesite fails (a high density mineral)
        // if the distance is smaller than 1e-4 Angsgtrom, use the limit values
        if(r>1e-4)
        {
          // add contribution to the energy
          (*Energy)-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                                ChargeA*ChargeB/r;

          DF=-(COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
             (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
             (r*rr));

          DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
              (((-3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(r*rr))+
              (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
             (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
        }
        else
        {
          (*Energy)-=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*2.0*Alpha[CurrentSystem]/sqrt(M_PI);
          DF=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*4.0*CUBE(Alpha[CurrentSystem])/(3.0*sqrt(M_PI));
          DDF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*8.0*CUBE(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])/(5.0*sqrt(M_PI));
        }

        if((index_i<0)&&(index_j<0)) continue;

        S.ax=-dr.x*dr.x;
        S.bx=-dr.y*dr.x;
        S.cx=-dr.z*dr.x;

        S.ay=-dr.x*dr.y;
        S.by=-dr.y*dr.y;
        S.cy=-dr.z*dr.y;

        S.az=-dr.x*dr.z;
        S.bz=-dr.y*dr.z;
        S.cz=-dr.z*dr.z;

        // add contribution to the first derivatives
        if(ComputeGradient)
        {
          if(index_i>=0)
          {
            Gradient[index_i]+=DF*dr.x;
            Gradient[index_i+1]+=DF*dr.y;
            Gradient[index_i+2]+=DF*dr.z;
          }

          if(index_j>=0)
          {
            Gradient[index_j]-=DF*dr.x;
            Gradient[index_j+1]-=DF*dr.y;
            Gradient[index_j+2]-=DF*dr.z;
          }

          GradientStrain(DF,Gradient,S);
        }

        StrainFirstDerivative->ax-=DF*S.ax;
        StrainFirstDerivative->bx-=DF*S.bx;
        StrainFirstDerivative->cx-=DF*S.cx;

        StrainFirstDerivative->ay-=DF*S.ay;
        StrainFirstDerivative->by-=DF*S.by;
        StrainFirstDerivative->cy-=DF*S.cy;

        StrainFirstDerivative->az-=DF*S.az;
        StrainFirstDerivative->bz-=DF*S.bz;
        StrainFirstDerivative->cz-=DF*S.cz;

        if(ComputeHessian)
        {
          Hessian.ax=DDF*dr.x*dr.x+DF; Hessian.bx=DDF*dr.y*dr.x;    Hessian.cx=DDF*dr.z*dr.x;
          Hessian.ay=DDF*dr.x*dr.y;    Hessian.by=DDF*dr.y*dr.y+DF; Hessian.cy=DDF*dr.z*dr.y;
          Hessian.az=DDF*dr.x*dr.z;    Hessian.bz=DDF*dr.y*dr.z;    Hessian.cz=DDF*dr.z*dr.z+DF;

          if(index_i>=0)
          { 
            HessianMatrix.element[index_i][index_i]+=Hessian.ax;
            HessianMatrix.element[index_i][index_i+1]+=Hessian.ay;
            HessianMatrix.element[index_i][index_i+2]+=Hessian.az;
            HessianMatrix.element[index_i+1][index_i+1]+=Hessian.by;
            HessianMatrix.element[index_i+1][index_i+2]+=Hessian.bz;
            HessianMatrix.element[index_i+2][index_i+2]+=Hessian.cz;
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j][index_j]+=Hessian.ax;
            HessianMatrix.element[index_j][index_j+1]+=Hessian.ay;
            HessianMatrix.element[index_j][index_j+2]+=Hessian.az;
            HessianMatrix.element[index_j+1][index_j+1]+=Hessian.by;
            HessianMatrix.element[index_j+1][index_j+2]+=Hessian.bz;
            HessianMatrix.element[index_j+2][index_j+2]+=Hessian.cz;
          }

          if((index_i>=0)&&(index_j>=0))
          {
            if(index_i<index_j)
            {
              HessianMatrix.element[index_i][index_j]-=Hessian.ax;
              HessianMatrix.element[index_i][index_j+1]-=Hessian.ay;
              HessianMatrix.element[index_i][index_j+2]-=Hessian.az;
              HessianMatrix.element[index_i+1][index_j]-=Hessian.ay;
              HessianMatrix.element[index_i+1][index_j+1]-=Hessian.by;
              HessianMatrix.element[index_i+1][index_j+2]-=Hessian.bz;
              HessianMatrix.element[index_i+2][index_j]-=Hessian.az;
              HessianMatrix.element[index_i+2][index_j+1]-=Hessian.bz;
              HessianMatrix.element[index_i+2][index_j+2]-=Hessian.cz;
            }
            else
            {
              HessianMatrix.element[index_j][index_i]-=Hessian.ax;
              HessianMatrix.element[index_j][index_i+1]-=Hessian.ay;
              HessianMatrix.element[index_j][index_i+2]-=Hessian.az;
              HessianMatrix.element[index_j+1][index_i]-=Hessian.ay;
              HessianMatrix.element[index_j+1][index_i+1]-=Hessian.by;
              HessianMatrix.element[index_j+1][index_i+2]-=Hessian.bz;
              HessianMatrix.element[index_j+2][index_i]-=Hessian.az;
              HessianMatrix.element[index_j+2][index_i+1]-=Hessian.bz;
              HessianMatrix.element[index_j+2][index_i+2]-=Hessian.cz;
            }
          }

          HessianAtomicPositionStrainExcluded(HessianMatrix,index_i,index_j,DF,DDF,dr);
          HessianAtomicStrainStrainExcluded(HessianMatrix,index_i,index_j,DF,DDF,dr);
        }
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        A=Components[Type].ExcludedIntraChargeCharge[i].A;
        B=Components[Type].ExcludedIntraChargeCharge[i].B;

        TypeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
        if(PseudoAtoms[TypeA].HasCharges)
        {
          TypeB=Adsorbates[CurrentSystem][m].Atoms[B].Type;
          if(PseudoAtoms[TypeB].HasCharges)
          {
            ChargeA=PseudoAtoms[TypeA].Charge;
            ChargeB=PseudoAtoms[TypeB].Charge;

            posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
            posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

            grpA=Components[Type].group[A];
            grpB=Components[Type].group[B];

            if(Components[Type].Groups[grpA].Rigid)
              index_i=Adsorbates[CurrentSystem][m].Groups[grpA].HessianIndex;
            else
              index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;

            if(Components[Type].Groups[grpB].Rigid)
              index_j=Adsorbates[CurrentSystem][m].Groups[grpB].HessianIndex;
            else
              index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            // add contribution to the energy
            (*Energy)-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                               ChargeA*ChargeB/r;

            DF=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
               (ErrorFunction(Alpha[CurrentSystem]*r)-2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
               (r*rr);

            DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                 (((-3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(r*rr))+
                 (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                 (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));

            if((index_i<0)&&(index_j<0)) continue;

            // skip the remainder of the loop if both atoms belong to the same rigid group
            if((grpA==grpB)&&(Components[Type].Groups[grpA].Rigid)) continue;

            S.ax=-dr.x*dr.x;
            S.bx=-dr.y*dr.x;
            S.cx=-dr.z*dr.x;

            S.ay=-dr.x*dr.y;
            S.by=-dr.y*dr.y;
            S.cy=-dr.z*dr.y;

            S.az=-dr.x*dr.z;
            S.bz=-dr.y*dr.z;
            S.cz=-dr.z*dr.z;

            // add contribution to the first derivatives
            if(ComputeGradient)
            {
              if(index_i>=0)
              {
                Gradient[index_i]+=DF*dr.x;
                Gradient[index_i+1]+=DF*dr.y;
                Gradient[index_i+2]+=DF*dr.z;
              }

              if(index_j>=0)
              {
                Gradient[index_j]-=DF*dr.x;
                Gradient[index_j+1]-=DF*dr.y;
                Gradient[index_j+2]-=DF*dr.z;
              }

              GradientStrain(DF,Gradient,S);
            }

            StrainFirstDerivative->ax-=DF*S.ax;
            StrainFirstDerivative->bx-=DF*S.bx;
            StrainFirstDerivative->cx-=DF*S.cx;

            StrainFirstDerivative->ay-=DF*S.ay;
            StrainFirstDerivative->by-=DF*S.by;
            StrainFirstDerivative->cy-=DF*S.cy;

            StrainFirstDerivative->az-=DF*S.az;
            StrainFirstDerivative->bz-=DF*S.bz;
            StrainFirstDerivative->cz-=DF*S.cz;

            if(ComputeHessian)
            {
              Hessian.ax=DDF*dr.x*dr.x+DF; Hessian.bx=DDF*dr.y*dr.x;    Hessian.cx=DDF*dr.z*dr.x;
              Hessian.ay=DDF*dr.x*dr.y;    Hessian.by=DDF*dr.y*dr.y+DF; Hessian.cy=DDF*dr.z*dr.y;
              Hessian.az=DDF*dr.x*dr.z;    Hessian.bz=DDF*dr.y*dr.z;    Hessian.cz=DDF*dr.z*dr.z+DF;

              if(index_i>=0)
              { 
                HessianMatrix.element[index_i][index_i]+=Hessian.ax;
                HessianMatrix.element[index_i][index_i+1]+=Hessian.ay;
                HessianMatrix.element[index_i][index_i+2]+=Hessian.az;
                HessianMatrix.element[index_i+1][index_i+1]+=Hessian.by;
                HessianMatrix.element[index_i+1][index_i+2]+=Hessian.bz;
                HessianMatrix.element[index_i+2][index_i+2]+=Hessian.cz;
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][index_j]+=Hessian.ax;
                HessianMatrix.element[index_j][index_j+1]+=Hessian.ay;
                HessianMatrix.element[index_j][index_j+2]+=Hessian.az;
                HessianMatrix.element[index_j+1][index_j+1]+=Hessian.by;
                HessianMatrix.element[index_j+1][index_j+2]+=Hessian.bz;
                HessianMatrix.element[index_j+2][index_j+2]+=Hessian.cz;
              }

              if((index_i>=0)&&(index_j>=0))
              {
                if(index_i<index_j)
                {
                  HessianMatrix.element[index_i][index_j]-=Hessian.ax;
                  HessianMatrix.element[index_i][index_j+1]-=Hessian.ay;
                  HessianMatrix.element[index_i][index_j+2]-=Hessian.az;
                  HessianMatrix.element[index_i+1][index_j]-=Hessian.ay;
                  HessianMatrix.element[index_i+1][index_j+1]-=Hessian.by;
                  HessianMatrix.element[index_i+1][index_j+2]-=Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j]-=Hessian.az;
                  HessianMatrix.element[index_i+2][index_j+1]-=Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j+2]-=Hessian.cz;
                }
                else
                {
                  HessianMatrix.element[index_j][index_i]-=Hessian.ax;
                  HessianMatrix.element[index_j][index_i+1]-=Hessian.ay;
                  HessianMatrix.element[index_j][index_i+2]-=Hessian.az;
                  HessianMatrix.element[index_j+1][index_i]-=Hessian.ay;
                  HessianMatrix.element[index_j+1][index_i+1]-=Hessian.by;
                  HessianMatrix.element[index_j+1][index_i+2]-=Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i]-=Hessian.az;
                  HessianMatrix.element[index_j+2][index_i+1]-=Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i+2]-=Hessian.cz;
                }
              }

              HessianAtomicPositionStrainExcluded(HessianMatrix,index_i,index_j,DF,DDF,dr);
              HessianAtomicStrainStrainExcluded(HessianMatrix,index_i,index_j,DF,DDF,dr);
            }
          }
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        A=Components[Type].ExcludedIntraChargeCharge[i].A;
        B=Components[Type].ExcludedIntraChargeCharge[i].B;

        TypeA=Cations[CurrentSystem][m].Atoms[A].Type;
        if(PseudoAtoms[TypeA].HasCharges)
        {
          TypeB=Cations[CurrentSystem][m].Atoms[B].Type;
          if(PseudoAtoms[TypeB].HasCharges)
          {
            ChargeA=PseudoAtoms[TypeA].Charge;
            ChargeB=PseudoAtoms[TypeB].Charge;

            posA=Cations[CurrentSystem][m].Atoms[A].Position;
            posB=Cations[CurrentSystem][m].Atoms[B].Position;

            grpA=Components[Type].group[A];
            grpB=Components[Type].group[B];

            if(Components[Type].Groups[grpA].Rigid)
              index_i=Cations[CurrentSystem][m].Groups[grpA].HessianIndex;
            else
              index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;

            if(Components[Type].Groups[grpB].Rigid)
              index_j=Cations[CurrentSystem][m].Groups[grpB].HessianIndex;
            else
              index_j=Cations[CurrentSystem][m].Atoms[B].HessianIndex;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            // add contribution to the energy
            (*Energy)-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                               ChargeA*ChargeB/r;

            DF=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
               (ErrorFunction(Alpha[CurrentSystem]*r)-2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
               (r*rr);

            DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                 (((-3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(r*rr))+
                 (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                 (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));

            if((index_i<0)&&(index_j<0)) continue;

            // skip the remainder of the loop if both atoms belong to the same rigid group
            if((grpA==grpB)&&(Components[Type].Groups[grpA].Rigid)) continue;

            S.ax=-dr.x*dr.x;
            S.bx=-dr.y*dr.x;
            S.cx=-dr.z*dr.x;

            S.ay=-dr.x*dr.y;
            S.by=-dr.y*dr.y;
            S.cy=-dr.z*dr.y;

            S.az=-dr.x*dr.z;
            S.bz=-dr.y*dr.z;
            S.cz=-dr.z*dr.z;

            // add contribution to the first derivatives
            if(ComputeGradient)
            {
              if(index_i>=0)
              {
                Gradient[index_i]+=DF*dr.x;
                Gradient[index_i+1]+=DF*dr.y;
                Gradient[index_i+2]+=DF*dr.z;
              }

              if(index_j>=0)
              {
                Gradient[index_j]-=DF*dr.x;
                Gradient[index_j+1]-=DF*dr.y;
                Gradient[index_j+2]-=DF*dr.z;
              }

              GradientStrain(DF,Gradient,S);
            }

            StrainFirstDerivative->ax-=DF*S.ax;
            StrainFirstDerivative->bx-=DF*S.bx;
            StrainFirstDerivative->cx-=DF*S.cx;

            StrainFirstDerivative->ay-=DF*S.ay;
            StrainFirstDerivative->by-=DF*S.by;
            StrainFirstDerivative->cy-=DF*S.cy;

            StrainFirstDerivative->az-=DF*S.az;
            StrainFirstDerivative->bz-=DF*S.bz;
            StrainFirstDerivative->cz-=DF*S.cz;

            if(ComputeHessian)
            {
              Hessian.ax=DDF*dr.x*dr.x+DF; Hessian.bx=DDF*dr.y*dr.x;    Hessian.cx=DDF*dr.z*dr.x;
              Hessian.ay=DDF*dr.x*dr.y;    Hessian.by=DDF*dr.y*dr.y+DF; Hessian.cy=DDF*dr.z*dr.y;
              Hessian.az=DDF*dr.x*dr.z;    Hessian.bz=DDF*dr.y*dr.z;    Hessian.cz=DDF*dr.z*dr.z+DF;

              if(index_i>=0)
              { 
                HessianMatrix.element[index_i][index_i]+=Hessian.ax;
                HessianMatrix.element[index_i][index_i+1]+=Hessian.ay;
                HessianMatrix.element[index_i][index_i+2]+=Hessian.az;
                HessianMatrix.element[index_i+1][index_i+1]+=Hessian.by;
                HessianMatrix.element[index_i+1][index_i+2]+=Hessian.bz;
                HessianMatrix.element[index_i+2][index_i+2]+=Hessian.cz;
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][index_j]+=Hessian.ax;
                HessianMatrix.element[index_j][index_j+1]+=Hessian.ay;
                HessianMatrix.element[index_j][index_j+2]+=Hessian.az;
                HessianMatrix.element[index_j+1][index_j+1]+=Hessian.by;
                HessianMatrix.element[index_j+1][index_j+2]+=Hessian.bz;
                HessianMatrix.element[index_j+2][index_j+2]+=Hessian.cz;
              }

              if((index_i>=0)&&(index_j>=0))
              {
                if(index_i<index_j)
                {
                  HessianMatrix.element[index_i][index_j]-=Hessian.ax;
                  HessianMatrix.element[index_i][index_j+1]-=Hessian.ay;
                  HessianMatrix.element[index_i][index_j+2]-=Hessian.az;
                  HessianMatrix.element[index_i+1][index_j]-=Hessian.ay;
                  HessianMatrix.element[index_i+1][index_j+1]-=Hessian.by;
                  HessianMatrix.element[index_i+1][index_j+2]-=Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j]-=Hessian.az;
                  HessianMatrix.element[index_i+2][index_j+1]-=Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j+2]-=Hessian.cz;
                }
                else
                {
                  HessianMatrix.element[index_j][index_i]-=Hessian.ax;
                  HessianMatrix.element[index_j][index_i+1]-=Hessian.ay;
                  HessianMatrix.element[index_j][index_i+2]-=Hessian.az;
                  HessianMatrix.element[index_j+1][index_i]-=Hessian.ay;
                  HessianMatrix.element[index_j+1][index_i+1]-=Hessian.by;
                  HessianMatrix.element[index_j+1][index_i+2]-=Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i]-=Hessian.az;
                  HessianMatrix.element[index_j+2][index_i+1]-=Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i+2]-=Hessian.cz;
                }
              }

              HessianAtomicPositionStrainExcluded(HessianMatrix,index_i,index_j,DF,DDF,dr);
              HessianAtomicStrainStrainExcluded(HessianMatrix,index_i,index_j,DF,DDF,dr);
            }
          }
        }
      }
    }
  }

  Uself_sum=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        Uself_sum+=SQR(PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Charge);
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      Uself_sum+=SQR(PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].Charge);

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
      Uself_sum+=SQR(PseudoAtoms[Cations[CurrentSystem][i].Atoms[j].Type].Charge);

  Uself_sum*=Alpha[CurrentSystem]/sqrt(M_PI);

  *Energy+=USum-COULOMBIC_CONVERSION_FACTOR*Uself_sum;

  return 0;
}
