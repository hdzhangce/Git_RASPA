/*****************************************************************************************************
    ewald_exclusion.c -  description
    --------------------------------
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
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "potentials.h"
#include "spectra.h"

REAL UHostHostChargeChargeExcludedFourier;
REAL UHostHostChargeBondDipoleExcludedFourier;
REAL UHostHostBondDipoleBondDipoleExcludedFourier;

REAL UAdsorbateAdsorbateChargeChargeExcludedFourier;
REAL UAdsorbateAdsorbateChargeBondDipoleExcludedFourier;
REAL UAdsorbateAdsorbateBondDipoleBondDipoleExcludedFourier;

REAL UCationCationChargeChargeExcludedFourier;
REAL UCationCationChargeBondDipoleExcludedFourier;
REAL UCationCationBondDipoleBondDipoleExcludedFourier;

// Calculates the intra-molecular correction in Real-space for the skipped
// excluded-pairs in Fourier space. Note the use of the error-function
// instead of the error-function-complement.
int EwaldFourierIntraCorrectionAdsorbate(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr;
  POINT posA,posB;
  VECTOR dr;

  UAdsorbateAdsorbateChargeChargeExcludedFourier=0.0;

  if(OmitAdsorbateAdsorbateCoulombInteractions)
    return 0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
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

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            // note: no cutoff
            r=sqrt(rr);
            UAdsorbateAdsorbateChargeChargeExcludedFourier-=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)/r;
          }
        }
      }
    }
  }
  return 0;
}

void EwaldFourierIntraCorrectionCation(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr;
  POINT posA,posB;
  VECTOR dr;

  UCationCationChargeChargeExcludedFourier=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
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

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            // note: no cutoff
            r=sqrt(rr);
            UCationCationChargeChargeExcludedFourier-=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)/r;
          }
        }
      }
    }
  }
}

// Calculates the intra-molecular correction in Real-space for the skipped
// excluded-pairs in Fourier space. Note the use of the error-function
// instead of the error-function-complement.

void EwaldFourierIntraForceCorrectionChargeChargeFramework(void)
{
  int i,A,B,NumberOfExcludedPairs;
  int typeA,typeB;
  REAL ChargeA,ChargeB;
  REAL r,rr,temp;
  VECTOR posA,posB,dr;

  UHostHostChargeChargeExcludedFourier=0.0;
  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    if(Framework[CurrentSystem].FrameworkModels[CurrentFramework]==FLEXIBLE)
    {
      NumberOfExcludedPairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[CurrentFramework];
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        A=Framework[CurrentSystem].ExcludedIntraChargeCharge[CurrentFramework][i].A;
        B=Framework[CurrentSystem].ExcludedIntraChargeCharge[CurrentFramework][i].B;

        typeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
        typeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

        ChargeA=PseudoAtoms[typeA].Charge;
        ChargeB=PseudoAtoms[typeB].Charge;

        posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
        posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        switch(ChargeMethod)
        {
          case EWALD:
            if(r>1e-4)
            {
              UHostHostChargeChargeExcludedFourier-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                                    ChargeA*ChargeB/r;

              temp=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                   (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                   (r*rr);
            }
            else
            {
              UHostHostChargeChargeExcludedFourier-=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*2.0*Alpha[CurrentSystem]/sqrt(M_PI);
              temp=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*4.0*CUBE(Alpha[CurrentSystem])/(3.0*sqrt(M_PI));
            }

            Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.x+=temp*dr.x;
            Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.y+=temp*dr.y;
            Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.z+=temp*dr.z;

            Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.x-=temp*dr.x;
            Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.y-=temp*dr.y;
            Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.z-=temp*dr.z;

            StrainDerivativeTensor[CurrentSystem].ax-=temp*dr.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx-=temp*dr.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx-=temp*dr.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay-=temp*dr.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by-=temp*dr.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy-=temp*dr.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az-=temp*dr.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz-=temp*dr.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz-=temp*dr.z*dr.z;
            break;
        }
      }
    }
  }
}

void EwaldFourierIntraForceCorrectionChargeBondDipoleFramework(void)
{
  int i,A,B,B1,B2,NumberOfExcludedPairs,f1;
  int TypeA;
  REAL ChargeA;
  REAL r,rr,ri2,temp,fac1,fac2,cosB,energy,DipoleMagnitudeB,length;
  REAL Bt0,Bt1,Bt2;
  VECTOR posA,posB,posB1,posB2,dr,dipoleB;
  VECTOR term,fa1,fb1,fb2;
  REAL_MATRIX3x3 v;

  UHostHostChargeBondDipoleExcludedFourier=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      NumberOfExcludedPairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeBondDipole[f1];
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        // the index of the bonddipole is the second index 'B'
        A=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraChargeBondDipole[f1][i].B;

        TypeA=Framework[CurrentSystem].Atoms[f1][A].Type;
        ChargeA=PseudoAtoms[TypeA].Charge;
        posA=Framework[CurrentSystem].Atoms[f1][A].Position;

        DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][B];
        B1=Framework[CurrentSystem].BondDipoles[f1][B].A;
        B2=Framework[CurrentSystem].BondDipoles[f1][B].B;
        posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
        posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        ri2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeB/length;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posB.x-posA.x;
        dr.y=posB.y-posA.y;
        dr.z=posB.z-posA.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        Bt0=-ErrorFunction(Alpha[CurrentSystem]*r)/r;
        Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
            -ErrorFunction(Alpha[CurrentSystem]*r)/(rr*r);
        Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
            4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
            -3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*r);

        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
        UHostHostChargeBondDipoleExcludedFourier-=energy;

        term.x=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
        term.y=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
        term.z=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

        fa1.x=term.x;
        fa1.y=term.y;
        fa1.z=term.z;
        Framework[CurrentSystem].Atoms[f1][A].Force.x+=fa1.x;
        Framework[CurrentSystem].Atoms[f1][A].Force.y+=fa1.y;
        Framework[CurrentSystem].Atoms[f1][A].Force.z+=fa1.z;

        fac1=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*DipoleMagnitudeB/length;
        fac2=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*cosB/(DipoleMagnitudeB*length);

        fb1.x=0.5*term.x+fac1*dr.x-fac2*dipoleB.x;
        fb1.y=0.5*term.y+fac1*dr.y-fac2*dipoleB.y;
        fb1.z=0.5*term.z+fac1*dr.z-fac2*dipoleB.z;

        fb2.x=0.5*term.x-fac1*dr.x+fac2*dipoleB.x;
        fb2.y=0.5*term.y-fac1*dr.y+fac2*dipoleB.y;
        fb2.z=0.5*term.z-fac1*dr.z+fac2*dipoleB.z;

        Framework[CurrentSystem].Atoms[f1][B1].Force.x-=fb1.x;
        Framework[CurrentSystem].Atoms[f1][B1].Force.y-=fb1.y;
        Framework[CurrentSystem].Atoms[f1][B1].Force.z-=fb1.z;

        Framework[CurrentSystem].Atoms[f1][B2].Force.x-=fb2.x;
        Framework[CurrentSystem].Atoms[f1][B2].Force.y-=fb2.y;
        Framework[CurrentSystem].Atoms[f1][B2].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }
}

void EwaldFourierIntraForceCorrectionBondDipoleBondDipoleFramework(void)
{
  int i,A,B,A1,A2,B1,B2,NumberOfExcludedPairs,f1;
  REAL ri2,rk2;
  REAL r,rr,temp,cosA,cosB,cosAB,energy,DipoleMagnitudeA,DipoleMagnitudeB,length;
  REAL Bt0,Bt1,Bt2,Bt3;
  VECTOR posA,posB,posA1,posA2,posB1,posB2,dr,dipoleA,dipoleB;
  VECTOR term,fb1,fb2,fa1,fa2,termA,termB;
  REAL_MATRIX3x3 v;

  UHostHostBondDipoleBondDipoleExcludedFourier=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      NumberOfExcludedPairs=Framework[CurrentSystem].NumberOfExcludedIntraBondDipoleBondDipole[f1];
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        // the index of the first bonddipole is the first index 'A'
        A=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].A;
        B=Framework[CurrentSystem].ExcludedIntraBondDipoleBondDipole[f1][i].B;

        DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[f1][A];
        A1=Framework[CurrentSystem].BondDipoles[f1][A].A;
        A2=Framework[CurrentSystem].BondDipoles[f1][A].B;
        posA1=Framework[CurrentSystem].Atoms[f1][A1].Position;
        posA2=Framework[CurrentSystem].Atoms[f1][A2].Position;
        dipoleA.x=posA2.x-posA1.x;
        dipoleA.y=posA2.y-posA1.y;
        dipoleA.z=posA2.z-posA1.z;
        dipoleA=ApplyBoundaryCondition(dipoleA);
        posA.x=posA1.x+0.5*dipoleA.x;
        posA.y=posA1.y+0.5*dipoleA.y;
        posA.z=posA1.z+0.5*dipoleA.z;
        ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
        length=sqrt(ri2);
        temp=DipoleMagnitudeA/length;
        dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

        DipoleMagnitudeB=Framework[CurrentSystem].BondDipoleMagnitude[f1][B];
        B1=Framework[CurrentSystem].BondDipoles[f1][B].A;
        B2=Framework[CurrentSystem].BondDipoles[f1][B].B;
        posB1=Framework[CurrentSystem].Atoms[f1][B1].Position;
        posB2=Framework[CurrentSystem].Atoms[f1][B2].Position;
        dipoleB.x=posB2.x-posB1.x;
        dipoleB.y=posB2.y-posB1.y;
        dipoleB.z=posB2.z-posB1.z;
        dipoleB=ApplyBoundaryCondition(dipoleB);
        posB.x=posB1.x+0.5*dipoleB.x;
        posB.y=posB1.y+0.5*dipoleB.y;
        posB.z=posB1.z+0.5*dipoleB.z;
        rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
        length=sqrt(rk2);
        temp=DipoleMagnitudeB/length;
        dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        r=sqrt(rr);

        // ErrorFunctionComplement replaced by ErrorFunction, and sign reversed (taking the sign changes of 'Erf' into account)
        // in the derivative of 'Erf' instead of 'Erfc' all terms without 'Erf' changes sign
        Bt0=-ErrorFunction(Alpha[CurrentSystem]*r)/r;
        Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
            -ErrorFunction(Alpha[CurrentSystem]*r)/(rr*r);
        Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
            4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
            -3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*r);
        Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
            20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
            8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
            -15.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);

        cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
        cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
        cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
        energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
        UHostHostBondDipoleBondDipoleExcludedFourier+=energy;

        term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
        term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
        term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

        termA.x=COULOMBIC_CONVERSION_FACTOR*(
                Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
                -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
                +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

        termA.y=COULOMBIC_CONVERSION_FACTOR*(
                Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
                -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
                +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

        termA.z=COULOMBIC_CONVERSION_FACTOR*(
                Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
                -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
                +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

        termB.x=COULOMBIC_CONVERSION_FACTOR*(
                Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
                -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
                +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

        termB.y=COULOMBIC_CONVERSION_FACTOR*(
                Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
                -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
                +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

        termB.z=COULOMBIC_CONVERSION_FACTOR*(
                Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
                -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
                +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

        fa1.x=0.5*term.x+termA.x;
        fa1.y=0.5*term.y+termA.y;
        fa1.z=0.5*term.z+termA.z;
        fa2.x=0.5*term.x-termA.x;
        fa2.y=0.5*term.y-termA.y;
        fa2.z=0.5*term.z-termA.z;

        fb1.x=-0.5*term.x+termB.x;
        fb1.y=-0.5*term.y+termB.y;
        fb1.z=-0.5*term.z+termB.z;
        fb2.x=-0.5*term.x-termB.x;
        fb2.y=-0.5*term.y-termB.y;
        fb2.z=-0.5*term.z-termB.z;

        Framework[CurrentSystem].Atoms[f1][A1].Force.x-=fa1.x;
        Framework[CurrentSystem].Atoms[f1][A1].Force.y-=fa1.y;
        Framework[CurrentSystem].Atoms[f1][A1].Force.z-=fa1.z;

        Framework[CurrentSystem].Atoms[f1][A2].Force.x-=fa2.x;
        Framework[CurrentSystem].Atoms[f1][A2].Force.y-=fa2.y;
        Framework[CurrentSystem].Atoms[f1][A2].Force.z-=fa2.z;

        Framework[CurrentSystem].Atoms[f1][B1].Force.x-=fb1.x;
        Framework[CurrentSystem].Atoms[f1][B1].Force.y-=fb1.y;
        Framework[CurrentSystem].Atoms[f1][B1].Force.z-=fb1.z;

        Framework[CurrentSystem].Atoms[f1][B2].Force.x-=fb2.x;
        Framework[CurrentSystem].Atoms[f1][B2].Force.y-=fb2.y;
        Framework[CurrentSystem].Atoms[f1][B2].Force.z-=fb2.z;

        // convert forces on atoms to molecular virial
        v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
        v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
        v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

        v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
        v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
        v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

        v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
        v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
        v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

        // the strain derivative
        StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
        StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
        StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

        StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
        StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
        StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

        StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
        StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
        StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
      }
    }
  }
}

int EwaldFourierIntraForceCorrectionChargeChargeAdsorbate(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr,temp;
  POINT posA,posB;
  VECTOR dr;

  UAdsorbateAdsorbateChargeChargeExcludedFourier=0.0;

  if(OmitAdsorbateAdsorbateCoulombInteractions)
    return 0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
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

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            // note: no cutoff
            UAdsorbateAdsorbateChargeChargeExcludedFourier-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                          ChargeA*ChargeB/r;

            temp=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                 (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                 (r*rr);

            Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=temp*dr.x;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=temp*dr.y;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=temp*dr.z;

            Adsorbates[CurrentSystem][m].Atoms[B].Force.x-=temp*dr.x;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.y-=temp*dr.y;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.z-=temp*dr.z;

            StrainDerivativeTensor[CurrentSystem].ax-=temp*dr.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx-=temp*dr.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx-=temp*dr.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay-=temp*dr.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by-=temp*dr.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy-=temp*dr.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az-=temp*dr.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz-=temp*dr.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz-=temp*dr.z*dr.z;
          }
        }
      }
    }
  }
  return 0;
}

int EwaldFourierIntraForceCorrectionChargeChargeCation(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr,temp;
  POINT posA,posB;
  VECTOR dr;

  UCationCationChargeChargeExcludedFourier=0.0;

  if(OmitCationCationCoulombInteractions)
    return 0;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
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

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            // note: no cutoff
            UCationCationChargeChargeExcludedFourier-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                          ChargeA*ChargeB/r;

            temp=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                 (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                 (r*rr);

            Cations[CurrentSystem][m].Atoms[A].Force.x+=temp*dr.x;
            Cations[CurrentSystem][m].Atoms[A].Force.y+=temp*dr.y;
            Cations[CurrentSystem][m].Atoms[A].Force.z+=temp*dr.z;

            Cations[CurrentSystem][m].Atoms[B].Force.x-=temp*dr.x;
            Cations[CurrentSystem][m].Atoms[B].Force.y-=temp*dr.y;
            Cations[CurrentSystem][m].Atoms[B].Force.z-=temp*dr.z;

            StrainDerivativeTensor[CurrentSystem].ax-=temp*dr.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx-=temp*dr.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx-=temp*dr.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay-=temp*dr.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by-=temp*dr.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy-=temp*dr.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az-=temp*dr.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz-=temp*dr.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz-=temp*dr.z*dr.z;
          }
        }
      }
    }
  }
  return 0;
}

int EwaldFourierIntraForceCorrectionChargeBondDipoleAdsorbate(void)
{
  int i,m,NumberOfExcludedPairs,Type;
  int A,B,B1,B2;
  REAL r,rr,temp,ChargeA,fac1,fac2;
  POINT posA,posB,posB1,posB2;
  VECTOR dr,dipoleB;
  REAL ri2,length,cosB,energy;
  REAL DipoleMagnitudeB;
  REAL Bt0,Bt1,Bt2;
  VECTOR term,fa1,fb1,fb2;
  REAL_MATRIX3x3 v;

  UAdsorbateAdsorbateChargeBondDipoleExcludedFourier=0.0;

  //if(OmitAdsorbateAdsorbateCoulombInteractions)
  //  return 0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<NumberOfExcludedPairs;i++)
    {
      A=Components[Type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[Type].ExcludedIntraChargeBondDipole[i].B;

      ChargeA=PseudoAtoms[Components[Type].Type[A]].Charge;
      posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
      posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      ri2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);
      switch(ChargeMethod)
      {
        case NONE:
          Bt0=Bt1=Bt2=0.0;
          break;
        case EWALD:
          Bt0=-ErrorFunction(Alpha[CurrentSystem]*r)/r;
          Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
              -ErrorFunction(Alpha[CurrentSystem]*r)/(rr*r);
          Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
              4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
              -3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*r);
          break;
        default:
          Bt0=Bt1=Bt2=0.0;
          break;
      }

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
      UAdsorbateAdsorbateChargeBondDipoleExcludedFourier-=energy;

      term.x=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
      term.y=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
      term.z=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

      fa1.x=term.x;
      fa1.y=term.y;
      fa1.z=term.z;
      Adsorbates[CurrentSystem][m].Atoms[A].Force.x+=fa1.x;
      Adsorbates[CurrentSystem][m].Atoms[A].Force.y+=fa1.y;
      Adsorbates[CurrentSystem][m].Atoms[A].Force.z+=fa1.z;

      fac1=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*DipoleMagnitudeB/length;
      fac2=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*cosB/(DipoleMagnitudeB*length);

      fb1.x=0.5*term.x+fac1*dr.x-fac2*dipoleB.x;
      fb1.y=0.5*term.y+fac1*dr.y-fac2*dipoleB.y;
      fb1.z=0.5*term.z+fac1*dr.z-fac2*dipoleB.z;

      fb2.x=0.5*term.x-fac1*dr.x+fac2*dipoleB.x;
      fb2.y=0.5*term.y-fac1*dr.y+fac2*dipoleB.y;
      fb2.z=0.5*term.z-fac1*dr.z+fac2*dipoleB.z;


      Adsorbates[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
      Adsorbates[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
      Adsorbates[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

      Adsorbates[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
      Adsorbates[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
      Adsorbates[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;


      // convert forces on atoms to molecular virial
      v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
      v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
      v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

      v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
      v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
      v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

      v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
      v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
      v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

      // the strain derivative
      StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
      StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
      StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

      StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
      StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
      StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

      StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
      StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
      StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;

    }
  }
  return 0;
}

int EwaldFourierIntraForceCorrectionChargeBondDipoleCation(void)
{
  int i,m,NumberOfExcludedPairs,Type;
  int A,B,B1,B2;
  REAL r,rr,temp,ChargeA,fac1,fac2;
  POINT posA,posB,posB1,posB2;
  VECTOR dr,dipoleB;
  REAL ri2,length,cosB,energy;
  REAL DipoleMagnitudeB;
  REAL Bt0,Bt1,Bt2;
  VECTOR term,fa1,fb1,fb2;
  REAL_MATRIX3x3 v;

  UCationCationChargeBondDipoleExcludedFourier=0.0;

  //if(OmitAdsorbateAdsorbateCoulombInteractions)
  //  return 0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeBondDipole;
    for(i=0;i<NumberOfExcludedPairs;i++)
    {
      A=Components[Type].ExcludedIntraChargeBondDipole[i].A;
      B=Components[Type].ExcludedIntraChargeBondDipole[i].B;

      ChargeA=PseudoAtoms[Components[Type].Type[A]].Charge;
      posA=Cations[CurrentSystem][m].Atoms[A].Position;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
      posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      ri2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posB.x-posA.x;
      dr.y=posB.y-posA.y;
      dr.z=posB.z-posA.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);
      switch(ChargeMethod)
      {
        case NONE:
          Bt0=Bt1=Bt2=0.0;
          break;
        case EWALD:
          Bt0=-ErrorFunction(Alpha[CurrentSystem]*r)/r;
          Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
              -ErrorFunction(Alpha[CurrentSystem]*r)/(rr*r);
          Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
              4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
              -3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*r);
          break;
        default:
          Bt0=Bt1=Bt2=0.0;
          break;
      }

      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*ChargeA*cosB);
      UCationCationChargeBondDipoleExcludedFourier-=energy;

      term.x=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
      term.y=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
      term.z=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

      fa1.x=term.x;
      fa1.y=term.y;
      fa1.z=term.z;
      Cations[CurrentSystem][m].Atoms[A].Force.x+=fa1.x;
      Cations[CurrentSystem][m].Atoms[A].Force.y+=fa1.y;
      Cations[CurrentSystem][m].Atoms[A].Force.z+=fa1.z;

      fac1=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*DipoleMagnitudeB/length;
      fac2=COULOMBIC_CONVERSION_FACTOR*ChargeA*Bt1*cosB/(DipoleMagnitudeB*length);

      fb1.x=0.5*term.x+fac1*dr.x-fac2*dipoleB.x;
      fb1.y=0.5*term.y+fac1*dr.y-fac2*dipoleB.y;
      fb1.z=0.5*term.z+fac1*dr.z-fac2*dipoleB.z;

      fb2.x=0.5*term.x-fac1*dr.x+fac2*dipoleB.x;
      fb2.y=0.5*term.y-fac1*dr.y+fac2*dipoleB.y;
      fb2.z=0.5*term.z-fac1*dr.z+fac2*dipoleB.z;


      Cations[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
      Cations[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
      Cations[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

      Cations[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
      Cations[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
      Cations[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;

      // convert forces on atoms to molecular virial
      v.ax=fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
      v.bx=fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
      v.cx=fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

      v.ay=fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
      v.by=fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
      v.cy=fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

      v.az=fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
      v.bz=fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
      v.cz=fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

      // the strain derivative
      StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
      StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
      StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

      StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
      StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
      StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

      StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
      StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
      StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;

    }
  }
  return 0;
}


int EwaldFourierIntraForceCorrectionBondDipoleBondDipoleAdsorbate(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type;
  int A1,A2,B1,B2;
  REAL r,rr,temp;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  VECTOR dr,dipoleA,dipoleB;
  REAL ri2,rk2,length,cosAB,cosA,cosB,energy;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL Bt0,Bt1,Bt2,Bt3;
  VECTOR term,fb1,fb2,fa1,fa2,termA,termB;
  REAL_MATRIX3x3 v;


  UAdsorbateAdsorbateBondDipoleBondDipoleExcludedFourier=0.0;

  if(OmitAdsorbateAdsorbateCoulombInteractions)
    return 0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<NumberOfExcludedPairs;i++)
    {
      A=Components[Type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[Type].ExcludedIntraBondDipoleBondDipole[i].B;

      A1=Components[Type].BondDipoles[A].A;
      A2=Components[Type].BondDipoles[A].B;
      DipoleMagnitudeA=Components[Type].BondDipoleMagnitude[A];
      posA1=Adsorbates[CurrentSystem][m].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][m].Atoms[A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      posB1=Adsorbates[CurrentSystem][m].Atoms[B1].Position;
      posB2=Adsorbates[CurrentSystem][m].Atoms[B2].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(rk2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      // ErrorFunctionComplement replaced by ErrorFunction, and sign reversed (taking the sign changes of 'Erf' into account)
      // in the derivative of 'Erf' instead of 'Erfc' all terms without 'Erf' changes sign
      Bt0=-ErrorFunction(Alpha[CurrentSystem]*r)/r;
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -ErrorFunction(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*r);
      Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
          20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -15.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
      UAdsorbateAdsorbateBondDipoleBondDipoleExcludedFourier+=energy;

      term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
      term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
      term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

      termA.x=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
              +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

      termA.y=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
              +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

      termA.z=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
              +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

      termB.x=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
              +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

      termB.y=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
              +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

      termB.z=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
              +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

      fa1.x=0.5*term.x+termA.x;
      fa1.y=0.5*term.y+termA.y;
      fa1.z=0.5*term.z+termA.z;
      fa2.x=0.5*term.x-termA.x;
      fa2.y=0.5*term.y-termA.y;
      fa2.z=0.5*term.z-termA.z;

      fb1.x=-0.5*term.x+termB.x;
      fb1.y=-0.5*term.y+termB.y;
      fb1.z=-0.5*term.z+termB.z;
      fb2.x=-0.5*term.x-termB.x;
      fb2.y=-0.5*term.y-termB.y;
      fb2.z=-0.5*term.z-termB.z;

      Adsorbates[CurrentSystem][m].Atoms[A1].Force.x-=fa1.x;
      Adsorbates[CurrentSystem][m].Atoms[A1].Force.y-=fa1.y;
      Adsorbates[CurrentSystem][m].Atoms[A1].Force.z-=fa1.z;

      Adsorbates[CurrentSystem][m].Atoms[A2].Force.x-=fa2.x;
      Adsorbates[CurrentSystem][m].Atoms[A2].Force.y-=fa2.y;
      Adsorbates[CurrentSystem][m].Atoms[A2].Force.z-=fa2.z;

      Adsorbates[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
      Adsorbates[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
      Adsorbates[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

      Adsorbates[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
      Adsorbates[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
      Adsorbates[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;

      // convert forces on atoms to molecular virial
      v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
      v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
      v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

      v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
      v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
      v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

      v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
      v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
      v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

      // the strain derivative
      StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
      StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
      StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

      StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
      StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
      StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

      StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
      StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
      StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
    }
  }
  return 0;
}

int EwaldFourierIntraForceCorrectionBondDipoleBondDipoleCation(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type;
  int A1,A2,B1,B2;
  REAL r,rr,temp;
  POINT posA,posB,posA1,posA2,posB1,posB2;
  VECTOR dr,dipoleA,dipoleB;
  REAL ri2,rk2,length,cosAB,cosA,cosB,energy;
  REAL DipoleMagnitudeA,DipoleMagnitudeB;
  REAL Bt0,Bt1,Bt2,Bt3;
  VECTOR term,fb1,fb2,fa1,fa2,termA,termB;
  REAL_MATRIX3x3 v;


  UCationCationBondDipoleBondDipoleExcludedFourier=0.0;

  if(OmitCationCationCoulombInteractions)
    return 0;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraBondDipoleBondDipole;
    for(i=0;i<NumberOfExcludedPairs;i++)
    {
      A=Components[Type].ExcludedIntraBondDipoleBondDipole[i].A;
      B=Components[Type].ExcludedIntraBondDipoleBondDipole[i].B;

      A1=Components[Type].BondDipoles[A].A;
      A2=Components[Type].BondDipoles[A].B;
      DipoleMagnitudeA=Components[Type].BondDipoleMagnitude[A];
      posA1=Cations[CurrentSystem][m].Atoms[A1].Position;
      posA2=Cations[CurrentSystem][m].Atoms[A2].Position;
      dipoleA.x=posA2.x-posA1.x;
      dipoleA.y=posA2.y-posA1.y;
      dipoleA.z=posA2.z-posA1.z;
      posA.x=posA1.x+0.5*dipoleA.x;
      posA.y=posA1.y+0.5*dipoleA.y;
      posA.z=posA1.z+0.5*dipoleA.z;
      ri2=SQR(dipoleA.x)+SQR(dipoleA.y)+SQR(dipoleA.z);
      length=sqrt(ri2);
      temp=DipoleMagnitudeA/length;
      dipoleA.x*=temp; dipoleA.y*=temp; dipoleA.z*=temp;

      B1=Components[Type].BondDipoles[B].A;
      B2=Components[Type].BondDipoles[B].B;
      DipoleMagnitudeB=Components[Type].BondDipoleMagnitude[B];
      posB1=Cations[CurrentSystem][m].Atoms[B1].Position;
      posB2=Cations[CurrentSystem][m].Atoms[B2].Position;
      dipoleB.x=posB2.x-posB1.x;
      dipoleB.y=posB2.y-posB1.y;
      dipoleB.z=posB2.z-posB1.z;
      posB.x=posB1.x+0.5*dipoleB.x;
      posB.y=posB1.y+0.5*dipoleB.y;
      posB.z=posB1.z+0.5*dipoleB.z;
      rk2=SQR(dipoleB.x)+SQR(dipoleB.y)+SQR(dipoleB.z);
      length=sqrt(rk2);
      temp=DipoleMagnitudeB/length;
      dipoleB.x*=temp; dipoleB.y*=temp; dipoleB.z*=temp;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      // ErrorFunctionComplement replaced by ErrorFunction, and sign reversed (taking the sign changes of 'Erf' into account)
      // in the derivative of 'Erf' instead of 'Erfc' all terms without 'Erf' changes sign
      Bt0=-ErrorFunction(Alpha[CurrentSystem]*r)/r;
      Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -ErrorFunction(Alpha[CurrentSystem]*r)/(rr*r);
      Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*r);
      Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
          20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
          8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
          -15.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);

      cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
      cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
      cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;
      energy=COULOMBIC_CONVERSION_FACTOR*(Bt1*cosAB-Bt2*cosA*cosB);
      UCationCationBondDipoleBondDipoleExcludedFourier+=energy;

      term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
      term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
      term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

      termA.x=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt2*cosA*cosB*dipoleA.x/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt1*DipoleMagnitudeA*dipoleB.x/sqrt(ri2)
              +Bt2*DipoleMagnitudeA*cosB*dr.x/sqrt(ri2));

      termA.y=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt2*cosA*cosB*dipoleA.y/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt1*DipoleMagnitudeA*dipoleB.y/sqrt(ri2)
              +Bt2*DipoleMagnitudeA*cosB*dr.y/sqrt(ri2));

      termA.z=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt2*cosA*cosB*dipoleA.z/(sqrt(ri2)*DipoleMagnitudeA)
              -Bt1*DipoleMagnitudeA*dipoleB.z/sqrt(ri2)
              +Bt2*DipoleMagnitudeA*cosB*dr.z/sqrt(ri2));

      termB.x=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt2*cosA*cosB*dipoleB.x/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt1*DipoleMagnitudeB*dipoleA.x/sqrt(rk2)
              +Bt2*DipoleMagnitudeB*cosA*dr.x/sqrt(rk2));

      termB.y=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt2*cosA*cosB*dipoleB.y/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt1*DipoleMagnitudeB*dipoleA.y/sqrt(rk2)
              +Bt2*DipoleMagnitudeB*cosA*dr.y/sqrt(rk2));

      termB.z=COULOMBIC_CONVERSION_FACTOR*(
              Bt1*cosAB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt2*cosA*cosB*dipoleB.z/(sqrt(rk2)*DipoleMagnitudeB)
              -Bt1*DipoleMagnitudeB*dipoleA.z/sqrt(rk2)
              +Bt2*DipoleMagnitudeB*cosA*dr.z/sqrt(rk2));

      fa1.x=0.5*term.x+termA.x;
      fa1.y=0.5*term.y+termA.y;
      fa1.z=0.5*term.z+termA.z;
      fa2.x=0.5*term.x-termA.x;
      fa2.y=0.5*term.y-termA.y;
      fa2.z=0.5*term.z-termA.z;

      fb1.x=-0.5*term.x+termB.x;
      fb1.y=-0.5*term.y+termB.y;
      fb1.z=-0.5*term.z+termB.z;
      fb2.x=-0.5*term.x-termB.x;
      fb2.y=-0.5*term.y-termB.y;
      fb2.z=-0.5*term.z-termB.z;

      Cations[CurrentSystem][m].Atoms[A1].Force.x-=fa1.x;
      Cations[CurrentSystem][m].Atoms[A1].Force.y-=fa1.y;
      Cations[CurrentSystem][m].Atoms[A1].Force.z-=fa1.z;

      Cations[CurrentSystem][m].Atoms[A2].Force.x-=fa2.x;
      Cations[CurrentSystem][m].Atoms[A2].Force.y-=fa2.y;
      Cations[CurrentSystem][m].Atoms[A2].Force.z-=fa2.z;

      Cations[CurrentSystem][m].Atoms[B1].Force.x-=fb1.x;
      Cations[CurrentSystem][m].Atoms[B1].Force.y-=fb1.y;
      Cations[CurrentSystem][m].Atoms[B1].Force.z-=fb1.z;

      Cations[CurrentSystem][m].Atoms[B2].Force.x-=fb2.x;
      Cations[CurrentSystem][m].Atoms[B2].Force.y-=fb2.y;
      Cations[CurrentSystem][m].Atoms[B2].Force.z-=fb2.z;

      // convert forces on atoms to molecular virial
      v.ax=fa1.x*(posA1.x-posA.x)+fa2.x*(posA2.x-posA.x)+fb1.x*(posB1.x-posB.x)+fb2.x*(posB2.x-posB.x);
      v.bx=fa1.y*(posA1.x-posA.x)+fa2.y*(posA2.x-posA.x)+fb1.y*(posB1.x-posB.x)+fb2.y*(posB2.x-posB.x);
      v.cx=fa1.z*(posA1.x-posA.x)+fa2.z*(posA2.x-posA.x)+fb1.z*(posB1.x-posB.x)+fb2.z*(posB2.x-posB.x);

      v.ay=fa1.x*(posA1.y-posA.y)+fa2.x*(posA2.y-posA.y)+fb1.x*(posB1.y-posB.y)+fb2.x*(posB2.y-posB.y);
      v.by=fa1.y*(posA1.y-posA.y)+fa2.y*(posA2.y-posA.y)+fb1.y*(posB1.y-posB.y)+fb2.y*(posB2.y-posB.y);
      v.cy=fa1.z*(posA1.y-posA.y)+fa2.z*(posA2.y-posA.y)+fb1.z*(posB1.y-posB.y)+fb2.z*(posB2.y-posB.y);

      v.az=fa1.x*(posA1.z-posA.z)+fa2.x*(posA2.z-posA.z)+fb1.x*(posB1.z-posB.z)+fb2.x*(posB2.z-posB.z);
      v.bz=fa1.y*(posA1.z-posA.z)+fa2.y*(posA2.z-posA.z)+fb1.y*(posB1.z-posB.z)+fb2.y*(posB2.z-posB.z);
      v.cz=fa1.z*(posA1.z-posA.z)+fa2.z*(posA2.z-posA.z)+fb1.z*(posB1.z-posB.z)+fb2.z*(posB2.z-posB.z);

      // the strain derivative
      StrainDerivativeTensor[CurrentSystem].ax+=term.x*dr.x+v.ax;
      StrainDerivativeTensor[CurrentSystem].bx+=term.y*dr.x+v.bx;
      StrainDerivativeTensor[CurrentSystem].cx+=term.z*dr.x+v.cx;

      StrainDerivativeTensor[CurrentSystem].ay+=term.x*dr.y+v.ay;
      StrainDerivativeTensor[CurrentSystem].by+=term.y*dr.y+v.by;
      StrainDerivativeTensor[CurrentSystem].cy+=term.z*dr.y+v.cy;

      StrainDerivativeTensor[CurrentSystem].az+=term.x*dr.z+v.az;
      StrainDerivativeTensor[CurrentSystem].bz+=term.y*dr.z+v.bz;
      StrainDerivativeTensor[CurrentSystem].cz+=term.z*dr.z+v.cz;
    }
  }
  return 0;
}


// Born-terms

void EwaldFourierIntraForceCorrectionFrameworkBornTerm(void)
{
  int i,A,B,NumberOfExcludedPairs;
  int typeA,typeB;
  REAL ChargeA,ChargeB;
  REAL r,rr,temp,DDF;
  VECTOR posA,posB,dr;

  UHostHostChargeChargeExcludedFourier=0.0;
  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    NumberOfExcludedPairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[CurrentFramework];
    for(i=0;i<NumberOfExcludedPairs;i++)
    {
      A=Framework[CurrentSystem].ExcludedIntraChargeCharge[CurrentFramework][i].A;
      B=Framework[CurrentSystem].ExcludedIntraChargeCharge[CurrentFramework][i].B;

      typeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
      typeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

      ChargeA=PseudoAtoms[typeA].Charge;
      ChargeB=PseudoAtoms[typeB].Charge;

      posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
      posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;

      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(rr);

      if(r>1e-4)
      {
        UHostHostChargeChargeExcludedFourier-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                                  ChargeA*ChargeB/r;

        temp=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
             (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
             (r*rr);

        DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
            (((-3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(r*rr))+
            (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
            (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));
      }
      else
      {
        UHostHostChargeChargeExcludedFourier-=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*2.0*Alpha[CurrentSystem]/sqrt(M_PI);
        temp=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*4.0*CUBE(Alpha[CurrentSystem])/(3.0*sqrt(M_PI));
        DDF=-COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*8.0*CUBE(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])/(5.0*sqrt(M_PI));
      }


      StrainDerivativeTensor[CurrentSystem].ax-=temp*dr.x*dr.x;
      StrainDerivativeTensor[CurrentSystem].bx-=temp*dr.y*dr.x;
      StrainDerivativeTensor[CurrentSystem].cx-=temp*dr.z*dr.x;

      StrainDerivativeTensor[CurrentSystem].ay-=temp*dr.x*dr.y;
      StrainDerivativeTensor[CurrentSystem].by-=temp*dr.y*dr.y;
      StrainDerivativeTensor[CurrentSystem].cy-=temp*dr.z*dr.y;

      StrainDerivativeTensor[CurrentSystem].az-=temp*dr.x*dr.z;
      StrainDerivativeTensor[CurrentSystem].bz-=temp*dr.y*dr.z;
      StrainDerivativeTensor[CurrentSystem].cz-=temp*dr.z*dr.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.x+=temp*dr.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.y+=temp*dr.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][A].Force.z+=temp*dr.z;

      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.x-=temp*dr.x;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.y-=temp*dr.y;
      Framework[CurrentSystem].Atoms[CurrentFramework][B].Force.z-=temp*dr.z;

      // add contribution to the born term
      AddContributionToBornTerm(DDF,dr);
    }
  }
}

void EwaldFourierIntraForceCorrectionAdsorbateBornTerm(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr,DF,DDF;
  POINT posA,posB;
  VECTOR dr;

  UAdsorbateAdsorbateChargeChargeExcludedFourier=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
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

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            // note: no cutoff
            UAdsorbateAdsorbateChargeChargeExcludedFourier-=COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)*
                          ChargeA*ChargeB/r;

            DF=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
               (ErrorFunction(Alpha[CurrentSystem]*r)-2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
               (r*rr);

            DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                   (((-3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(r*rr))+
                   (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                   (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));

            Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=DF*dr.x;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=DF*dr.y;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=DF*dr.z;

            Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=DF*dr.x;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=DF*dr.y;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=DF*dr.z;

            StrainDerivativeTensor[CurrentSystem].ax+=DF*dr.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx+=DF*dr.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx+=DF*dr.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay+=DF*dr.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by+=DF*dr.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy+=DF*dr.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az+=DF*dr.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz+=DF*dr.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz+=DF*dr.z*dr.z;

            // add contribution to the born term
            AddContributionToBornTerm(DDF,dr);
          }
        }
      }
    }
  }
}

void EwaldFourierIntraForceCorrectionCationBornTerm(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr,temp,DDF;
  POINT posA,posB;
  VECTOR dr;

  UCationCationChargeChargeExcludedFourier=0.0;
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
      NumberOfExcludedPairs=Components[Cations[CurrentSystem][m].Type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        Type=Cations[CurrentSystem][m].Type;
        A=Components[Type].ExcludedIntraChargeCharge[i].A;
        B=Components[Type].ExcludedIntraChargeCharge[i].B;

        TypeA=Cations[CurrentSystem][m].Atoms[A].Type;
        if(PseudoAtoms[Type].HasCharges)
        {
          TypeB=Cations[CurrentSystem][m].Atoms[B].Type;
          if(PseudoAtoms[TypeB].HasCharges)
          {
            ChargeA=PseudoAtoms[TypeA].Charge;
            ChargeB=PseudoAtoms[TypeB].Charge;
            posA=Cations[CurrentSystem][m].Atoms[A].Position;
            posB=Cations[CurrentSystem][m].Atoms[B].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

            // note: no cutoff
            r=sqrt(rr);
            UCationCationChargeChargeExcludedFourier-=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*ErrorFunction(Alpha[CurrentSystem]*r)/r;

            temp=COULOMBIC_CONVERSION_FACTOR*ChargeA*ChargeB*
                 (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                 (r*rr);

            DDF=ChargeA*ChargeB*COULOMBIC_CONVERSION_FACTOR*
                   (((-3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(r*rr))+
                   (4.0*CUBE(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))+
                   (6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem])*rr)/(sqrt(M_PI)*rr)))/(rr));

            Cations[CurrentSystem][m].Atoms[A].Force.x+=temp*dr.x;
            Cations[CurrentSystem][m].Atoms[A].Force.y+=temp*dr.y;
            Cations[CurrentSystem][m].Atoms[A].Force.z+=temp*dr.z;

            Cations[CurrentSystem][m].Atoms[B].Force.x-=temp*dr.x;
            Cations[CurrentSystem][m].Atoms[B].Force.y-=temp*dr.y;
            Cations[CurrentSystem][m].Atoms[B].Force.z-=temp*dr.z;

            StrainDerivativeTensor[CurrentSystem].ax-=temp*dr.x*dr.x;
            StrainDerivativeTensor[CurrentSystem].bx-=temp*dr.y*dr.x;
            StrainDerivativeTensor[CurrentSystem].cx-=temp*dr.z*dr.x;

            StrainDerivativeTensor[CurrentSystem].ay-=temp*dr.x*dr.y;
            StrainDerivativeTensor[CurrentSystem].by-=temp*dr.y*dr.y;
            StrainDerivativeTensor[CurrentSystem].cy-=temp*dr.z*dr.y;

            StrainDerivativeTensor[CurrentSystem].az-=temp*dr.x*dr.z;
            StrainDerivativeTensor[CurrentSystem].bz-=temp*dr.y*dr.z;
            StrainDerivativeTensor[CurrentSystem].cz-=temp*dr.z*dr.z;

            // add contribution to the born term
            AddContributionToBornTerm(DDF,dr);
          }
        }
      }
    }
  }
}

void EwaldFourierIntraElectricFieldCorrectionChargeChargeFramework(void)
{
  int i,A,B,NumberOfExcludedPairs,f1;
  int typeA,typeB;
  REAL ChargeA,ChargeB;
  REAL r,rr,temp;
  REAL force_factor_A,force_factor_B;
  VECTOR posA,posB,dr;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    NumberOfExcludedPairs=Framework[CurrentSystem].NumberOfExcludedIntraChargeCharge[f1];
    for(i=0;i<NumberOfExcludedPairs;i++)
    {
      A=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i].A;
      B=Framework[CurrentSystem].ExcludedIntraChargeCharge[f1][i].B;

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

      temp=COULOMBIC_CONVERSION_FACTOR*
           (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
           (r*rr);
      force_factor_A=temp*ChargeB;
      force_factor_B=temp*ChargeA;

      Framework[CurrentSystem].Atoms[f1][A].ElectricField.x+=force_factor_A*dr.x;
      Framework[CurrentSystem].Atoms[f1][A].ElectricField.y+=force_factor_A*dr.y;
      Framework[CurrentSystem].Atoms[f1][A].ElectricField.z+=force_factor_A*dr.z;

      Framework[CurrentSystem].Atoms[f1][B].ElectricField.x-=force_factor_B*dr.x;
      Framework[CurrentSystem].Atoms[f1][B].ElectricField.y-=force_factor_B*dr.y;
      Framework[CurrentSystem].Atoms[f1][B].ElectricField.z-=force_factor_B*dr.z;
    }
  }
}

int EwaldFourierIntraElectricFieldCorrectionChargeChargeAdsorbate(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr,temp;
  REAL force_factor_A,force_factor_B;
  POINT posA,posB;
  VECTOR dr;

  if(OmitAdsorbateAdsorbateCoulombInteractions)
    return 0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
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

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            temp=COULOMBIC_CONVERSION_FACTOR*
                 (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                 (r*rr);
            force_factor_A=temp*ChargeB;
            force_factor_B=temp*ChargeA;

            Adsorbates[CurrentSystem][m].Atoms[A].ElectricField.x+=force_factor_A*dr.x;
            Adsorbates[CurrentSystem][m].Atoms[A].ElectricField.y+=force_factor_A*dr.y;
            Adsorbates[CurrentSystem][m].Atoms[A].ElectricField.z+=force_factor_A*dr.z;

            Adsorbates[CurrentSystem][m].Atoms[B].ElectricField.x-=force_factor_B*dr.x;
            Adsorbates[CurrentSystem][m].Atoms[B].ElectricField.y-=force_factor_B*dr.y;
            Adsorbates[CurrentSystem][m].Atoms[B].ElectricField.z-=force_factor_B*dr.z;
          }
        }
      }
    }
  }
  return 0;
}

int EwaldFourierIntraElectricFieldCorrectionChargeChargeCation(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL force_factor_A,force_factor_B;
  REAL r,rr,temp;
  POINT posA,posB;
  VECTOR dr;

  if(OmitCationCationCoulombInteractions)
    return 0;

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
      NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        A=Components[Type].ExcludedIntraChargeCharge[i].A;
        B=Components[Type].ExcludedIntraChargeCharge[i].B;

        TypeA=Cations[CurrentSystem][m].Atoms[A].Type;
        ChargeA=PseudoAtoms[TypeA].Charge;
        if(PseudoAtoms[TypeA].HasCharges)
        {
          TypeB=Cations[CurrentSystem][m].Atoms[B].Type;
          if(PseudoAtoms[TypeB].HasCharges)
          {
            ChargeB=PseudoAtoms[TypeB].Charge;
            posA=Cations[CurrentSystem][m].Atoms[A].Position;
            posB=Cations[CurrentSystem][m].Atoms[B].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            temp=COULOMBIC_CONVERSION_FACTOR*
                 (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                 (r*rr);
            force_factor_A=temp*ChargeB;
            force_factor_B=temp*ChargeA;

            Cations[CurrentSystem][m].Atoms[A].ElectricField.x+=force_factor_A*dr.x;
            Cations[CurrentSystem][m].Atoms[A].ElectricField.y+=force_factor_A*dr.y;
            Cations[CurrentSystem][m].Atoms[A].ElectricField.z+=force_factor_A*dr.z;

            Cations[CurrentSystem][m].Atoms[B].ElectricField.x-=force_factor_B*dr.x;
            Cations[CurrentSystem][m].Atoms[B].ElectricField.y-=force_factor_B*dr.y;
            Cations[CurrentSystem][m].Atoms[B].ElectricField.z-=force_factor_B*dr.z;
          }
        }
      }
    }
  }
  return 0;
}

int EwaldFourierIntraElectricFieldCorrectionChargeChargeMC(int New, int excl_ads,int excl_cation)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr,temp;
  REAL force_factor_A,force_factor_B;
  POINT posA,posB;
  VECTOR dr;

  if(OmitAdsorbateAdsorbateCoulombInteractions)
    return 0;

  if(!OmitAdsorbateAdsorbatePolarization)
  {
    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      if(m!=excl_ads)
      {
        Type=Adsorbates[CurrentSystem][m].Type;
        if(Components[Type].HasCharges)
        {
          energy_sum=0.0;
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

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(rr);

                temp=COULOMBIC_CONVERSION_FACTOR*
                     (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                     (r*rr);
                force_factor_A=temp*ChargeB;
                force_factor_B=temp*ChargeA;

                Adsorbates[CurrentSystem][m].Atoms[A].ElectricField.x+=force_factor_A*dr.x;
                Adsorbates[CurrentSystem][m].Atoms[A].ElectricField.y+=force_factor_A*dr.y;
                Adsorbates[CurrentSystem][m].Atoms[A].ElectricField.z+=force_factor_A*dr.z;

                Adsorbates[CurrentSystem][m].Atoms[B].ElectricField.x-=force_factor_B*dr.x;
                Adsorbates[CurrentSystem][m].Atoms[B].ElectricField.y-=force_factor_B*dr.y;
                Adsorbates[CurrentSystem][m].Atoms[B].ElectricField.z-=force_factor_B*dr.z;
              }
            }
          }
        }
      }
    }
  }

  if(!OmitCationCationPolarization)
  {
    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      if(m!=excl_cation)
      {
        Type=Cations[CurrentSystem][m].Type;
        if(Components[Type].HasCharges)
        {
          energy_sum=0.0;
          NumberOfExcludedPairs=Components[Type].NumberOfExcludedIntraChargeCharge;
          for(i=0;i<NumberOfExcludedPairs;i++)
          {
            A=Components[Type].ExcludedIntraChargeCharge[i].A;
            B=Components[Type].ExcludedIntraChargeCharge[i].B;

            TypeA=Cations[CurrentSystem][m].Atoms[A].Type;
            ChargeA=PseudoAtoms[TypeA].Charge;
            if(PseudoAtoms[TypeA].HasCharges)
            {
              TypeB=Cations[CurrentSystem][m].Atoms[B].Type;
              if(PseudoAtoms[TypeB].HasCharges)
              {
                ChargeB=PseudoAtoms[TypeB].Charge;
                posA=Cations[CurrentSystem][m].Atoms[A].Position;
                posB=Cations[CurrentSystem][m].Atoms[B].Position;

                dr.x=posA.x-posB.x;
                dr.y=posA.y-posB.y;
                dr.z=posA.z-posB.z;
                rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(rr);

                temp=COULOMBIC_CONVERSION_FACTOR*
                     (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                     (r*rr);
                force_factor_A=temp*ChargeB;
                force_factor_B=temp*ChargeA;

                Cations[CurrentSystem][m].Atoms[A].ElectricField.x+=force_factor_A*dr.x;
                Cations[CurrentSystem][m].Atoms[A].ElectricField.y+=force_factor_A*dr.y;
                Cations[CurrentSystem][m].Atoms[A].ElectricField.z+=force_factor_A*dr.z;

                Cations[CurrentSystem][m].Atoms[B].ElectricField.x-=force_factor_B*dr.x;
                Cations[CurrentSystem][m].Atoms[B].ElectricField.y-=force_factor_B*dr.y;
                Cations[CurrentSystem][m].Atoms[B].ElectricField.z-=force_factor_B*dr.z;
              }
            }
          }
        }
      }
    }
  }

  if(New)
  {
    if((Components[CurrentComponent].ExtraFrameworkMolecule&&(!OmitCationCationPolarization))||
      ((!Components[CurrentComponent].ExtraFrameworkMolecule)&&(!OmitAdsorbateAdsorbatePolarization)))
    {
      NumberOfExcludedPairs=Components[CurrentComponent].NumberOfExcludedIntraChargeCharge;
      for(i=0;i<NumberOfExcludedPairs;i++)
      {
        A=Components[CurrentComponent].ExcludedIntraChargeCharge[i].A;
        B=Components[CurrentComponent].ExcludedIntraChargeCharge[i].B;

        TypeA=Components[CurrentSystem].Type[A];
        if(PseudoAtoms[TypeA].HasCharges)
        {
          TypeB=Components[CurrentSystem].Type[B];
          if(PseudoAtoms[TypeB].HasCharges)
          {
            ChargeA=PseudoAtoms[TypeA].Charge;
            ChargeB=PseudoAtoms[TypeB].Charge;
            posA=TrialPosition[CurrentSystem][A];
            posB=TrialPosition[CurrentSystem][B];

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            temp=COULOMBIC_CONVERSION_FACTOR*
                 (2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI)-ErrorFunction(Alpha[CurrentSystem]*r))/
                 (r*rr);
            force_factor_A=temp*ChargeB;
            force_factor_B=temp*ChargeA;

            ElectricFieldAtTrialPosition[CurrentSystem][A].x+=force_factor_A*dr.x;
            ElectricFieldAtTrialPosition[CurrentSystem][A].y+=force_factor_A*dr.y;
            ElectricFieldAtTrialPosition[CurrentSystem][A].z+=force_factor_A*dr.z;

            ElectricFieldAtTrialPosition[CurrentSystem][B].x-=force_factor_B*dr.x;
            ElectricFieldAtTrialPosition[CurrentSystem][B].y-=force_factor_B*dr.y;
            ElectricFieldAtTrialPosition[CurrentSystem][B].z-=force_factor_B*dr.z;
          }
        }
      }
    }
  }
  return 0;
}

int EwaldFourierIntraForceCorrectionChargeInducedDipoleAdsorbate(void)
{
  int i,m,A,B,NumberOfExcludedPairs,Type,TypeA,TypeB;
  REAL ChargeA,ChargeB,energy_sum;
  REAL r,rr;
  REAL cosA,cosB;
  REAL Bt0,Bt1,Bt2,Bt3;
  VECTOR dipoleA,dipoleB,term;
  POINT posA,posB;
  VECTOR dr;

  if(OmitAdsorbateAdsorbateCoulombInteractions)
    return 0;

  if(OmitAdsorbateAdsorbatePolarization)
    return 0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    if(Components[Type].HasCharges)
    {
      energy_sum=0.0;
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
            dipoleA=Adsorbates[CurrentSystem][m].Atoms[A].InducedDipole;
            dipoleB=Adsorbates[CurrentSystem][m].Atoms[B].InducedDipole;
            ChargeA=PseudoAtoms[TypeA].Charge;
            ChargeB=PseudoAtoms[TypeB].Charge;
            posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
            posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(rr);

            switch(ChargeMethod)
            {
              case NONE:
                Bt0=Bt1=Bt2=0.0;
                break;
              case TRUNCATED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case SHIFTED_COULOMB:
                Bt0=1.0/(r);
                Bt1=1.0/(r*rr);
                Bt2=3.0/(r*rr*rr);
                break;
              case EWALD:
              default:
                Bt0=-ErrorFunction(Alpha[CurrentSystem]*r)/r;
                Bt1=2.0*Alpha[CurrentSystem]*exp(-SQR(-Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    -ErrorFunction(Alpha[CurrentSystem]*r)/(rr*r);
                Bt2=6.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                    4.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    -3.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*r);
                Bt3=30.0*Alpha[CurrentSystem]*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr*rr)+
                    20.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr*rr)+
                     8.0*Alpha[CurrentSystem]*SQR(Alpha[CurrentSystem])*SQR(Alpha[CurrentSystem])*exp(-SQR(Alpha[CurrentSystem]*r))/(sqrt(M_PI)*rr)+
                    -15.0*ErrorFunction(Alpha[CurrentSystem]*r)/(rr*rr*rr*r);
                break;
            }

            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;

            term.x=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.x-Bt1*dipoleA.x);
            term.y=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.y-Bt1*dipoleA.y);
            term.z=COULOMBIC_CONVERSION_FACTOR*ChargeB*(Bt2*cosA*dr.z-Bt1*dipoleA.z);

            Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=term.x;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=term.y;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=term.z;

            Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=term.x;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=term.y;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=term.z;

            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

            term.x=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.x-Bt1*dipoleB.x);
            term.y=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.y-Bt1*dipoleB.y);
            term.z=COULOMBIC_CONVERSION_FACTOR*ChargeA*(Bt2*cosB*dr.z-Bt1*dipoleB.z);

            Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=term.x;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=term.y;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=term.z;

            Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=term.x;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=term.y;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=term.z;

/*
            cosAB=dipoleA.x*dipoleB.x+dipoleA.y*dipoleB.y+dipoleA.z*dipoleB.z;
            cosA=dipoleA.x*dr.x+dipoleA.y*dr.y+dipoleA.z*dr.z;
            cosB=dipoleB.x*dr.x+dipoleB.y*dr.y+dipoleB.z*dr.z;

            term.x=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.x-Bt2*(cosAB*dr.x+cosB*dipoleA.x+cosA*dipoleB.x));
            term.y=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.y-Bt2*(cosAB*dr.y+cosB*dipoleA.y+cosA*dipoleB.y));
            term.z=COULOMBIC_CONVERSION_FACTOR*(Bt3*(cosA*cosB)*dr.z-Bt2*(cosAB*dr.z+cosB*dipoleA.z+cosA*dipoleB.z));

            Adsorbates[CurrentSystem][m].Atoms[A].Force.x-=term.x;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.y-=term.y;
            Adsorbates[CurrentSystem][m].Atoms[A].Force.z-=term.z;

            Adsorbates[CurrentSystem][m].Atoms[B].Force.x+=term.x;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.y+=term.y;
            Adsorbates[CurrentSystem][m].Atoms[B].Force.z+=term.z;
*/

          }
        }
      }
    }
  }
  return 0;
}
