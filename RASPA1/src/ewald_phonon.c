/*****************************************************************************************************
    ewald_phonon.c -  description
    -----------------------------
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


// Remark 10 February 2011: Check and fix atoms without charge
int CalculateEwaldFourierPhonon(VECTOR k,REAL *Energy,REAL* Gradient,COMPLEX_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainFirstDerivative,int ComputeGradient,int ComputeHessian)
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
  VECTOR drA,drB;
  COMPLEX phase_factor;
  REAL dot_product;

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
            index1=0;
            if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
            {
              for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
              {
                for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
                {
                  if(PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].HasCharges)
                  {
                    drA=Framework[CurrentSystem].Atoms[f1][i].Position;

                    index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;
                    index_i2=-1;

                    f1_I=2.0*COULOMBIC_CONVERSION_FACTOR*(4.0*M_PI/Volume[CurrentSystem])*exp_term*
                         (-Ck[index1].im*Cksum.re+Ck[index1].re*Cksum.im);
                    f2_I=2.0*COULOMBIC_CONVERSION_FACTOR*exp_term*(4.0*M_PI/Volume[CurrentSystem])*
                         (-Ck[index1].re*Cksum.re-Ck[index1].im*Cksum.im);

                    if(index_i>=0)
                    {
                      HessianMatrix.element[index_i][index_i].re+=f2_I*Rk.x*Rk.x;
                      HessianMatrix.element[index_i][index_i+1].re+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i][index_i+2].re+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+1][index_i+1].re+=f2_I*Rk.y*Rk.y;
                      HessianMatrix.element[index_i+1][index_i+2].re+=f2_I*Rk.y*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+2].re+=f2_I*Rk.z*Rk.z;

                      HessianMatrix.element[index_i+1][index_i].re+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i+2][index_i].re+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+1].re+=f2_I*Rk.y*Rk.z;

                      index2=0;
                      for(f2=0;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
                      {
                        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f2];j++)
                        {
                          if(PseudoAtoms[Framework[CurrentSystem].Atoms[f2][j].Type].HasCharges)
                          {
                            drB=Framework[CurrentSystem].Atoms[f2][j].Position;
                            index_j=Framework[CurrentSystem].Atoms[f2][j].HessianIndex;
                            index_j2=-1;

                            // form e(-k.r)
                            dr.x=drA.x-drB.x;
                            dr.y=drA.y-drB.y;
                            dr.z=drA.z-drB.z;
                            dr=ApplyBoundaryCondition(dr);

                            dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
                            phase_factor.re=cos(dot_product);
                            phase_factor.im=-sin(dot_product);

                            if((index_i>=0)&&(index_j>=0))
                            {
                              if(index_i<=index_j)
                              {
                                f2_IJ=2.0*COULOMBIC_CONVERSION_FACTOR*exp_term*(4.0*M_PI/Volume[CurrentSystem])*
                                      (-Ck[index2].re*Ck[index1].re-Ck[index2].im*Ck[index1].im);
                                HessianMatrix.element[index_i][index_j].re-=phase_factor.re*f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j].im-=phase_factor.im*f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+1].im-=phase_factor.im*f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i][index_j+2].im-=phase_factor.im*f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j].im-=phase_factor.im*f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+1].im-=phase_factor.im*f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+1][index_j+2].im-=phase_factor.im*f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j].im-=phase_factor.im*f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+1].im-=phase_factor.im*f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*Rk.z;
                                HessianMatrix.element[index_i+2][index_j+2].im-=phase_factor.im*f2_IJ*Rk.z*Rk.z;
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

                              if(index_i<=index_j)
                              {
                                if(index_j>=0)
                                {
                                  // TODO
                                  HessianMatrix.element[index_i][index_j].re-=phase_factor.re*f2_IJ*Rk.x*Rk.x;
                                  HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*f2_IJ*Rk.x*Rk.y;
                                  HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*f2_IJ*Rk.x*Rk.z;
                                  HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*f2_IJ*Rk.y*Rk.x;
                                  HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*Rk.y;
                                  HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.y*Rk.z;
                                  HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*f2_IJ*Rk.z*Rk.x;
                                  HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.z*Rk.y;
                                  HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*Rk.z;
                                }

                                // com of I with orientation of J
                                if(Components[TypeMolB].Groups[jg].Rigid)
                                {
                                  if((index_i>=0)&&(index_j2>=0))
                                  {
                                    HessianMatrix.element[index_i][index_j2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.x;
                                    HessianMatrix.element[index_i][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.y;
                                    HessianMatrix.element[index_i][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.z;
                                    HessianMatrix.element[index_i+1][index_j2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.x;
                                    HessianMatrix.element[index_i+1][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.y;
                                    HessianMatrix.element[index_i+1][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.z;
                                    HessianMatrix.element[index_i+2][index_j2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.x;
                                    HessianMatrix.element[index_i+2][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.y;
                                    HessianMatrix.element[index_i+2][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.z;
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

                              if(index_i<=index_j)
                              {
                                if((index_i>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i][index_j].re-=phase_factor.re*f2_IJ*Rk.x*Rk.x;
                                  HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*f2_IJ*Rk.x*Rk.y;
                                  HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*f2_IJ*Rk.x*Rk.z;
                                  HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*f2_IJ*Rk.y*Rk.x;
                                  HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*Rk.y;
                                  HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.y*Rk.z;
                                  HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*f2_IJ*Rk.z*Rk.x;
                                  HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.z*Rk.y;
                                  HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*Rk.z;
                                }

                                // com of I with orientation of J
                                if(Components[TypeMolB].Groups[jg].Rigid)
                                {
                                  if((index_i>=0)&&(index_j2>=0))
                                  {
                                    HessianMatrix.element[index_i][index_j2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.x;
                                    HessianMatrix.element[index_i][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.y;
                                    HessianMatrix.element[index_i][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.z;
                                    HessianMatrix.element[index_i+1][index_j2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.x;
                                    HessianMatrix.element[index_i+1][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.y;
                                    HessianMatrix.element[index_i+1][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.z;
                                    HessianMatrix.element[index_i+2][index_j2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.x;
                                    HessianMatrix.element[index_i+2][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.y;
                                    HessianMatrix.element[index_i+2][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.z;
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
                      HessianMatrix.element[index_i][index_i].re+=f2_I*Rk.x*Rk.x;
                      HessianMatrix.element[index_i][index_i+1].re+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i][index_i+2].re+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+1][index_i+1].re+=f2_I*Rk.y*Rk.y;
                      HessianMatrix.element[index_i+1][index_i+2].re+=f2_I*Rk.y*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+2].re+=f2_I*Rk.z*Rk.z;

                      HessianMatrix.element[index_i+1][index_i].re+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i+2][index_i].re+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+1].re+=f2_I*Rk.y*Rk.z;
                    }


                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      pos=Components[TypeMolA].Positions[i];

                      if(index_i2>=0)
                      {
                        HessianMatrix.element[index_i2][index_i2].re+=f2_I*SQR(dot_product_i.x);
                        HessianMatrix.element[index_i2+1][index_i2+1].re+=f2_I*SQR(dot_product_i.y);
                        HessianMatrix.element[index_i2+2][index_i2+2].re+=f2_I*SQR(dot_product_i.z);
                        HessianMatrix.element[index_i2][index_i2+1].re+=f2_I*dot_product_i.x*dot_product_i.y;
                        HessianMatrix.element[index_i2][index_i2+2].re+=f2_I*dot_product_i.x*dot_product_i.z;
                        HessianMatrix.element[index_i2+1][index_i2+2].re+=f2_I*dot_product_i.y*dot_product_i.z;
                      }

                      if((index_i>=0)&&(index_i2>=0))
                      {
                        HessianMatrix.element[index_i][index_i2].re+=f2_I*Rk.x*dot_product_i.x;
                        HessianMatrix.element[index_i][index_i2+1].re+=f2_I*Rk.x*dot_product_i.y;
                        HessianMatrix.element[index_i][index_i2+2].re+=f2_I*Rk.x*dot_product_i.z;
                        HessianMatrix.element[index_i+1][index_i2].re+=f2_I*Rk.y*dot_product_i.x;
                        HessianMatrix.element[index_i+1][index_i2+1].re+=f2_I*Rk.y*dot_product_i.y;
                        HessianMatrix.element[index_i+1][index_i2+2].re+=f2_I*Rk.y*dot_product_i.z;
                        HessianMatrix.element[index_i+2][index_i2].re+=f2_I*Rk.z*dot_product_i.x;
                        HessianMatrix.element[index_i+2][index_i2+1].re+=f2_I*Rk.z*dot_product_i.y;
                        HessianMatrix.element[index_i+2][index_i2+2].re+=f2_I*Rk.z*dot_product_i.z;
                      }

                      if(index_i2>=0)
                      {
                        HessianMatrix.element[index_i2][index_i2].re+=f1_I*dot_product_AX;
                        HessianMatrix.element[index_i2+1][index_i2+1].re+=f1_I*dot_product_BY;
                        HessianMatrix.element[index_i2+2][index_i2+2].re+=f1_I*dot_product_CZ;
                        HessianMatrix.element[index_i2][index_i2+1].re+=f1_I*dot_product_AY;
                        HessianMatrix.element[index_i2][index_i2+2].re+=f1_I*dot_product_AZ;
                        HessianMatrix.element[index_i2+1][index_i2+2].re+=f1_I*dot_product_BZ;
                      }
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

                            // form e(-k.r)
                            if(Components[TypeMolA].Groups[ig].Rigid) drA=comA; else drA=posA; 
                            if(Components[TypeMolB].Groups[jg].Rigid) drB=comB; else drB=posB; 
                            dr.x=drA.x-drB.x;
                            dr.y=drA.y-drB.y;
                            dr.z=drA.z-drB.z;
                            dr=ApplyBoundaryCondition(dr);

                            dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
                            phase_factor.re=cos(dot_product);
                            phase_factor.im=-sin(dot_product);

                            if(index_i<=index_j)
                            {
                              if((index_i>=0)&&(index_j>=0))
                              {
                                HessianMatrix.element[index_i][index_j].re-=phase_factor.re*f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j].im-=phase_factor.im*f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+1].im-=phase_factor.im*f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i][index_j+2].im-=phase_factor.im*f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j].im-=phase_factor.im*f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+1].im-=phase_factor.im*f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+1][index_j+2].im-=phase_factor.im*f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j].im-=phase_factor.im*f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+1].im-=phase_factor.im*f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*Rk.z;
                                HessianMatrix.element[index_i+2][index_j+2].im-=phase_factor.im*f2_IJ*Rk.z*Rk.z;
                              }

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i][index_j2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.x;
                                  HessianMatrix.element[index_i][index_j2].im-=phase_factor.im*f2_IJ*Rk.x*dot_product_j.x;
                                  HessianMatrix.element[index_i][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.y;
                                  HessianMatrix.element[index_i][index_j2+1].im-=phase_factor.im*f2_IJ*Rk.x*dot_product_j.y;
                                  HessianMatrix.element[index_i][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.z;
                                  HessianMatrix.element[index_i][index_j2+2].im-=phase_factor.im*f2_IJ*Rk.x*dot_product_j.z;
                                  HessianMatrix.element[index_i+1][index_j2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.x;
                                  HessianMatrix.element[index_i+1][index_j2].im-=phase_factor.im*f2_IJ*Rk.y*dot_product_j.x;
                                  HessianMatrix.element[index_i+1][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.y;
                                  HessianMatrix.element[index_i+1][index_j2+1].im-=phase_factor.im*f2_IJ*Rk.y*dot_product_j.y;
                                  HessianMatrix.element[index_i+1][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.z;
                                  HessianMatrix.element[index_i+1][index_j2+2].im-=phase_factor.im*f2_IJ*Rk.y*dot_product_j.z;
                                  HessianMatrix.element[index_i+2][index_j2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.x;
                                  HessianMatrix.element[index_i+2][index_j2].im-=phase_factor.im*f2_IJ*Rk.z*dot_product_j.x;
                                  HessianMatrix.element[index_i+2][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.y;
                                  HessianMatrix.element[index_i+2][index_j2+1].im-=phase_factor.im*f2_IJ*Rk.z*dot_product_j.y;
                                  HessianMatrix.element[index_i+2][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.z;
                                  HessianMatrix.element[index_i+2][index_j2+2].im-=phase_factor.im*f2_IJ*Rk.z*dot_product_j.z;
                                }
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.x;
                                  HessianMatrix.element[index_i2][index_j].im-=phase_factor.im*f2_IJ*Rk.x*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.y;
                                  HessianMatrix.element[index_i2+1][index_j].im-=phase_factor.im*f2_IJ*Rk.x*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.z;
                                  HessianMatrix.element[index_i2+2][index_j].im-=phase_factor.im*f2_IJ*Rk.x*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.x;
                                  HessianMatrix.element[index_i2][index_j+1].im-=phase_factor.im*f2_IJ*Rk.y*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.y;
                                  HessianMatrix.element[index_i2+1][index_j+1].im-=phase_factor.im*f2_IJ*Rk.y*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.z;
                                  HessianMatrix.element[index_i2+2][index_j+1].im-=phase_factor.im*f2_IJ*Rk.y*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.x;
                                  HessianMatrix.element[index_i2][index_j+2].im-=phase_factor.im*f2_IJ*Rk.z*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.y;
                                  HessianMatrix.element[index_i2+1][index_j+2].im-=phase_factor.im*f2_IJ*Rk.z*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.z;
                                  HessianMatrix.element[index_i2+2][index_j+2].im-=phase_factor.im*f2_IJ*Rk.z*dot_product_i.z;
                                }
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.x;
                                  HessianMatrix.element[index_i2][index_j2].im-=phase_factor.im*f2_IJ*dot_product_i.x*dot_product_j.x;
                                  HessianMatrix.element[index_i2][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.y;
                                  HessianMatrix.element[index_i2][index_j2+1].im-=phase_factor.im*f2_IJ*dot_product_i.x*dot_product_j.y;
                                  HessianMatrix.element[index_i2][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.z;
                                  HessianMatrix.element[index_i2][index_j2+2].im-=phase_factor.im*f2_IJ*dot_product_i.x*dot_product_j.z;
                                  HessianMatrix.element[index_i2+1][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.x;
                                  HessianMatrix.element[index_i2+1][index_j2].im-=phase_factor.im*f2_IJ*dot_product_i.y*dot_product_j.x;
                                  HessianMatrix.element[index_i2+1][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.y;
                                  HessianMatrix.element[index_i2+1][index_j2+1].im-=phase_factor.im*f2_IJ*dot_product_i.y*dot_product_j.y;
                                  HessianMatrix.element[index_i2+1][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.z;
                                  HessianMatrix.element[index_i2+1][index_j2+2].im-=phase_factor.im*f2_IJ*dot_product_i.y*dot_product_j.z;
                                  HessianMatrix.element[index_i2+2][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.x;
                                  HessianMatrix.element[index_i2+2][index_j2].im-=phase_factor.im*f2_IJ*dot_product_i.z*dot_product_j.x;
                                  HessianMatrix.element[index_i2+2][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.y;
                                  HessianMatrix.element[index_i2+2][index_j2+1].im-=phase_factor.im*f2_IJ*dot_product_i.z*dot_product_j.y;
                                  HessianMatrix.element[index_i2+2][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.z;
                                  HessianMatrix.element[index_i2+2][index_j2+2].im-=phase_factor.im*f2_IJ*dot_product_i.z*dot_product_j.z;
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

                            if(index_i<=index_j)
                            {
                              if((index_i>=0)&&(index_j>=0))
                              {
                                HessianMatrix.element[index_i][index_j].re-=phase_factor.re*f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*Rk.z;
                              }

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i][index_j2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.x;
                                  HessianMatrix.element[index_i][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.y;
                                  HessianMatrix.element[index_i][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.z;
                                  HessianMatrix.element[index_i+1][index_j2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.x;
                                  HessianMatrix.element[index_i+1][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.y;
                                  HessianMatrix.element[index_i+1][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.z;
                                  HessianMatrix.element[index_i+2][index_j2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.x;
                                  HessianMatrix.element[index_i+2][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.y;
                                  HessianMatrix.element[index_i+2][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.z;
                                }
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.z;
                                }
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.x;
                                  HessianMatrix.element[index_i2][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.y;
                                  HessianMatrix.element[index_i2][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.z;
                                  HessianMatrix.element[index_i2+1][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.x;
                                  HessianMatrix.element[index_i2+1][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.y;
                                  HessianMatrix.element[index_i2+1][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.z;
                                  HessianMatrix.element[index_i2+2][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.x;
                                  HessianMatrix.element[index_i2+2][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.y;
                                  HessianMatrix.element[index_i2+2][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.z;
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
                      HessianMatrix.element[index_i][index_i].re+=f2_I*Rk.x*Rk.x;
                      HessianMatrix.element[index_i][index_i+1].re+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i][index_i+2].re+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+1][index_i+1].re+=f2_I*Rk.y*Rk.y;
                      HessianMatrix.element[index_i+1][index_i+2].re+=f2_I*Rk.y*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+2].re+=f2_I*Rk.z*Rk.z;

                      HessianMatrix.element[index_i+1][index_i].re+=f2_I*Rk.x*Rk.y;
                      HessianMatrix.element[index_i+2][index_i].re+=f2_I*Rk.x*Rk.z;
                      HessianMatrix.element[index_i+2][index_i+1].re+=f2_I*Rk.y*Rk.z;
                    }

                    if(Components[TypeMolA].Groups[ig].Rigid)
                    {
                      pos=Components[TypeMolA].Positions[i];

                      if(index_i2>=0)
                      {
                        HessianMatrix.element[index_i2][index_i2].re+=f2_I*dot_product_i.x*dot_product_i.x;
                        HessianMatrix.element[index_i2+1][index_i2+1].re+=f2_I*dot_product_i.y*dot_product_i.y;
                        HessianMatrix.element[index_i2+2][index_i2+2].re+=f2_I*dot_product_i.z*dot_product_i.z;
                        HessianMatrix.element[index_i2][index_i2+1].re+=f2_I*dot_product_i.x*dot_product_i.y;
                        HessianMatrix.element[index_i2][index_i2+2].re+=f2_I*dot_product_i.x*dot_product_i.z;
                        HessianMatrix.element[index_i2+1][index_i2+2].re+=f2_I*dot_product_i.y*dot_product_i.z;
                      }

                      if((index_i>=0)&&(index_i2>=0))
                      { 
                        HessianMatrix.element[index_i][index_i2].re+=f2_I*Rk.x*dot_product_i.x;
                        HessianMatrix.element[index_i][index_i2+1].re+=f2_I*Rk.x*dot_product_i.y;
                        HessianMatrix.element[index_i][index_i2+2].re+=f2_I*Rk.x*dot_product_i.z;
                        HessianMatrix.element[index_i+1][index_i2].re+=f2_I*Rk.y*dot_product_i.x;
                        HessianMatrix.element[index_i+1][index_i2+1].re+=f2_I*Rk.y*dot_product_i.y;
                        HessianMatrix.element[index_i+1][index_i2+2].re+=f2_I*Rk.y*dot_product_i.z;
                        HessianMatrix.element[index_i+2][index_i2].re+=f2_I*Rk.z*dot_product_i.x;
                        HessianMatrix.element[index_i+2][index_i2+1].re+=f2_I*Rk.z*dot_product_i.y;
                        HessianMatrix.element[index_i+2][index_i2+2].re+=f2_I*Rk.z*dot_product_i.z;
                      }

                      if(index_i2>=0)
                      {
                        pos=Components[TypeMolA].Positions[i];
                        HessianMatrix.element[index_i2][index_i2].re+=f1_I*dot_product_AX;
                        HessianMatrix.element[index_i2+1][index_i2+1].re+=f1_I*dot_product_BY;
                        HessianMatrix.element[index_i2+2][index_i2+2].re+=f1_I*dot_product_CZ;
                        HessianMatrix.element[index_i2][index_i2+1].re+=f1_I*dot_product_AY;
                        HessianMatrix.element[index_i2][index_i2+2].re+=f1_I*dot_product_AZ;
                        HessianMatrix.element[index_i2+1][index_i2+2].re+=f1_I*dot_product_BZ;
                      }
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

                            if(index_i<=index_j)
                            {
                              if((index_i>=0)&(index_j>=0))
                              { 
                                HessianMatrix.element[index_i][index_j].re-=phase_factor.re*f2_IJ*Rk.x*Rk.x;
                                HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*f2_IJ*Rk.x*Rk.y;
                                HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*f2_IJ*Rk.x*Rk.z;
                                HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*f2_IJ*Rk.y*Rk.x;
                                HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*Rk.y;
                                HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.y*Rk.z;
                                HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*f2_IJ*Rk.z*Rk.x;
                                HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.z*Rk.y;
                                HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*Rk.z;
                              }

                              // com of I with orientation of J
                              if(Components[TypeMolB].Groups[jg].Rigid)
                              {
                                if((index_i>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i][index_j2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.x;
                                  HessianMatrix.element[index_i][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.y;
                                  HessianMatrix.element[index_i][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_j.z;
                                  HessianMatrix.element[index_i+1][index_j2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.x;
                                  HessianMatrix.element[index_i+1][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.y;
                                  HessianMatrix.element[index_i+1][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_j.z;
                                  HessianMatrix.element[index_i+2][index_j2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.x;
                                  HessianMatrix.element[index_i+2][index_j2+1].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.y;
                                  HessianMatrix.element[index_i+2][index_j2+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_j.z;
                                }
                              }

                              // orientation of I with com of J
                              if(Components[TypeMolA].Groups[ig].Rigid)
                              {
                                if((index_i2>=0)&&(index_j>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j].re-=phase_factor.re*f2_IJ*Rk.x*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+1].re-=phase_factor.re*f2_IJ*Rk.y*dot_product_i.z;
                                  HessianMatrix.element[index_i2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.x;
                                  HessianMatrix.element[index_i2+1][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.y;
                                  HessianMatrix.element[index_i2+2][index_j+2].re-=phase_factor.re*f2_IJ*Rk.z*dot_product_i.z;
                                }
                              }

                              // orientation of I with orientation of J
                              if((Components[TypeMolA].Groups[ig].Rigid)&&(Components[TypeMolB].Groups[jg].Rigid))
                              {
                                if((index_i2>=0)&&(index_j2>=0))
                                {
                                  HessianMatrix.element[index_i2][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.x;
                                  HessianMatrix.element[index_i2][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.y;
                                  HessianMatrix.element[index_i2][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.x*dot_product_j.z;
                                  HessianMatrix.element[index_i2+1][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.x;
                                  HessianMatrix.element[index_i2+1][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.y;
                                  HessianMatrix.element[index_i2+1][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.y*dot_product_j.z;
                                  HessianMatrix.element[index_i2+2][index_j2].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.x;
                                  HessianMatrix.element[index_i2+2][index_j2+1].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.y;
                                  HessianMatrix.element[index_i2+2][index_j2+2].re-=phase_factor.re*f2_IJ*dot_product_i.z*dot_product_j.z;
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

        // form e(-k.r)
        dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
        phase_factor.re=cos(dot_product);
        phase_factor.im=-sin(dot_product);

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
            HessianMatrix.element[index_i][index_i].re+=Hessian.ax;
            HessianMatrix.element[index_i][index_i+1].re+=Hessian.ay;
            HessianMatrix.element[index_i][index_i+2].re+=Hessian.az;
            HessianMatrix.element[index_i+1][index_i+1].re+=Hessian.by;
            HessianMatrix.element[index_i+1][index_i+2].re+=Hessian.bz;
            HessianMatrix.element[index_i+2][index_i+2].re+=Hessian.cz;
          }

          if(index_j>=0)
          {
            HessianMatrix.element[index_j][index_j].re+=Hessian.ax;
            HessianMatrix.element[index_j][index_j+1].re+=Hessian.ay;
            HessianMatrix.element[index_j][index_j+2].re+=Hessian.az;
            HessianMatrix.element[index_j+1][index_j+1].re+=Hessian.by;
            HessianMatrix.element[index_j+1][index_j+2].re+=Hessian.bz;
            HessianMatrix.element[index_j+2][index_j+2].re+=Hessian.cz;
          }

          if((index_i>=0)&&(index_j>=0))
          {
            if(index_i<index_j)
            {
              HessianMatrix.element[index_i][index_j].re-=phase_factor.re*Hessian.ax;
              HessianMatrix.element[index_i][index_j].im-=phase_factor.im*Hessian.ax;
              HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*Hessian.ay;
              HessianMatrix.element[index_i][index_j+1].im-=phase_factor.im*Hessian.ay;
              HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*Hessian.az;
              HessianMatrix.element[index_i][index_j+2].im-=phase_factor.im*Hessian.az;
              HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*Hessian.ay;
              HessianMatrix.element[index_i+1][index_j].im-=phase_factor.im*Hessian.ay;
              HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*Hessian.by;
              HessianMatrix.element[index_i+1][index_j+1].im-=phase_factor.im*Hessian.by;
              HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*Hessian.bz;
              HessianMatrix.element[index_i+1][index_j+2].im-=phase_factor.im*Hessian.bz;
              HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*Hessian.az;
              HessianMatrix.element[index_i+2][index_j].im-=phase_factor.im*Hessian.az;
              HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*Hessian.bz;
              HessianMatrix.element[index_i+2][index_j+1].im-=phase_factor.im*Hessian.bz;
              HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*Hessian.cz;
              HessianMatrix.element[index_i+2][index_j+2].im-=phase_factor.im*Hessian.cz;
            }
            else
            {
              HessianMatrix.element[index_j][index_i].re-=phase_factor.re*Hessian.ax;
              HessianMatrix.element[index_j][index_i].im-=phase_factor.im*Hessian.ax;
              HessianMatrix.element[index_j][index_i+1].re-=phase_factor.re*Hessian.ay;
              HessianMatrix.element[index_j][index_i+1].im-=phase_factor.im*Hessian.ay;
              HessianMatrix.element[index_j][index_i+2].re-=phase_factor.re*Hessian.az;
              HessianMatrix.element[index_j][index_i+2].im-=phase_factor.im*Hessian.az;
              HessianMatrix.element[index_j+1][index_i].re-=phase_factor.re*Hessian.ay;
              HessianMatrix.element[index_j+1][index_i].im-=phase_factor.im*Hessian.ay;
              HessianMatrix.element[index_j+1][index_i+1].re-=phase_factor.re*Hessian.by;
              HessianMatrix.element[index_j+1][index_i+1].im-=phase_factor.im*Hessian.by;
              HessianMatrix.element[index_j+1][index_i+2].re-=phase_factor.re*Hessian.bz;
              HessianMatrix.element[index_j+1][index_i+2].im-=phase_factor.im*Hessian.bz;
              HessianMatrix.element[index_j+2][index_i].re-=phase_factor.re*Hessian.az;
              HessianMatrix.element[index_j+2][index_i].im-=phase_factor.im*Hessian.az;
              HessianMatrix.element[index_j+2][index_i+1].re-=phase_factor.re*Hessian.bz;
              HessianMatrix.element[index_j+2][index_i+1].im-=phase_factor.im*Hessian.bz;
              HessianMatrix.element[index_j+2][index_i+2].re-=phase_factor.re*Hessian.cz;
              HessianMatrix.element[index_j+2][index_i+2].im-=phase_factor.im*Hessian.cz;
            }
          }
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

            // form e(-k.r)
            dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
            phase_factor.re=cos(dot_product);
            phase_factor.im=-sin(dot_product);

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
                HessianMatrix.element[index_i][index_i].re+=Hessian.ax;
                HessianMatrix.element[index_i][index_i+1].re+=Hessian.ay;
                HessianMatrix.element[index_i][index_i+2].re+=Hessian.az;
                HessianMatrix.element[index_i+1][index_i+1].re+=Hessian.by;
                HessianMatrix.element[index_i+1][index_i+2].re+=Hessian.bz;
                HessianMatrix.element[index_i+2][index_i+2].re+=Hessian.cz;
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][index_j].re+=Hessian.ax;
                HessianMatrix.element[index_j][index_j+1].re+=Hessian.ay;
                HessianMatrix.element[index_j][index_j+2].re+=Hessian.az;
                HessianMatrix.element[index_j+1][index_j+1].re+=Hessian.by;
                HessianMatrix.element[index_j+1][index_j+2].re+=Hessian.bz;
                HessianMatrix.element[index_j+2][index_j+2].re+=Hessian.cz;
              }

              if((index_i>=0)&&(index_j>=0))
              {
                if(index_i<index_j)
                {
                  HessianMatrix.element[index_i][index_j].re-=phase_factor.re*Hessian.ax;
                  HessianMatrix.element[index_i][index_j].im-=phase_factor.im*Hessian.ax;
                  HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_i][index_j+1].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_i][index_j+2].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_i+1][index_j].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*Hessian.by;
                  HessianMatrix.element[index_i+1][index_j+1].im-=phase_factor.im*Hessian.by;
                  HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_i+1][index_j+2].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_i+2][index_j].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j+1].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*Hessian.cz;
                  HessianMatrix.element[index_i+2][index_j+2].im-=phase_factor.im*Hessian.cz;
                }
                else
                {
                  HessianMatrix.element[index_j][index_i].re-=phase_factor.re*Hessian.ax;
                  HessianMatrix.element[index_j][index_i].im-=phase_factor.im*Hessian.ax;
                  HessianMatrix.element[index_j][index_i+1].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_j][index_i+1].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_j][index_i+2].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_j][index_i+2].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_j+1][index_i].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_j+1][index_i].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_j+1][index_i+1].re-=phase_factor.re*Hessian.by;
                  HessianMatrix.element[index_j+1][index_i+1].im-=phase_factor.im*Hessian.by;
                  HessianMatrix.element[index_j+1][index_i+2].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_j+1][index_i+2].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_j+2][index_i].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_j+2][index_i+1].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i+1].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i+2].re-=phase_factor.re*Hessian.cz;
                  HessianMatrix.element[index_j+2][index_i+2].im-=phase_factor.im*Hessian.cz;
                }
              }
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

            // form e(-k.r)
            dot_product=k.x*dr.x+k.y*dr.y+k.z*dr.z;
            phase_factor.re=cos(dot_product);
            phase_factor.im=-sin(dot_product);

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
                HessianMatrix.element[index_i][index_i].re+=Hessian.ax;
                HessianMatrix.element[index_i][index_i+1].re+=Hessian.ay;
                HessianMatrix.element[index_i][index_i+2].re+=Hessian.az;
                HessianMatrix.element[index_i+1][index_i+1].re+=Hessian.by;
                HessianMatrix.element[index_i+1][index_i+2].re+=Hessian.bz;
                HessianMatrix.element[index_i+2][index_i+2].re+=Hessian.cz;
              }

              if(index_j>=0)
              {
                HessianMatrix.element[index_j][index_j].re+=Hessian.ax;
                HessianMatrix.element[index_j][index_j+1].re+=Hessian.ay;
                HessianMatrix.element[index_j][index_j+2].re+=Hessian.az;
                HessianMatrix.element[index_j+1][index_j+1].re+=Hessian.by;
                HessianMatrix.element[index_j+1][index_j+2].re+=Hessian.bz;
                HessianMatrix.element[index_j+2][index_j+2].re+=Hessian.cz;
              }

              if((index_i>=0)&&(index_j>=0))
              {
                if(index_i<index_j)
                {
                  HessianMatrix.element[index_i][index_j].re-=phase_factor.re*Hessian.ax;
                  HessianMatrix.element[index_i][index_j].im-=phase_factor.im*Hessian.ax;
                  HessianMatrix.element[index_i][index_j+1].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_i][index_j+1].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_i][index_j+2].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_i][index_j+2].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_i+1][index_j].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_i+1][index_j].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_i+1][index_j+1].re-=phase_factor.re*Hessian.by;
                  HessianMatrix.element[index_i+1][index_j+1].im-=phase_factor.im*Hessian.by;
                  HessianMatrix.element[index_i+1][index_j+2].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_i+1][index_j+2].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_i+2][index_j].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_i+2][index_j+1].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j+1].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_i+2][index_j+2].re-=phase_factor.re*Hessian.cz;
                  HessianMatrix.element[index_i+2][index_j+2].im-=phase_factor.im*Hessian.cz;
                }
                else
                {
                  HessianMatrix.element[index_j][index_i].re-=phase_factor.re*Hessian.ax;
                  HessianMatrix.element[index_j][index_i].im-=phase_factor.im*Hessian.ax;
                  HessianMatrix.element[index_j][index_i+1].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_j][index_i+1].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_j][index_i+2].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_j][index_i+2].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_j+1][index_i].re-=phase_factor.re*Hessian.ay;
                  HessianMatrix.element[index_j+1][index_i].im-=phase_factor.im*Hessian.ay;
                  HessianMatrix.element[index_j+1][index_i+1].re-=phase_factor.re*Hessian.by;
                  HessianMatrix.element[index_j+1][index_i+1].im-=phase_factor.im*Hessian.by;
                  HessianMatrix.element[index_j+1][index_i+2].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_j+1][index_i+2].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i].re-=phase_factor.re*Hessian.az;
                  HessianMatrix.element[index_j+2][index_i].im-=phase_factor.im*Hessian.az;
                  HessianMatrix.element[index_j+2][index_i+1].re-=phase_factor.re*Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i+1].im-=phase_factor.im*Hessian.bz;
                  HessianMatrix.element[index_j+2][index_i+2].re-=phase_factor.re*Hessian.cz;
                  HessianMatrix.element[index_j+2][index_i+2].im-=phase_factor.im*Hessian.cz;
                }
              }
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
