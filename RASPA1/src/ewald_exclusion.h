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



#ifndef EWALD_EXCLUSION_H
#define EWALD_EXCLUSION_H

#include "simulation.h"
#include "vector.h"
#include "utils.h"

extern REAL UHostHostChargeChargeExcludedFourier;
extern REAL UHostHostChargeBondDipoleExcludedFourier;
extern REAL UHostHostBondDipoleBondDipoleExcludedFourier;

extern REAL UAdsorbateAdsorbateChargeChargeExcludedFourier;
extern REAL UAdsorbateAdsorbateChargeBondDipoleExcludedFourier;
extern REAL UAdsorbateAdsorbateBondDipoleBondDipoleExcludedFourier;

extern REAL UCationCationChargeChargeExcludedFourier;
extern REAL UCationCationChargeBondDipoleExcludedFourier;
extern REAL UCationCationBondDipoleBondDipoleExcludedFourier;

int EwaldFourierIntraCorrectionAdsorbate(void);
int EwaldFourierIntraCorrectionCation(void);
void EwaldFourierIntraForceCorrectionChargeChargeFramework(void);
void EwaldFourierIntraForceCorrectionChargeBondDipoleFramework(void);
void EwaldFourierIntraForceCorrectionBondDipoleBondDipoleFramework(void);

int EwaldFourierIntraForceCorrectionChargeChargeAdsorbate(void);
int EwaldFourierIntraForceCorrectionChargeChargeCation(void);

int EwaldFourierIntraForceCorrectionChargeBondDipoleAdsorbate(void);
int EwaldFourierIntraForceCorrectionChargeBondDipoleCation(void);

int EwaldFourierIntraForceCorrectionBondDipoleBondDipoleAdsorbate(void);
int EwaldFourierIntraForceCorrectionBondDipoleBondDipoleCation(void);

void EwaldFourierIntraForceCorrectionFrameworkBornTerm(void);
void EwaldFourierIntraForceCorrectionAdsorbateBornTerm(void);
void EwaldFourierIntraForceCorrectionCationBornTerm(void);

void EwaldFourierIntraElectricFieldCorrectionChargeChargeFramework(void);
int EwaldFourierIntraElectricFieldCorrectionChargeChargeAdsorbate(void);
int EwaldFourierIntraElectricFieldCorrectionChargeChargeCation(void);

int EwaldFourierIntraElectricFieldCorrectionChargeChargeMC(int New, int excl_ads,int int_cation);

int EwaldFourierIntraForceCorrectionChargeInducedDipoleAdsorbate(void);

#endif
