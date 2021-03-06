/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'input.c' is part of RASPA.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include "molecule.h"
#include "framework.h"
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "input.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "ewald.h"
#include "grids.h"
#include "spacegroup.h"
#include "recrossing.h"
#include "movies.h"
#include "utils.h"
#include "sample.h"
#include "spectra.h"
#include "rigid.h"
#include "integration.h"
#include "statistics.h"
#include "thermo_baro_stats.h"
#include "equations_of_state.h"
#include "inter_energy.h"
#include "minimization.h"
#include "scattering_factors.h"
#include "framework_force.h"
#include "warnings.h"
#include "internal_force.h"
#include "status.h"
#include "charge_equilibration.h"

int EwaldAutomatic;
static int *read_frameworks;
static int *InitializeBox;

int IsActiveAtomType(int type)
{
  int i;

  for(i=0;i<NumberOfActiveAtomTypes;i++)
   if(ReturnPseudoAtomNumber(ActiveAtomTypes[i])==type) return TRUE;
  return FALSE;
}

int IsFixedAtomType(int type)
{
  int i;

  for(i=0;i<NumberOfFixedAtomTypes;i++)
   if(ReturnPseudoAtomNumber(FixedAtomTypes[i])==type) return TRUE;
  return FALSE;
}

void CheckConstraintInputFramework(FRAMEWORK_COMPONENT* framework,int framework_nr,int atom_nr)
{
  if(framework[CurrentSystem].NumberOfFrameworks==0)
  {
    printf("Error: framework constraint definition, framework does not exists\n");
    exit(0);
  }
  if(framework_nr>=Framework[CurrentSystem].NumberOfFrameworks)
  {
    printf("Error: framework constraint definition, framework id %d does not exists\n",framework_nr);
    exit(0);
  }
  if(atom_nr>=Framework[CurrentSystem].NumberOfAtoms[framework_nr])
  {
    printf("Error: framework constraint definition, framework-atom %d of framework %d does not exists\n",atom_nr,framework_nr);
    exit(0);
  } 
}

void CheckConstraintInputAdsorbate(ADSORBATE_MOLECULE** adsorbates,int molecule_nr,int atom_nr)
{
  if(!adsorbates[CurrentSystem])
  {
    printf("Error: adsorbate constraint definition, adsorbate-list does not exists\n");
    exit(0);
  }
  if(molecule_nr>=NumberOfAdsorbateMolecules[CurrentSystem])
  {
    printf("Error: adsorbate constraint definition, adsorbate-molecule id %d does not exists\n",molecule_nr);
    exit(0);
  }
  if(atom_nr>=adsorbates[CurrentSystem][molecule_nr].NumberOfAtoms)
  {
    printf("Error: adsorbate constraint definition, adsorbate-atom %d of molecule %d does not exists\n",atom_nr,molecule_nr);
    exit(0);
  }
}

void CheckConstraintInputCation(CATION_MOLECULE** cations,int molecule_nr,int atom_nr)
{
  if(!cations[CurrentSystem])
  {
    printf("Error: cation constraint definition, cation-list does not exists\n");
    exit(0);
  }
  if(molecule_nr>=NumberOfCationMolecules[CurrentSystem])
  {
    printf("Error: cation constraint definition, cation-molecule id %d does not exists\n",molecule_nr);
    exit(0);
  }
  if(atom_nr>=Cations[CurrentSystem][molecule_nr].NumberOfAtoms)
  {
    printf("Error: cations constraint definition, cation-atom %d of molecule %d does not exists\n",atom_nr,molecule_nr);
    exit(0);
  }
}



// Reads at most 'length-1' characters from a file into a buffer
// The character sequence is zero-terminated
// more characters then the line are discarded
char *ReadLine(char *buffer, size_t length, FILE *file)
{
  char *p;
  size_t last;

  if((p=fgets(buffer,length,file)))
  {
    last=strlen(buffer)-1;

    if(buffer[last]=='\n') 
      buffer[last]='\0'; // discard the trailing newline
    else 
    {
      fscanf(file,"%*[^\n]");
      (void) fgetc(file); // discard the newline
    }
  }
  return p;
}

// read the input file and parses the options
// all string are converted to lowercase
// leading spaces are removed
// line-input after a parsed command is ignored
int ReadInputFile(char *inputfilename)
{
  int i,j,k,l,m,n,index,nr_sites;
  int Restart,RestartStyle,Type;
  int StartingBead,Swapable,A1,B1;
  FILE *FilePtr;
  char line[16384],string[1024];
  char keyword[1024],arguments[16384],firstargument[1024],*arg_pointer;
  REAL det,r;
  REAL A,B,C,tempd;
  int temp_int;
  double Wavelength;
  int f1,f2;
  VECTOR posf1,posf2,dr,s;
  int intinput1,intinput2,intinput3,intinput4,intinput5,intinput6,intinput7,intinput8,intinput9,intinput10,intinput11,intinput12;
  char charinput1,charinput2,charinput3,charinput4,charinput5,charinput6;
  REAL realinput1,realinput2;
  int LineNumber;
  int CurrentPrism,CurrentCylinder,CurrentSphere;

  NumberOfFixedAtomTypes=0;
  NumberOfActiveAtomTypes=0;
  FixedAtomTypes=NULL;
  ActiveAtomTypes=NULL;

  CurrentPrism=0;
  CurrentCylinder=0;
  CurrentSphere=0;

  //NumberOfFixedFrameworkAtoms=0;
  //FixedFrameworkAtoms=NULL;

  // not set, means get the random seed from the system time
  seed=0lu;

  UseReducedUnits=FALSE;
  Dimension=3;

  EBCBMC=FALSE;
  UseEBCBMCEnergyFirstBead=TRUE;
  HermiteInterpolation=FALSE;
  NumberOfEBCBMCBins=10000;
  UseTabularGrid=FALSE;
  UseDynamicGrid=FALSE;
  UseDelaunayGrid=FALSE;
  Swapable=FALSE;
  Type=0;
  MaxNumberOfBeads=0;
  MaxNumberOfBonds=0;
  MaxNumberOfUreyBradleys=0;
  MaxNumberOfTorsions=0;
  MaxNumberOfImproperTorsions=0;
  MaxNumberOfOutOfPlanes=0;
  MaxNumberOfIntraVDW=0;
  MaxNumberOfIntraChargeCharge=0;
  MaxNumberOfIntraChargeBondDipole=0;
  MaxNumberOfIntraBondDipoleBondDipole=0;

  CurrentIsothermPressure=0;
  NumberOfIsothermPressures=0;

  ComputeHessian=HESSIAN_ANALYTICAL;
  ComputeLambda=LAMBDA_METHOD_2;
  MinimizationVariables=CARTESIAN;
  UseSymmetryInMinimization=FALSE;
  RemoveTranslationFromHessian=FALSE;
  RemoveRotationFromHessian=FALSE;
  MaximumStepLength=0.3;
  MaximumStepLengthInput=0.3;
  MinimizationConvergenceFactor=1.0;

  NumberOfElementsDistanceHistogram=200;
  MaxRangeDistanceHistogram=12.0;

  WriteVTKGrids=FALSE;

  MovieScale=1.0;

  WritePseudoAtomsToOutput=TRUE;
  WriteForcefieldToOutput=TRUE;
  WriteMoleculeDefinitionToOutput=TRUE;

  Restart=FALSE;
  RestartStyle=RASPA_STYLE;
  Wavelength=0.0;

  PrintPropertiesEvery=10000;
  PrintEvery=5000;
  OptimizeAcceptenceEvery=5000;

  // default values for the restart-file
  ContinueAfterCrash=FALSE;
  WriteBinaryRestartFileEvery=0;

  // default values for equations-of-state
  MultiComponentMixingRules=VAN_DER_WAALS_MIXING_RULES;
  EquationOfState=PENG_ROBINSON;

  // default values for CBMC
  BiasingMethod=LJ_BIASING;
  MinimumInnerCycles=20;
  NumberOfTrialPositions=10;
  NumberOfTrialPositionsForTheFirstBead=10;
  NumberOfTrialPositionsTorsion=100;
  NumberOfTrialMovesPerOpenBead=150;
  EnergyOverlapCriteria=1.0e5;
  MinimumRosenbluthFactor=1.0e-150;
  TargetAccRatioSmallMCScheme=0.4;
  TargetAccRatioTranslation=0.5;

  BlockGridPockets=TRUE;
  BlockGridPores=FALSE;
  BlockEnergyGridOverlapCriteria=EnergyOverlapCriteria/5.0;

  // default value for the cut-off
  CutOffVDW=12.0;
  InverseCutOffVDW=1.0/CutOffVDW;
  CutOffVDWSquared=SQR(CutOffVDW);

  CutOffChargeCharge=12.0;
  CutOffChargeChargeSquared=SQR(CutOffChargeCharge);

  CutOffChargeBondDipole=12.0;
  CutOffChargeBondDipoleSquared=SQR(CutOffChargeBondDipole);

  CutOffBondDipoleBondDipole=12.0;
  CutOffBondDipoleBondDipoleSquared=SQR(CutOffBondDipoleBondDipole);

  CutOffIons=2.0;

  CutOffVDWSwitch=0.0;
  CutOffChargeChargeSwitch=0.0;
  CutOffChargeBondDipoleSwitch=0.0;
  CutOffBondDipoleBondDipoleSwitch=0.0;

  // tail-corrections do not apply on shifted potentials
  ShiftPotentials=TRUE;
  TailCorrections=FALSE;

  // default values for the Ewald-summation
  ChargeMethod=EWALD;
  SymmetrizeFrameworkCharges=TRUE;
  ChargeFromChargeEquilibration=FALSE;
  ChargeEquilibrationMethod=WILMER_SNURR;
  EwaldPrecision=1e-6;
  EwaldAutomatic=TRUE;
  DielectricConstantOfTheMedium=1.0;

  OmitAdsorbateAdsorbateVDWInteractions=FALSE;
  OmitAdsorbateAdsorbateCoulombInteractions=FALSE;
  OmitCationCationVDWInteractions=FALSE;
  OmitCationCationCoulombInteractions=FALSE;
  OmitEwaldFourier=FALSE;

  ComputePolarization=FALSE;
  BackPolarization=FALSE;
  NumberOfBackPolarizationSteps=7;
  OmitIntraFrameworkPolarization=TRUE;
  OmitAdsorbateAdsorbatePolarization=TRUE;
  OmitCationCationPolarization=FALSE;


  SampleEveryInfraRed=10;

  NumberOfBlockElementsMSDOrderN=25;
  MaxNumberOfBlocksMSDOrderN=25;
  ComputeIndividualMSDOrderN=FALSE;
  ComputeMSDOrderNPerPseudoAtom=FALSE;

  NumberOfBlockElementsVACFOrderN=10;
  MaxNumberOfBlocksVACFOrderN=5000;
  ComputeIndividualVACFOrderN=FALSE;
  ComputeVACFOrderNPerPseudoAtom=FALSE;

  NumberOfBlockElementsRVACFOrderN=10;
  MaxNumberOfBlocksRVACFOrderN=5000;
  ComputeIndividualRVACFOrderN=FALSE;
  ComputeRVACFOrderNPerPseudoAtom=FALSE;

  NumberOfBlockElementsMolecularOrientationOrderN=10;
  MaxNumberOfBlocksMolecularOrientationOrderN=5000;

  NumberOfBuffersMSD=20;
  BufferLengthMSD=5000;

  NumberOfBuffersVACF=20;
  BufferLengthVACF=5000;

  //
  ParallelMolFractionComponentA=0;
  ParallelMolFractionComponentB=1;

	// default values for the spacing in Angstrom of grid points
	SpacingVDWGrid=0.15;
	SpacingCoulombGrid=0.15;
	SpacingDelaunayGrid=0.15;
	SpacingDynamicGrid=0.15;

  AsymmetricIons=FALSE;
  CorrectNetChargeOnPseudoAtom=-1;

  DensityProfile3DVTKGridPoints.x=150;
  DensityProfile3DVTKGridPoints.y=150;
  DensityProfile3DVTKGridPoints.z=150;

  WidomParticleInsertionComponent=0;

  ReinitializeVelocities=TRUE;
  QuenchCoreShellVelocities=FALSE;

  ComputePrincipleMomentsOfInertia=FALSE;

  TargetAccRatioVolumeChange=0.5;
  TargetAccRatioBoxShapeChange=0.5;
  TargetAccRatioGibbsVolumeChange=0.5;

  NumberOfHybridNVESteps=5;
  NumberOfHybridNPHSteps=5;
  NumberOfHybridNPHPRSteps=5;

  therm_baro_stats.ThermostatChainLength=3;
  therm_baro_stats.BarostatChainLength=3;

  therm_baro_stats.NumberOfYoshidaSuzukiSteps=5;
  therm_baro_stats.NumberOfRespaSteps=5;

  therm_baro_stats.time_scale_parameter_thermostat=0.15;
  therm_baro_stats.time_scale_parameter_barostat=0.15;

  RemoveBondNeighboursFromLongRangeInteraction=TRUE;
  RemoveBendNeighboursFromLongRangeInteraction=TRUE;
  RemoveTorsionNeighboursFromLongRangeInteraction=FALSE;

  Remove12NeighboursFromChargeChargeInteraction=FALSE;
  Remove13NeighboursFromChargeChargeInteraction=FALSE;
  Remove14NeighboursFromChargeChargeInteraction=FALSE;

  Remove11NeighboursFromChargeBondDipoleInteraction=FALSE;
  Remove12NeighboursFromChargeBondDipoleInteraction=FALSE;
  Remove13NeighboursFromChargeBondDipoleInteraction=FALSE;
  Remove14NeighboursFromChargeBondDipoleInteraction=FALSE;

  Remove12NeighboursFromBondDipoleBondDipoleInteraction=FALSE;
  Remove13NeighboursFromBondDipoleBondDipoleInteraction=FALSE;
  Remove14NeighboursFromBondDipoleBondDipoleInteraction=FALSE;

  InternalFrameworkLennardJonesInteractions=TRUE;

  ComputePowderDiffractionPattern=FALSE;

  ComputeNormalModes=FALSE;
  CorrectNormalModesForConstraints=TRUE;
  MinimumMode=0;
  MaximumMode=20;
  ModeResolution=100;

  MinimizationMethod=BAKER_MINIMIZATION;
  MinimizationType=SEARCH_LOCAL_MINIMUM;
  MinimizationPotentialMethod=ANALYTICALLY;

  MaximumNumberOfMinimizationSteps=10000;
  UseGradientInLineMinimization=TRUE;
  RMSGradientTolerance=1e-6;
  MaxGradientTolerance=1e-6;

  ProbabilityFrameworkChangeMove=0.0;
  ProbabilityVolumeChangeMove=0.0;
  ProbabilityBoxShapeChangeMove=0.0;
  ProbabilityParallelTemperingMove=0.0;
  ProbabilityHyperParallelTemperingMove=0.0;
  ProbabilityParallelMolFractionMove=0.0;
  ProbabilityGibbsVolumeChangeMove=0.0;
  ProbabilityHybridNVEMove=0.0;
  ProbabilityHybridNPHMove=0.0;
  ProbabilityHybridNPHPRMove=0.0;
  ProbabilityFrameworkChangeMove=0.0;
  ProbabilityFrameworkShiftMove=0.0;

  // default time step is 0.5 fs
  DeltaT=0.0005;

  // powder diffraction defaults
  Diffraction.Type=XRAY_DIFFRACTION;
  Diffraction.RadiationType=COPPER_RADIATION;
  Diffraction.lambda_type=DIFFRACTION_SINGLE;
  Diffraction.u=0.01;
  Diffraction.v=-0.001;
  Diffraction.w=0.002;
  Diffraction.asym=0.5;

  Diffraction.two_theta_min=3.0;
  Diffraction.two_theta_max=50.0;
  Diffraction.two_theta_step=0.01;
  Diffraction.PeakShape=DIFFRACTION_PSEUDO_VOIGT;

  TargetAccRatioTranslation=0.5;

  ComputeRattleSteps=FALSE;
  DistanceConstraintType=DISTANCE_R_SQUARED;
  BendConstraintType=THETA;
  DihedralConstraintType=PHI;
  ImproperDihedralConstraintType=PHI;
  InversionBendConstraintType=CHI;
  OutOfPlaneConstraintType=DISTANCE_R;

  PrintFrameworkBondStatus=TRUE;
  PrintFrameworkUreyBradleyStatus=TRUE;
  PrintFrameworkBendStatus=TRUE;
  PrintFrameworkInversionBendStatus=TRUE;
  PrintFrameworkTorsionStatus=TRUE;
  PrintFrameworkImproperTorsionStatus=TRUE;
  PrintFrameworkBondBondStatus=TRUE;
  PrintFrameworkBondBendStatus=TRUE;
  PrintFrameworkBendBendStatus=TRUE;
  PrintFrameworkBendTorsionStatus=TRUE;
  PrintFrameworkIntraVDWStatus=TRUE;
  PrintFrameworkIntraChargeChargeStatus=TRUE;
  PrintFrameworkIntraChargeBondDipoleStatus=TRUE;
  PrintFrameworkIntraBondDipoleBondDipoleStatus=TRUE;

  PrintAdsorbateBondStatus=TRUE;
  PrintAdsorbateUreyBradleyStatus=TRUE;
  PrintAdsorbateBendStatus=TRUE;
  PrintAdsorbateTorsionStatus=TRUE;
  PrintAdsorbateImproperTorsionStatus=TRUE;
  PrintAdsorbateBondBondStatus=TRUE;
  PrintAdsorbateBondBendStatus=TRUE;
  PrintAdsorbateBendBendStatus=TRUE;
  PrintAdsorbateBondTorsionStatus=TRUE;
  PrintAdsorbateBendTorsionStatus=TRUE;
  PrintAdsorbateIntraVDWStatus=TRUE;
  PrintAdsorbateIntraChargeChargeStatus=TRUE;
  PrintAdsorbateIntraChargeBondDipoleStatus=TRUE;
  PrintAdsorbateIntraBondDipoleBondDipoleStatus=TRUE;

  PrintCationBondStatus=TRUE;
  PrintCationUreyBradleyStatus=TRUE;
  PrintCationBendStatus=TRUE;
  PrintCationTorsionStatus=TRUE;
  PrintCationImproperTorsionStatus=TRUE;
  PrintCationBondBondStatus=TRUE;
  PrintCationBondBendStatus=TRUE;
  PrintCationBendBendStatus=TRUE;
  PrintCationBondTorsionStatus=TRUE;
  PrintCationBendTorsionStatus=TRUE;
  PrintCationIntraVDWStatus=TRUE;
  PrintCationIntraChargeChargeStatus=TRUE;
  PrintCationIntraChargeBondDipoleStatus=TRUE;
  PrintCationIntraBondDipoleBondDipoleStatus=TRUE;

  PrintInterVDWStatus=TRUE;
  PrintInterChargeChargeStatus=TRUE;
  PrintInterChargeBondDipoleStatus=TRUE;
  PrintInterBondDipoleBondDipoleStatus=TRUE;

  PrintFrameworkAdsorbateVDWStatus=TRUE;
  PrintFrameworkAdsorbateChargeChargeStatus=TRUE;
  PrintFrameworkAdsorbateChargeBondDipoleStatus=TRUE;
  PrintFrameworkAdsorbateBondDipoleBondDipoleStatus=TRUE;

  PrintFrameworkCationVDWStatus=TRUE;
  PrintFrameworkCationChargeChargeStatus=TRUE;
  PrintFrameworkCationChargeBondDipoleStatus=TRUE;
  PrintFrameworkCationBondDipoleBondDipoleStatus=TRUE;


  // first pass to get NumberOfSystems and NumberOfComponents etc.
  // also loof for a binary restart-file, in that case immediately return
  if(!(FilePtr=fopen(inputfilename,"r")))
  {
    printf("Error opening input-file '%s' (routine int ReadInputFile(char *inputfilename))\n",inputfilename);
    exit(1); 
  }


  NumberOfSystems=0;
  NumberOfComponents=0;
  while(ReadLine(line,16384,FilePtr))
  {
    // extract first word
    strcpy(keyword,"keyword");
    sscanf(line,"%s%[^\n]",keyword,arguments);
    sscanf(arguments,"%s",firstargument);

    if(strcasecmp("ReducedUnits",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) UseReducedUnits=TRUE;
      if(strcasecmp("no",firstargument)==0) UseReducedUnits=FALSE;
    }

    if(strcasecmp("Dimension",keyword)==0) sscanf(arguments,"%d",&Dimension);

    // read statements for the number of systems
    if(strcasecmp("Box",keyword)==0)
    {
      if(sscanf(arguments,"%d",&CurrentSystem)) 
      {
        if(CurrentSystem<0)
        {
          printf("Input error in 'Box'-command: the number of the system is incorrect (%d given, 0 <= system-number\n",CurrentSystem);
          exit(0);
        }
        NumberOfSystems++;
      }
    }

    if(strcasecmp("BoxMatrix",keyword)==0)
    {
      if(sscanf(arguments,"%d",&CurrentSystem)) 
      {
        if(CurrentSystem<0)
        {
          printf("Input error in 'BoxMatrix'-command: the number of the system is incorrect (%d given, 0 <= system-number\n",CurrentSystem);
          exit(0);
        }
        NumberOfSystems++;
      }
    }

    if(strcasecmp("Framework",keyword)==0)
    {
      if(sscanf(arguments,"%d",&CurrentSystem)) 
      {
        if(CurrentSystem<0)
        {
          printf("Input error in 'Framework'-command: the number of the system is incorrect (%d given, 0 <= system-number\n",CurrentSystem);
          exit(0);
        }
        NumberOfSystems++;
      }
    }

    // read statement for the number of components
    if(strcasecmp("Component",keyword)==0)
    {
      if(sscanf(arguments,"%d %*s %s",&CurrentComponent,string)) 
      {
        if(CurrentComponent<0)
        {
          printf("Input error in 'Component'-command: the number of the component is incorrect (%d given, 0 <= component-number)\n",CurrentComponent);
          exit(0);
        }
        NumberOfComponents++;
      }
    }

    // if restarted from a binary restart file we  can skip everything, and read
    // the full system status from that binary restart-file
    if(strcasecmp("ContinueAfterCrash",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0)
      {
        ContinueAfterCrash=TRUE;

        // no need to reed the rest of the input file
        // instead read the binary status file
        fclose(FilePtr);
        ReadBinaryRestartFiles();
        return 0;
      }
    }
  }
  fclose(FilePtr);

  // set units, either reduced or real units
  SetSimulationUnits();

  // 'NumberOfSystems' and 'NumberOfComponents' are known
  // allocate memory for the components that depends on 'NumberOfSystems' and
  // 'NumberOfComponents', and set some defaults if needed

  // allocate memory
  AllocateSimulationMemory();
  AllocateWarningsMemory();
  AllocateEquationOfStateMemory();
  AllocateMovieMemory();
  AllocateRecrossingMemory();
  AllocateMinimizationMemory();

  read_frameworks=(int*)calloc(NumberOfSystems,sizeof(REAL));
  InitializeBox=(int*)calloc(NumberOfSystems,sizeof(REAL));

  for(i=0;i<NumberOfSystems;i++)
  {
    AdsorbateFixedInitialization[i]=FREE;
    CationFixedInitialization[i]=FREE;
  }

  for(i=0;i<NumberOfSystems;i++)
    for(j=0;j<NumberOfComponents;j++)
      ComputeFugacityCoefficient[i][j]=TRUE;

  Framework=(FRAMEWORK_COMPONENT*)calloc(NumberOfSystems,sizeof(FRAMEWORK_COMPONENT));

  crystallographic_stats=(CRYSTALLOGRAPHIC_STATISTICS*)calloc(NumberOfSystems,sizeof(CRYSTALLOGRAPHIC_STATISTICS));

  for(i=0;i<NumberOfSystems;i++)
  {
    Framework[i].NumberOfFrameworks=0;
    Framework[i].SurfaceAreaSamplingPointsPerShere=5;
    strcpy(Framework[i].SurfaceAreaProbeAtom,"N_n2");
    Framework[i].SurfaceAreaProbeDistance=pow(2.0,1.0/6.0);
  }

  UnitCellSize=(VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));
  DelaunayUnitCellSize = (VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));
  NumberOfUnitCells=(INT_VECTOR3*)calloc(NumberOfSystems,sizeof(INT_VECTOR3));
  UnitCellBox=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  InverseUnitCellBox=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  Lowenstein=(int*)calloc(NumberOfSystems,sizeof(int));

  NumberOfIsothermPressures=1;
  for(i=0;i<NumberOfSystems;i++)
  {
    read_frameworks[i]=FALSE;
    BarrierNormal[i].x=0.0;
    BarrierNormal[i].y=0.0;
    BarrierNormal[i].z=0.0;
    MaxBarrierTime[i]=10.0;
    NumberOfVelocities[i]=5;
    BoundaryCondition[i]=UNINITIALIZED_BOUNDARY_CONDITION;
    Ensemble[i]=NVT;
    InitEnsemble[i]=NVE;
    RunEnsemble[i]=NVE;

    Box[i].ax=Box[i].by=Box[i].cz=25.0;
    NumberOfUnitCells[i].x=1;
    NumberOfUnitCells[i].y=1;
    NumberOfUnitCells[i].z=1;
    AlphaAngle[i]=90.0*M_PI/180.0;
    BetaAngle[i]=90.0*M_PI/180.0;
    GammaAngle[i]=90.0*M_PI/180.0;
    Lowenstein[i]=TRUE;
  }

  if(NumberOfComponents>0)
    Components=(COMPONENT*)calloc(NumberOfComponents,sizeof(COMPONENT));

  MaxNumberOfAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  TotalAdsorbateMoleculesAdded=(int*)calloc(NumberOfSystems,sizeof(int));

  MaxNumberOfCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));

  NumberOfAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));

  NumberOfAtomsPerSystem=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfChargesPerSystem=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfBondDipolesPerSystem=(int*)calloc(NumberOfSystems,sizeof(int));

  for(i=0;i<NumberOfComponents;i++)
  {
    Components[i].IdealGasRosenbluthWeight=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].IdealGasTotalEnergy=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].PartialPressure=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].FugacityCoefficient=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].BulkFluidDensity=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].Compressibility=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MolFraction=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].AmountOfExcessMolecules=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    Components[i].CreateNumberOfMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].NumberOfMolecules=(int*)calloc(NumberOfSystems,sizeof(int));

    Components[i].FractionalMolecule=(int*)calloc(NumberOfSystems,sizeof(int));

    Components[i].MOLEC_PER_UC_TO_MOL_PER_KG=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOLEC_PER_UC_TO_CC_STP_G=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOLEC_PER_UC_TO_CC_STP_CC=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOL_PER_KG_TO_CC_STP_G=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOL_PER_KG_TO_CC_STP_CC=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    Components[i].IdentityChanges=(int*)calloc(NumberOfComponents,sizeof(int));
    Components[i].GibbsIdentityChanges=(int*)calloc(NumberOfComponents,sizeof(int));

    Components[i].BlockPockets=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].ComputeFreeEnergyProfile=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].BlockPocketsFilename=(char(*)[256])calloc(NumberOfSystems,sizeof(char[256]));
    Components[i].NumberOfBlockCenters=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].BlockDistance=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].BlockCenters=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));

    for(j=0;j<NumberOfSystems;j++)
    {
      Components[i].BlockPockets[j]=FALSE;
      Components[i].ComputeFreeEnergyProfile[j]=FALSE;
      Components[i].IdealGasRosenbluthWeight[j]=1.0;
      Components[i].IdealGasTotalEnergy[j]=0.0;
      Components[i].FugacityCoefficient[j]=-1.0;
      Components[i].MolFraction[j]=1.0;
    }

    strcpy(Components[i].MoleculeDefinition,"TraPPE");
    strcpy(Components[i].Name,"");
    Components[i].Swapable=FALSE;
    Components[i].Widom=FALSE;

    // initialize box restriction to false (the full box is selected)
    Components[i].RestrictMoves=FALSE;

    // initialize the first box to all, the other 3 boxes to empty
    Components[i].RestrictMovesToBox=FALSE;
    Components[i].BoxAxisABC_Min.x=Components[i].BoxAxisABC_Min.y=Components[i].BoxAxisABC_Min.z=0.0;
    Components[i].BoxAxisABC_Max.x=Components[i].BoxAxisABC_Max.y=Components[i].BoxAxisABC_Max.z=1.0;
    Components[i].BoxAxisABC_Min2.x=Components[i].BoxAxisABC_Min2.y=Components[i].BoxAxisABC_Min2.z=1.0;
    Components[i].BoxAxisABC_Max2.x=Components[i].BoxAxisABC_Max2.y=Components[i].BoxAxisABC_Max2.z=0.0;
    Components[i].BoxAxisABC_Min3.x=Components[i].BoxAxisABC_Min3.y=Components[i].BoxAxisABC_Min3.z=1.0;
    Components[i].BoxAxisABC_Max3.x=Components[i].BoxAxisABC_Max3.y=Components[i].BoxAxisABC_Max3.z=0.0;
    Components[i].BoxAxisABC_Min4.x=Components[i].BoxAxisABC_Min4.y=Components[i].BoxAxisABC_Min4.z=1.0;
    Components[i].BoxAxisABC_Max4.x=Components[i].BoxAxisABC_Max4.y=Components[i].BoxAxisABC_Max4.z=0.0;

    Components[i].RestrictMovesToPrisms=FALSE;
    for(j=0;j<MAX_NUMBER_OF_PRISMS;j++)
    {
      Components[i].RestrictMovesToPrism[j]=FALSE;
      Components[i].RestrictPrismABC_Min[j].x=0.0;
      Components[i].RestrictPrismABC_Max[j].x=1.0;
    }

    Components[i].RestrictMovesToCylinders=FALSE;
    for(j=0;j<MAX_NUMBER_OF_CYLINDERS;j++)
    {
      Components[i].RestrictMovesToCylinder[j]=FALSE;
      Components[i].RestrictCylinderABC_Min[j].x=0.0;
      Components[i].RestrictCylinderABC_Max[j].x=1.0;
      Components[i].RestrictCylinderABC_Min[j].y=0.0;
      Components[i].RestrictCylinderABC_Max[j].y=1.0;
      Components[i].RestrictCylinderABC_Min[j].z=0.0;
      Components[i].RestrictCylinderABC_Max[j].z=1.0;
      Components[i].RestrictCylinderCenter[j].x=0.0;
      Components[i].RestrictCylinderCenter[j].y=0.5;
      Components[i].RestrictCylinderCenter[j].z=0.5;
      Components[i].RestrictCylinderDirection[j]=X_DIR;
      Components[i].RestrictCylinderRadius[j]=5.0;
    }

    Components[i].RestrictMovesToSpheres=FALSE;
    for(j=0;j<MAX_NUMBER_OF_SPHERES;j++)
    {
      Components[i].RestrictMovesToSphere[j]=FALSE;
      Components[i].RestrictSphereCenter[j].x=0.5;
      Components[i].RestrictSphereCenter[j].y=0.5;
      Components[i].RestrictSphereCenter[j].z=0.5;
      Components[i].RestrictSphereRadius[j]=5.0;
    }

    Components[i].RuizMonteroFactor=1.0;
    Components[i].UmbrellaFactor=1.0;

    Components[i].Intra14VDWScalingValue=1.0;
    Components[i].Intra14ChargeChargeScalingValue=1.0;

    Components[i].AnisotropicType=ANISOTROPIC_MID_POINT;
    Components[i].BiasingDirection=A_DIR;
  }

  Alpha=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  kvec=(INT_VECTOR3*)calloc(NumberOfSystems,sizeof(INT_VECTOR3));
  NumberOfKVectors=(int*)calloc(NumberOfSystems,sizeof(int));
  ReciprocalCutOffSquared=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  // second pass to get the number of frameworks per system
  if(!(FilePtr=fopen(inputfilename,"r")))
  {
    printf("Error opening input-file '%s'\n",inputfilename);
    exit(1);
  }

  CurrentSystem=0;
  while(ReadLine(line,16384,FilePtr))
  {
   // extract first word
    strcpy(keyword,"keyword");
    sscanf(line,"%s%[^\n]",keyword,arguments);
    sscanf(arguments,"%s",firstargument);

    if(strcasecmp("Framework",keyword)==0) 
    {
      sscanf(arguments,"%d",&CurrentSystem);

      if(CurrentSystem>=NumberOfSystems)
      {
        printf("Input error in 'Framework'-command: the number of the system is incorrect (%d given, 0 <= system-number < %d)\n",CurrentSystem,NumberOfSystems);
        exit(0);
      }
    }
    if(strcasecmp("FrameworkName",keyword)==0)
      Framework[CurrentSystem].NumberOfFrameworks++;

    if(strcasecmp("ExternalPressure",keyword)==0)
    {
      arg_pointer=arguments;
      NumberOfIsothermPressures=0;
      while(sscanf(arg_pointer,"%lf%n",&tempd,&n)==1)
      {
        arg_pointer+=n;
        NumberOfIsothermPressures++;
      }
    }
  }
  fclose(FilePtr);

  AllocateThermoBaroStatMemory();
  AllocateSampleMemory();

  NumberOfIsothermPressures=1;
  for(i=0;i<NumberOfSystems;i++)
  {
    ReciprocalCutOffSquared[i]=-1.0;

    therm_baro_stats.ExternalTemperature[i]=298.0;

    therm_baro_stats.ExternalSurfaceTension[i].x=0.0;
    therm_baro_stats.ExternalSurfaceTension[i].x=0.0;
    therm_baro_stats.ExternalSurfaceTension[i].x=0.0;

    Movies[i]=FALSE;
    WriteMoviesEvery[i]=1000;

    // sampling the radial distribution function (RDF)
    ComputeRDF[i]=FALSE;
    WriteRDFEvery[i]=5000;
    RDFHistogramSize[i]=500;
    RDFRange[i]=12.0;

    // sampling the number-of-molecules histogram
    ComputeNumberOfMoleculesHistogram[i]=FALSE;
    WriteNumberOfMoleculesHistogramEvery[i]=5000;
    NumberOfMoleculesRange[i]=500.0;
    NumberOfMoleculesHistogramSize[i]=1000;

    // sampling position histograms/free energies
    ComputePositionHistogram[i]=FALSE;
    WritePositionHistogramEvery[i]=5000;
    PositionHistogramSize[i]=1000;
    PositionHistogramMappingType[i]=ABC_MAPPING;

    // samples the free energy profiles in a,b,c directions
    ComputeFreeEnergyProfile[i]=FALSE;
    WriteFreeEnergyProfileEvery[i]=5000;
    FreeEnergyHistogramSize[i]=1800;
    FreeEnergyMappingType[i]=ABC_MAPPING;

    // sampling the pore-size distribution (PSD)
    ComputePSDHistogram[i]=FALSE;
    WritePSDHistogramEvery[i]=1000;
    PSDHistogramSize[i]=200;
    PSDRange[i]=10.0;

    // sampling the end-to-end histograms
    ComputeEndToEndDistanceHistogram[i]=FALSE;
    WriteEndToEndDistanceHistogramEvery[i]=5000;
    EndToEndHistogramSize[i]=1000;
    EndToEndRange[i]=30.0;

    // sampling the energy histogram
    WriteEnergyHistogramEvery[i]=5000;
    EnergyHistogramLowerLimit[i]=-1e4;
    EnergyHistogramUpperLimit[i]=0;
    EnergyHistogramSize[i]=1000;

    // sampling the thermodynamic factor
    ComputeThermoDynamicFactor[i]=FALSE;
    WriteThermoDynamicFactorEvery[i]=5000;

    // sampling the inter-framework spacing histogram
    ComputeFrameworkSpacingHistogram[i]=FALSE;
    WriteFrameworkSpacingHistogramEvery[i]=5000;
    FrameworkSpacingHistogramSize[i]=1000;
    FrameworkSpacingRange[i]=10.0;

    // sampling histograms of the residence times
    ComputeResidenceTimes[i]=FALSE;
    WriteResidenceTimesEvery[i]=5000;
    ResidenceTimesHistogramSize[i]=200;
    RangeResidenceTimes[i]=500.0; // 0.5 nanosecond


    // sampling histograms of the distance between 2 selected atoms
    ComputeDistanceHistograms[i]=FALSE;
    WriteDistanceHistogramsEvery[i]=5000;

    // sampling histograms of the bend angle between 3 selected atoms
    ComputeBendAngleHistograms[i]=FALSE;
    WriteBendAngleHistogramsEvery[i]=5000;

    // sampling histograms of the angle between two planes (each formed by 3 chosen atoms)
    ComputeDihedralAngleHistograms[i]=FALSE;
    WriteDihedralAngleHistogramsEvery[i]=5000;

    // sampling histograms of the dihedral angle between 4 selected atoms
    ComputeAngleBetweenPlanesHistograms[i]=FALSE;
    WriteAngleBetweenPlanesHistogramsEvery[i]=5000;

    // sampling molecular properties (bond distance, bend angle, dihedral angle)
    ComputeMoleculeProperties[i]=FALSE;
    WriteMoleculePropertiesEvery[i]=5000;
    BondLengthHistogramSize[i]=1000;
    BendAngleHistogramSize[i]=1000;
    DihedralHistogramSize[i]=1000;
    BondLengthRange[i]=4.0;
    BendAngleRange[i]=180.0;
    DihedralRange[i]=360.0;

    // sampling the IR spectra (spacings: 2048, 4196, 8192, 16384, 32768 points)
    ComputeInfraRedSpectra[i]=FALSE;
    WriteInfraRedSpectraEvery[i]=5000;

    // sampling the mean-squared displacement using a modified order-N algorithm
    ComputeMSDOrderN[i]=FALSE;
    WriteMSDOrderNEvery[i]=5000;
    SampleMSDOrderNEvery[i]=1;

    // sampling the velocity autocorrelation function using a modified order-N algorithm
    ComputeVACFOrderN[i]=FALSE;
    WriteVACFOrderNEvery[i]=5000;
    SampleVACFOrderNEvery[i]=5;

    // sampling of the rotational velocity autocorrelation function using a modified order-N algorithm
    ComputeRVACFOrderN[i]=FALSE;
    WriteRVACFOrderNEvery[i]=5000;
    SampleRVACFOrderNEvery[i]=5;

    // sampling of the molecular orientation autocorrelation function using a modified order-N algorithm
    ComputeMolecularOrientationOrderN[i]=FALSE;
    WriteMolecularOrientationOrderNEvery[i]=5000;
    SampleMolecularOrientationOrderNEvery[i]=5;

    // sampling the mean-square displacement function using a conventional algorithm
    ComputeMSD[i]=FALSE;
    WriteMSDEvery[i]=5000;
    SampleMSDEvery[i]=5;

    // sampling the velocity autocorrelation function using a conventional algorithm
    ComputeVACF[i]=FALSE;
    WriteVACFEvery[i]=5000;
    SampleVACFEvery[i]=5;

    // sampling the 3D histograms of position (i.e. 3D free energy)
    ComputeDensityProfile3DVTKGrid[i]=FALSE;
    WriteDensityProfile3DVTKGridEvery[i]=5000;

    // samples the cation sites and adsorption sites
    ComputeCationAndAdsorptionSites[i]=FALSE;
    WriteCationAndAdsorptionSitesEvery[i]=5000;

    // samples initial configurations for the tranmission coefficient (dcTST)
    WritedcTSTSnapShotsToFile[i]=FALSE;
    WritedcTSTSnapShotsEvery[i]=100;

    // samples pressure and stress
    ComputeMolecularPressure[i]=FALSE;
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    if(Framework[i].NumberOfFrameworks>0)
    {
      Framework[i].NumberOfCoreShells=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

      Framework[i].Atoms=(ATOM**)calloc(Framework[i].NumberOfFrameworks,sizeof(ATOM*));
      Framework[i].NumberOfAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfFixedAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfFreeAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfCharges=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfUnitCellAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].FrameworkProbability=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));

      Framework[i].AtomsAsymmetric=(FRAMEWORK_ASYMMETRIC_ATOM**)calloc(Framework[i].NumberOfFrameworks,sizeof(FRAMEWORK_ASYMMETRIC_ATOM*));
      Framework[i].NumberOfAsymmetricAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

      Framework[i].Name=(char(*)[256])calloc(Framework[i].NumberOfFrameworks,sizeof(char[256]));
      Framework[i].FrameworkModels=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].ShiftUnitCell=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
      Framework[i].Asymmetric=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].SpaceGroupIdentifier=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].CalculateSpaceGroup=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].RemoveAtomNumberCodeFromLabel=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].AddAtomNumberCodeToLabel=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].InputFileType=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].SurfaceAreas=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));
      Framework[i].FrameworkDensityPerComponent=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));
      Framework[i].FrameworkMassPerComponent=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));

      Framework[i].NumberOfCitations=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].CitationInformation=(CITATION_INFORMATION**)calloc(Framework[i].NumberOfFrameworks,sizeof(CITATION_INFORMATION*));

      Framework[i].Connectivity=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].Neighbours=(int***)calloc(Framework[i].NumberOfFrameworks,sizeof(int**));

      // allocate first dimension of bonds
      Framework[i].NumberOfBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].Bonds=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));
      Framework[i].BondType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].BondArguments=(REAL(**)[MAX_BOND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of bond-dipoles
      Framework[i].NumberOfBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].BondDipoles=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));
      Framework[i].BondDipoleMagnitude=(REAL**)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL*));

      // allocate first dimension of Urey-Bradleys
      Framework[i].NumberOfUreyBradleys=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfUreyBradleys=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].UreyBradleys=(TRIPLE**)calloc(Framework[i].NumberOfFrameworks,sizeof(TRIPLE*));
      Framework[i].UreyBradleyType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].UreyBradleyArguments=(REAL(**)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of bends
      Framework[i].NumberOfBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].Bends=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].BendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].BendArguments=(REAL(**)[MAX_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of inversion-bends
      Framework[i].NumberOfInversionBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfInversionBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].InversionBends=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].InversionBendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].InversionBendArguments=(REAL(**)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of torsions
      Framework[i].NumberOfTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].Torsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].TorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].TorsionArguments=(REAL(**)[MAX_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of improper torsions
      Framework[i].NumberOfImproperTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfImproperTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].ImproperTorsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].ImproperTorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].ImproperTorsionArguments=(REAL(**)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of out-of-plane
      Framework[i].NumberOfOutOfPlanes=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfOutOfPlanes=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].OutOfPlanes=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].OutOfPlaneType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].OutOfPlaneArguments=(REAL(**)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of bond-bond
      Framework[i].NumberOfBondBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBondBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].BondBonds=(TRIPLE**)calloc(Framework[i].NumberOfFrameworks,sizeof(TRIPLE*));
      Framework[i].BondBondType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].BondBondArguments=(REAL(**)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of bond-bend
      Framework[i].NumberOfBondBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBondBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].BondBends=(TRIPLE**)calloc(Framework[i].NumberOfFrameworks,sizeof(TRIPLE*));
      Framework[i].BondBendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].BondBendArguments=(REAL(**)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of bend-bend
      Framework[i].NumberOfBendBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBendBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].BendBends=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].BendBendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].BendBendArguments=(REAL(**)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of stretch-torsion
      Framework[i].NumberOfBondTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBondTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].BondTorsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].BondTorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].BondTorsionArguments=(REAL(**)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));

      // allocate first dimension of bend-torsion
      Framework[i].NumberOfBendTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfBendTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].BendTorsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
      Framework[i].BendTorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].BendTorsionArguments=(REAL(**)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));

      Framework[i].NumberOfIntraVDW=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfExcludedIntraVDW=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

      Framework[i].NumberOfIntraCharges=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfIntraChargeBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfIntraBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));


      Framework[i].NumberOfExcludedIntraChargeCharge=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfExcludedIntraChargeCharge=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].ExcludedIntraChargeCharge=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));

      Framework[i].NumberOfExcludedIntraChargeBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfExcludedIntraChargeBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].ExcludedIntraChargeBondDipole=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));

      Framework[i].NumberOfExcludedIntraBondDipoleBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfExcludedIntraBondDipoleBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].ExcludedIntraBondDipoleBondDipole=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));

      FrameworkFixedInitialization[i]=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      NumberOfFixedFrameworkAtoms[i]=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      NumberOfActiveFrameworkAtoms[i]=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      FixedFrameworkAtoms[i]=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      ActiveFrameworkAtoms[i]=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));

      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        FrameworkFixedInitialization[i][j]=-1;
        Framework[i].SpaceGroupIdentifier[j]=1;
        Framework[i].FrameworkProbability[j]=1.0/Framework[i].NumberOfFrameworks;
        Framework[i].RemoveAtomNumberCodeFromLabel[j]=FALSE;
        Framework[i].AddAtomNumberCodeToLabel[j]=FALSE;
        Framework[i].InputFileType[j]=CIF_FORMAT;
        Framework[i].FrameworkProbability[j]=1.0;
      }

      Framework[i].PoreSizeDistributionProbeDistance=1.0;
      Framework[i].AnisotropicType=ANISOTROPIC_MID_POINT;
      Framework[i].ForceSpaceGroupDetection=FALSE;
      Framework[i].ReadCIFAsCartesian=FALSE;
    }
    else
    {
      Framework[i].Name=(char(*)[256])calloc(1,sizeof(char[256]));
      Framework[i].NumberOfAtoms=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfCharges=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfAtoms[0]=0;
      Framework[i].NumberOfCharges[0]=0;
      Framework[i].NumberOfAsymmetricAtoms=(int*)calloc(1,sizeof(int));
      Framework[i].RemoveAtomNumberCodeFromLabel=(int*)calloc(1,sizeof(int));
      Framework[i].AddAtomNumberCodeToLabel=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfAsymmetricAtoms[0]=0;
      Framework[i].FrameworkModels=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

      Framework[i].PoreSizeDistributionProbeDistance=1.0;
      Framework[i].AnisotropicType=ANISOTROPIC_MID_POINT;

      Framework[i].FrameworkModel=NONE;
    }
  }

  // allocate memory
  AllocateStatisticsMemory();

  // final pass, most memory is now already allocated
  if(!(FilePtr=fopen(inputfilename,"r")))
  {
    printf("Error opening input-file '%s'\n",inputfilename);
    exit(1); 
  }

  CurrentComponent=0;
  LineNumber=0;
  while(ReadLine(line,16384,FilePtr))
  {
   LineNumber++;

   // extract first word
    strcpy(keyword,"keyword");
    sscanf(line,"%s%[^\n]",keyword,arguments);
    sscanf(arguments,"%s",firstargument);

    // read the simulation type
    if(strcasecmp("SimulationType",keyword)==0)
    {
      if(strcasecmp("MonteCarlo",firstargument)==0) SimulationType=MONTE_CARLO;
      if(strcasecmp("MC",firstargument)==0) SimulationType=MONTE_CARLO;
      if(strcasecmp("MolecularDynamics",firstargument)==0) SimulationType=MOLECULAR_DYNAMICS;
      if(strcasecmp("MD",firstargument)==0) SimulationType=MOLECULAR_DYNAMICS;
      if(strcasecmp("Spectra",firstargument)==0) SimulationType=SPECTRA;
      if(strcasecmp("PhononDispersion",firstargument)==0) SimulationType=PHONON_DISPERSION;
      if(strcasecmp("Minimization",firstargument)==0) SimulationType=MINIMIZATION;
      if(strcasecmp("GlobalMinimization",firstargument)==0) SimulationType=GLOBAL_MINIMIZATION;
      if(strcasecmp("Visualization",firstargument)==0) SimulationType=VISUALIZATION;
      if((strcasecmp("PoreSizeDistribution",firstargument)==0)||
         (strcasecmp("PSD",firstargument)==0)) SimulationType=PORE_SIZE_DISTRIBUTION;
      if(strcasecmp("SurfaceArea",firstargument)==0) SimulationType=SURFACE_AREA;
      if(strcasecmp("BarrierRecrossing",firstargument)==0) SimulationType=BARRIER_RECROSSING;
      if(strcasecmp("RuizMontero",firstargument)==0) SimulationType=BARRIER_RUIZ_MONTERO;
      if(strcasecmp("FreeEnergy",firstargument)==0) SimulationType=FREE_ENERGY;
      if(strcasecmp("MakeGrid",firstargument)==0) SimulationType=MAKE_GRID;
      if(strcasecmp("MakeASCIGrid",firstargument)==0) SimulationType=MAKE_ASCI_GRID;
      if(strcasecmp("Numerical",firstargument)==0) SimulationType=NUMERICAL;
      if(strcasecmp("PlotPotential",firstargument)==0) SimulationType=PLOT_POTENTIAL;
      if(strcasecmp("Status",firstargument)==0) SimulationType=STATUS;
      if(strcasecmp("GenerateFramework",firstargument)==0) SimulationType=GENERATE_FRAMEWORK;
    }

    if(strcasecmp("RandomSeed",keyword)==0) sscanf(arguments,"%lu",&seed);

    // read simulation parameters
    if(strcasecmp("TimeStep",keyword)==0) sscanf(arguments,"%lf",&DeltaT);
    if(strcasecmp("CutOff",keyword)==0) sscanf(arguments,"%lf",&CutOffVDW);
    if(strcasecmp("CutOffVDW",keyword)==0) sscanf(arguments,"%lf",&CutOffVDW);
    if(strcasecmp("CutOffVDWSwitch",keyword)==0) sscanf(arguments,"%lf",&CutOffVDWSwitch);
    if(strcasecmp("CutOffCoulomb",keyword)==0) sscanf(arguments,"%lf",&CutOffChargeCharge);
    if(strcasecmp("CutOffCoulombSwitch",keyword)==0) sscanf(arguments,"%lf",&CutOffChargeChargeSwitch);
    if(strcasecmp("CutOffChargeCharge",keyword)==0) sscanf(arguments,"%lf",&CutOffChargeCharge);
    if(strcasecmp("CutOffChargeChargeSwitch",keyword)==0) sscanf(arguments,"%lf",&CutOffChargeChargeSwitch);
    if(strcasecmp("CutOffChargeBondDipole",keyword)==0) sscanf(arguments,"%lf",&CutOffChargeBondDipole);
    if(strcasecmp("CutOffChargeBondDipoleSwitch",keyword)==0) sscanf(arguments,"%lf",&CutOffChargeBondDipoleSwitch);
    if(strcasecmp("CutOffBondDipoleBondDipole",keyword)==0) sscanf(arguments,"%lf",&CutOffBondDipoleBondDipole);
    if(strcasecmp("CutOffBondDipoleBondDipoleSwitch",keyword)==0) sscanf(arguments,"%lf",&CutOffBondDipoleBondDipoleSwitch);
    if(strcasecmp("CutOffIons",keyword)==0) sscanf(arguments,"%lf",&CutOffIons);

    // read the amount of cycles
    if(strcasecmp("NumberOfCycles",keyword)==0) sscanf(arguments,"%lld",&NumberOfCycles);
    if(strcasecmp("NumberOfInitializationCycles",keyword)==0) sscanf(arguments,"%lld",&NumberOfInitializationCycles);
    if(strcasecmp("NumberOfEquilibrationCycles",keyword)==0) sscanf(arguments,"%lld",&NumberOfEquilibrationCycles);
    if(strcasecmp("NumberOfVelocityScalingCycles",keyword)==0) sscanf(arguments,"%lld",&NumberOfVelocityScalingCycles);


    // read restart-options
    if(strcasecmp("RestartFile",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Restart=TRUE;
      if(strcasecmp("no",firstargument)==0) Restart=FALSE;
    }
    if(strcasecmp("ContinueAfterCrash",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) 
      {
        ContinueAfterCrash=TRUE;

        // no need to reed the rest of the input file
        // instead read the binary status file
        fclose(FilePtr); 
        ReadBinaryRestartFiles();
        return 0;
      }
      if(strcasecmp("no",firstargument)==0) ContinueAfterCrash=FALSE;
    }
    if(strcasecmp("WriteBinaryRestartFileEvery",keyword)==0) sscanf(arguments,"%d",&WriteBinaryRestartFileEvery);
    if(strcasecmp("RestartStyle",keyword)==0)
    {
      if(strcasecmp("BIGMAC",firstargument)==0) RestartStyle=BIGMAC_STYLE;
      if(strcasecmp("RASPA",firstargument)==0) RestartStyle=RASPA_STYLE;
    }
    if(strcasecmp("ReinitializeVelocities",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ReinitializeVelocities=TRUE;
      if(strcasecmp("no",firstargument)==0) ReinitializeVelocities=FALSE;
    }  


    // read charge parameters
    if(strcasecmp("ChargeMethod",keyword)==0)
    {
      if(strcasecmp("None",firstargument)==0) ChargeMethod=NONE;
      if(strcasecmp("Ewald",firstargument)==0) ChargeMethod=EWALD;
      if(strcasecmp("Coulomb",firstargument)==0) ChargeMethod=TRUNCATED_COULOMB;
      if(strcasecmp("WolfsMethod",firstargument)==0) ChargeMethod=WOLFS_METHOD;
      if(strcasecmp("WolfsMethodDamped",firstargument)==0) ChargeMethod=WOLFS_METHOD_DAMPED;
      if(strcasecmp("WolfsMethodDampedFG",firstargument)==0) ChargeMethod=WOLFS_METHOD_DAMPED_FG;
      if(strcasecmp("CoulombTruncated",firstargument)==0) ChargeMethod=TRUNCATED_COULOMB;
      if(strcasecmp("CoulombShifted",firstargument)==0) ChargeMethod=SHIFTED_COULOMB;
      if(strcasecmp("CoulombSmoothed",firstargument)==0) ChargeMethod=SMOOTHED_COULOMB;
    }
    if(strcasecmp("ChargeFromChargeEquilibration",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ChargeFromChargeEquilibration=TRUE;
      if(strcasecmp("no",firstargument)==0) ChargeFromChargeEquilibration=FALSE;
    }  
    if(strcasecmp("ChargeEquilibrationMethod",keyword)==0)
    {
      if(strcasecmp("RAPPE_GODDARD",firstargument)==0) ChargeEquilibrationMethod=RAPPE_GODDARD;
      if(strcasecmp("WILMER_SNURR",firstargument)==0) ChargeEquilibrationMethod=WILMER_SNURR;
      if(strcasecmp("DING_YAZAYDIN_DUBBELDAM",firstargument)==0) ChargeEquilibrationMethod=DING_YAZAYDIN_DUBBELDAM;
    }  
    if(strcasecmp("SymmetrizeFrameworkCharges",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) SymmetrizeFrameworkCharges=TRUE;
      if(strcasecmp("no",firstargument)==0) SymmetrizeFrameworkCharges=FALSE;
    }
    if(strcasecmp("UseChargesFromMOLFile",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) UseChargesFromMOLFile=TRUE;
      if(strcasecmp("no",firstargument)==0) UseChargesFromMOLFile=FALSE;
    }
    if(strcasecmp("UseChargesFromCIFFile",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) UseChargesFromCIFFile=TRUE;
      if(strcasecmp("no",firstargument)==0) UseChargesFromCIFFile=FALSE;
    }
    if(strcasecmp("EwaldPrecision",keyword)==0) sscanf(arguments,"%lf",&EwaldPrecision),EwaldAutomatic=TRUE;
    if(strcasecmp("EwaldParameters",keyword)==0) sscanf(arguments,"%lf %d %d %d",&Alpha[CurrentSystem],
       &kvec[CurrentSystem].x,&kvec[CurrentSystem].y,&kvec[CurrentSystem].z),EwaldAutomatic=FALSE;
    if(strcasecmp("ReciprocalCutOff",keyword)==0) 
    {
      sscanf(arguments,"%lf",&realinput1);
      ReciprocalCutOffSquared[CurrentSystem]=SQR(realinput1);
    }
    if(strcasecmp("NoReciprocalCutOff",keyword)==0) ReciprocalCutOffSquared[CurrentSystem]=DBL_MAX;
    if(strcasecmp("DielectricConstantOfTheMedium",keyword)==0) sscanf(arguments,"%lf",&DielectricConstantOfTheMedium);

    // read some printing options
    if(strcasecmp("PrintEvery",keyword)==0) sscanf(arguments,"%d",&PrintEvery);
    if(strcasecmp("PrintPropertiesEvery",keyword)==0) sscanf(arguments,"%d",&PrintPropertiesEvery);
    if(strcasecmp("Output",keyword)==0)
    {
      if(strcasecmp("screen",firstargument)==0) Output=OUTPUT_TO_SCREEN;
      if(strcasecmp("file",firstargument)==0) Output=OUTPUT_TO_FILE;
    }
    if(strcasecmp("FileNameAppend",keyword)==0)
      sscanf(arguments,"%s",FileNameAppend);
    if(strcasecmp("PrintForcefieldToOutput",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) WriteForcefieldToOutput=TRUE;
      if(strcasecmp("no",firstargument)==0) WriteForcefieldToOutput=FALSE;
    }
    if(strcasecmp("PrintPseudoAtomsToOutput",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) WritePseudoAtomsToOutput=TRUE;
      if(strcasecmp("no",firstargument)==0) WritePseudoAtomsToOutput=FALSE;
    }
    if(strcasecmp("PrintMoleculeDefinitionToOutput",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) WriteMoleculeDefinitionToOutput=TRUE;
      if(strcasecmp("no",firstargument)==0) WriteMoleculeDefinitionToOutput=FALSE;
    }

    // read the forcefield definition
    if(strcasecmp("ForceField",keyword)==0) sscanf(firstargument,"%s",ForceField);
    if(strcasecmp("ComputePolarization",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputePolarization=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputePolarization=FALSE;
    }
    if(strcasecmp("BackPolarization",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) BackPolarization=TRUE;
      if(strcasecmp("no",firstargument)==0) BackPolarization=FALSE;
    }
    if(strcasecmp("NumberOfBackPolarizationSteps",keyword)==0) sscanf(arguments,"%d",&NumberOfBackPolarizationSteps);
    if(strcasecmp("OmitInterMolecularInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) 
      {
        OmitInterMolecularInteractions=TRUE;
        OmitAdsorbateAdsorbateVDWInteractions=OmitAdsorbateAdsorbateCoulombInteractions=TRUE;
        OmitAdsorbateCationVDWInteractions=OmitAdsorbateCationCoulombInteractions=TRUE;
        OmitCationCationVDWInteractions=OmitCationCationCoulombInteractions=TRUE;
      }
      if(strcasecmp("no",firstargument)==0)
      {
        OmitInterMolecularInteractions=FALSE;
        OmitAdsorbateAdsorbateVDWInteractions=OmitAdsorbateAdsorbateCoulombInteractions=FALSE;
        OmitAdsorbateCationVDWInteractions=OmitAdsorbateCationCoulombInteractions=FALSE;
        OmitCationCationVDWInteractions=OmitCationCationCoulombInteractions=FALSE;
      }
    }
    if(strcasecmp("OmitInterMolecularVDWInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0)
      {
        OmitInterMolecularVDWInteractions=TRUE;
        OmitAdsorbateAdsorbateVDWInteractions=TRUE;
        OmitAdsorbateCationVDWInteractions=TRUE;
        OmitCationCationVDWInteractions=TRUE;
      }
      if(strcasecmp("no",firstargument)==0)
      {
        OmitInterMolecularVDWInteractions=FALSE;
        OmitAdsorbateAdsorbateVDWInteractions=FALSE;
        OmitAdsorbateCationVDWInteractions=FALSE;
        OmitCationCationVDWInteractions=FALSE;
      }
    }
    if(strcasecmp("OmitInterMolecularCoulombInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0)
      {
        OmitInterMolecularCoulombInteractions=TRUE;
        OmitAdsorbateAdsorbateCoulombInteractions=TRUE;
        OmitAdsorbateCationCoulombInteractions=TRUE;
        OmitCationCationCoulombInteractions=TRUE;
      }
      if(strcasecmp("no",firstargument)==0) 
      {
        OmitInterMolecularCoulombInteractions=FALSE;
        OmitAdsorbateAdsorbateCoulombInteractions=FALSE;
        OmitAdsorbateCationCoulombInteractions=FALSE;
        OmitCationCationCoulombInteractions=FALSE;
      }
    }
    if(strcasecmp("OmitAdsorbateAdsorbateVDWInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitAdsorbateAdsorbateVDWInteractions=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitAdsorbateAdsorbateVDWInteractions=FALSE;
    }
    if(strcasecmp("OmitAdsorbateAdsorbateCoulombInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitAdsorbateAdsorbateCoulombInteractions=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitAdsorbateAdsorbateCoulombInteractions=FALSE;
    }
    if(strcasecmp("OmitAdsorbateCationVDWInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitAdsorbateCationVDWInteractions=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitAdsorbateCationVDWInteractions=FALSE;
    }
    if(strcasecmp("OmitAdsorbateCationCoulombInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitAdsorbateCationCoulombInteractions=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitAdsorbateCationCoulombInteractions=FALSE;
    }
    if(strcasecmp("OmitCationCationVDWInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitCationCationVDWInteractions=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitCationCationVDWInteractions=FALSE;
    }
    if(strcasecmp("OmitCationCationCoulombInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitCationCationCoulombInteractions=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitCationCationCoulombInteractions=FALSE;
    }
    if(strcasecmp("OmitAdsorbateAdsorbatePolarization",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitAdsorbateAdsorbatePolarization=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitAdsorbateAdsorbatePolarization=FALSE;
    }
    if(strcasecmp("OmitAdsorbateCationPolarization",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitAdsorbateCationPolarization=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitAdsorbateCationPolarization=FALSE;
    }
    if(strcasecmp("OmitCationCationPolarization",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitCationCationPolarization=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitCationCationPolarization=FALSE;
    }
    if(strcasecmp("OmitIntraFrameworkPolarization",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitIntraFrameworkPolarization=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitIntraFrameworkPolarization=FALSE;
    }
    if(strcasecmp("OmitEwaldFourier",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) OmitEwaldFourier=TRUE;
      if(strcasecmp("no",firstargument)==0) OmitEwaldFourier=FALSE;
    }
    if(strcasecmp("InternalFrameworkLennardJonesInteractions",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) InternalFrameworkLennardJonesInteractions=TRUE;
      if(strcasecmp("no",firstargument)==0) InternalFrameworkLennardJonesInteractions=FALSE;
    }

    if((strcasecmp("Substitute",keyword)==0)||(strcasecmp("Modify",keyword)==0))
    {
      // allocate additional memory for this modification rule
      SubstitutionSingleFrameworkAtom=(char(*)[3][256])realloc(SubstitutionSingleFrameworkAtom,(NumberOfSingleSubstitutionRules+1)*sizeof(char[3][256]));
      SubstitutionSingleFrameworkAtomTypes=(int(*)[3])realloc(SubstitutionSingleFrameworkAtomTypes,(NumberOfSingleSubstitutionRules+1)*sizeof(int[3]));

      // scan the arguments from the input
      sscanf(arguments,"%s %s %s",
          SubstitutionSingleFrameworkAtom[NumberOfSingleSubstitutionRules][0],
          SubstitutionSingleFrameworkAtom[NumberOfSingleSubstitutionRules][1],
          SubstitutionSingleFrameworkAtom[NumberOfSingleSubstitutionRules][2]);

      // increase the amount of modifcation rules by 1
      NumberOfSingleSubstitutionRules++;
    }

    if((strcasecmp("RandomlySubstitute",keyword)==0)||(strcasecmp("RandomlyModify",keyword)==0))
    {
      // allocate additional memory for this modification rule
      SubstitutionFrameworkAtoms=(char(*)[3][256])realloc(SubstitutionFrameworkAtoms,(NumberOfSubstitutionRules+1)*sizeof(char[3][256]));
      SubstitutionFrameworkAtomTypes=(int(*)[3])realloc(SubstitutionFrameworkAtomTypes,(NumberOfSubstitutionRules+1)*sizeof(int[3]));

      // scan the arguments from the input
      sscanf(arguments,"%s %s %s",
          SubstitutionFrameworkAtoms[NumberOfSubstitutionRules][0],
          SubstitutionFrameworkAtoms[NumberOfSubstitutionRules][1],
          SubstitutionFrameworkAtoms[NumberOfSubstitutionRules][2]);

      // increase the amount of modifcation rules by 1
      NumberOfSubstitutionRules++;
    }

    if((strcasecmp("ModifyOxgensConnectedToAluminium",keyword)==0)||(strcasecmp("ModifyOxygensConnectedToAluminium",keyword)==0))
    {
      // allocate additional memory for this modification rule
      ModificationRuleType=(int*)realloc(ModificationRuleType,(NumberOfModificationRules+1)*sizeof(int));
      ModifyFrameworkAtoms=(char(*)[6][256])realloc(ModifyFrameworkAtoms,(NumberOfModificationRules+1)*sizeof(char[6][256]));
      ModifyFrameworkAtomTypes=(int(*)[6])realloc(ModifyFrameworkAtomTypes,(NumberOfModificationRules+1)*sizeof(int[6]));
      ModificationRuleType[NumberOfModificationRules]=MODIFY_FRAMEWORKATOM_CONNECTED_TO;
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][3],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][4],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][5],"wild card");

      // scan the arguments from the input
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][0],"O");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][1],"Oa");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][2],"Al");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][3],"wild card");

      // increase the amount of modifcation rules by 1
      NumberOfModificationRules++;
    }
    if(strcasecmp("ModifyFrameworkAtomConnectedTo",keyword)==0)
    {
      // allocate additional memory for this modification rule
      ModificationRuleType=(int*)realloc(ModificationRuleType,(NumberOfModificationRules+1)*sizeof(int));
      ModifyFrameworkAtoms=(char(*)[6][256])realloc(ModifyFrameworkAtoms,(NumberOfModificationRules+1)*sizeof(char[6][256]));
      ModifyFrameworkAtomTypes=(int(*)[6])realloc(ModifyFrameworkAtomTypes,(NumberOfModificationRules+1)*sizeof(int[6]));
      ModificationRuleType[NumberOfModificationRules]=MODIFY_FRAMEWORKATOM_CONNECTED_TO;
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][2],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][3],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][4],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][5],"wild card");

      // scan the arguments from the input
      sscanf(arguments,"%s %s %s %s",
          ModifyFrameworkAtoms[NumberOfModificationRules][0],
          ModifyFrameworkAtoms[NumberOfModificationRules][1],
          ModifyFrameworkAtoms[NumberOfModificationRules][2],
          ModifyFrameworkAtoms[NumberOfModificationRules][3]);

      // increase the amount of modifcation rules by 1
      NumberOfModificationRules++;
    }
    if(strcasecmp("ModifyFrameworkDimer",keyword)==0)
    {
      // allocate additional memory for this modification rule
      ModificationRuleType=(int*)realloc(ModificationRuleType,(NumberOfModificationRules+1)*sizeof(int));
      ModifyFrameworkAtoms=(char(*)[6][256])realloc(ModifyFrameworkAtoms,(NumberOfModificationRules+1)*sizeof(char[6][256]));
      ModifyFrameworkAtomTypes=(int(*)[6])realloc(ModifyFrameworkAtomTypes,(NumberOfModificationRules+1)*sizeof(int[6]));
      ModificationRuleType[NumberOfModificationRules]=MODIFY_FRAMEWORKATOM_DIMER;
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][0],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][1],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][2],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][3],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][4],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][5],"wild card");

      // scan the arguments from the input
      sscanf(arguments,"%s %s %s %s",
          ModifyFrameworkAtoms[NumberOfModificationRules][0],
          ModifyFrameworkAtoms[NumberOfModificationRules][1],
          ModifyFrameworkAtoms[NumberOfModificationRules][2],
          ModifyFrameworkAtoms[NumberOfModificationRules][3]);

      // increase the amount of modifcation rules by 1
      NumberOfModificationRules++;
    }
    if(strcasecmp("ModifyFrameworkTriple",keyword)==0)
    {
      // allocate additional memory for this modification rule
      ModificationRuleType=(int*)realloc(ModificationRuleType,(NumberOfModificationRules+1)*sizeof(int));
      ModifyFrameworkAtoms=(char(*)[6][256])realloc(ModifyFrameworkAtoms,(NumberOfModificationRules+1)*sizeof(char[6][256]));
      ModifyFrameworkAtomTypes=(int(*)[6])realloc(ModifyFrameworkAtomTypes,(NumberOfModificationRules+1)*sizeof(int[6]));
      ModificationRuleType[NumberOfModificationRules]=MODIFY_FRAMEWORKATOM_TRIPLE;
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][0],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][1],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][2],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][3],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][4],"wild card");
      strcpy(ModifyFrameworkAtoms[NumberOfModificationRules][5],"wild card");

      // scan the arguments from the input
      sscanf(arguments,"%s %s %s %s %s %s",
          ModifyFrameworkAtoms[NumberOfModificationRules][0],
          ModifyFrameworkAtoms[NumberOfModificationRules][1],
          ModifyFrameworkAtoms[NumberOfModificationRules][2],
          ModifyFrameworkAtoms[NumberOfModificationRules][3],
          ModifyFrameworkAtoms[NumberOfModificationRules][4],
          ModifyFrameworkAtoms[NumberOfModificationRules][5]);

      // increase the amount of modifcation rules by 1
      NumberOfModificationRules++;
    }
    if(strcasecmp("ForbiddenFrameworkConnectivity",keyword)==0)
    {
      // allocate additional memory for this modification rule
      ForbiddenConnectivityAtoms=(char(*)[3][256])realloc(ForbiddenConnectivityAtoms,(NumberOfForbiddenConnectivityRules+1)*sizeof(char[3][256]));
      ForbiddenConnectivityTypes=(int(*)[3])realloc(ForbiddenConnectivityTypes,(NumberOfForbiddenConnectivityRules+1)*sizeof(int[3]));

      // scan the arguments from the input
      sscanf(arguments,"%s %s %s",
          ForbiddenConnectivityAtoms[NumberOfForbiddenConnectivityRules][0],
          ForbiddenConnectivityAtoms[NumberOfForbiddenConnectivityRules][1],
          ForbiddenConnectivityAtoms[NumberOfForbiddenConnectivityRules][2]);

      // increase the amount of modifcation rules by 1
      NumberOfForbiddenConnectivityRules++;
    }
    if(strcasecmp("ConfirmPseudoAtomWithoutVDWInteraction",keyword)==0)
    {
      PseudoAtomsWithoutVDWInteraction=(char(*)[32])realloc(PseudoAtomsWithoutVDWInteraction,
          (NumberOfPseudoAtomsWithoutVDWInteraction+1)*sizeof(char[32]));
      sscanf(arguments,"%s",PseudoAtomsWithoutVDWInteraction[NumberOfPseudoAtomsWithoutVDWInteraction]);
      NumberOfPseudoAtomsWithoutVDWInteraction++;
    }

    if(strcasecmp("RemoveBondNeighboursFromLongRangeInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) RemoveBondNeighboursFromLongRangeInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) RemoveBondNeighboursFromLongRangeInteraction=FALSE;
    }
    if(strcasecmp("RemoveBendNeighboursFromLongRangeInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) RemoveBendNeighboursFromLongRangeInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) RemoveBendNeighboursFromLongRangeInteraction=FALSE;
    }
    if(strcasecmp("RemoveTorsionNeighboursFromLongRangeInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) RemoveTorsionNeighboursFromLongRangeInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) RemoveTorsionNeighboursFromLongRangeInteraction=FALSE;
    }
    if(strcasecmp("Remove12NeighboursFromChargeChargeInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove12NeighboursFromChargeChargeInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove12NeighboursFromChargeChargeInteraction=FALSE;
    }
    if(strcasecmp("Remove13NeighboursFromChargeChargeInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove13NeighboursFromChargeChargeInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove13NeighboursFromChargeChargeInteraction=FALSE;
    }
    if(strcasecmp("Remove14NeighboursFromChargeChargeInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove14NeighboursFromChargeChargeInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove14NeighboursFromChargeChargeInteraction=FALSE;
    }
    if(strcasecmp("Remove11NeighboursFromChargeBondDipoleInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove11NeighboursFromChargeBondDipoleInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove11NeighboursFromChargeBondDipoleInteraction=FALSE;
    }
    if(strcasecmp("Remove12NeighboursFromChargeBondDipoleInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove12NeighboursFromChargeBondDipoleInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove12NeighboursFromChargeBondDipoleInteraction=FALSE;
    }
    if(strcasecmp("Remove13NeighboursFromChargeBondDipoleInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove13NeighboursFromChargeBondDipoleInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove13NeighboursFromChargeBondDipoleInteraction=FALSE;
    }
    if(strcasecmp("Remove14NeighboursFromChargeBondDipoleInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove14NeighboursFromChargeBondDipoleInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove14NeighboursFromChargeBondDipoleInteraction=FALSE;
    }
    if(strcasecmp("Remove12NeighboursFromBondDipoleBondDipoleInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove12NeighboursFromBondDipoleBondDipoleInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove12NeighboursFromBondDipoleBondDipoleInteraction=FALSE;
    }
    if(strcasecmp("Remove13NeighboursFromBondDipoleBondDipoleInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove13NeighboursFromBondDipoleBondDipoleInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove13NeighboursFromBondDipoleBondDipoleInteraction=FALSE;
    }
    if(strcasecmp("Remove14NeighboursFromBondDipoleBondDipoleInteraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Remove14NeighboursFromBondDipoleBondDipoleInteraction=TRUE;
      if(strcasecmp("no",firstargument)==0) Remove14NeighboursFromBondDipoleBondDipoleInteraction=FALSE;
    }


    // read thermo/barostat options
    if(strcasecmp("ExternalTemperature",keyword)==0) 
    {
      sscanf(arguments,"%lf",&therm_baro_stats.ExternalTemperature[CurrentSystem]);
      Beta[CurrentSystem]=1.0/(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]);
    }
    if(strcasecmp("ExternalPressure",keyword)==0)
    {
      arg_pointer=arguments;
      NumberOfIsothermPressures=0;
      while(sscanf(arg_pointer,"%lf%n",&therm_baro_stats.ExternalPressure[CurrentSystem][NumberOfIsothermPressures],&n)==1)
      {
        arg_pointer+=n;
        NumberOfIsothermPressures++;
      }
      for(i=0;i<NumberOfIsothermPressures;i++)
        therm_baro_stats.ExternalPressure[CurrentSystem][i]/=PRESSURE_CONVERSION_FACTOR; // convert from Pascal
    }
    if(strcasecmp("ThermostatChainLength",keyword)==0) sscanf(arguments,"%d",&therm_baro_stats.ThermostatChainLength);
    if(strcasecmp("BarostatChainLength",keyword)==0) sscanf(arguments,"%d",&therm_baro_stats.BarostatChainLength);
    if(strcasecmp("NumberOfYoshidaSuzukiSteps",keyword)==0) sscanf(arguments,"%d",&therm_baro_stats.NumberOfYoshidaSuzukiSteps);
    if(strcasecmp("TimeScaleParameterThermostat",keyword)==0) 
      sscanf(arguments,"%lf",&therm_baro_stats.time_scale_parameter_thermostat);
    if(strcasecmp("TimeScaleParameterBarostat",keyword)==0) 
      sscanf(arguments,"%lf",&therm_baro_stats.time_scale_parameter_barostat);

    // read MD ensembles
    if(strcasecmp("Ensemble",keyword)==0) 
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          if(strcasecmp("NVE",keyword)==0)
            Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NVE;
          if(strcasecmp("NVT",keyword)==0)
            Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NVT;
          if(strcasecmp("NPT",keyword)==0)
            Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPT;
          if(strcasecmp("NPH",keyword)==0)
            Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPH;
          if(strcasecmp("NPTPR",keyword)==0)
            Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPTPR;
          if(strcasecmp("NPHPR",keyword)==0)
            Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPHPR;
          if(strcasecmp("NPTPR_ISO",keyword)==0) {Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPTPR; NPTPRCellType[index]=ISOTROPIC;}
          if(strcasecmp("NPTPR_ANISO",keyword)==0) {Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPTPR; NPTPRCellType[index]=ANISOTROPIC;}
          if(strcasecmp("NPTPR_MONO",keyword)==0) {Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPTPR; NPTPRCellType[index]=MONOCLINIC;}
          if(strcasecmp("NPTPR_UPPER",keyword)==0) {Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPTPR; NPTPRCellType[index]=REGULAR_UPPER_TRIANGLE;}
          if(strcasecmp("NPTPR_MONO_UPPER",keyword)==0) {Ensemble[index]=InitEnsemble[index]=RunEnsemble[index]=NPTPR; NPTPRCellType[index]=MONOCLINIC_UPPER_TRIANGLE;}

          index++;
        }
      }while(index<NumberOfSystems);
    }
    if(strcasecmp("InitEnsemble",keyword)==0)
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          if(strcasecmp("NVE",keyword)==0) InitEnsemble[index]=NVE;
          if(strcasecmp("NVT",keyword)==0) InitEnsemble[index]=NVT;
          if(strcasecmp("NPT",keyword)==0) InitEnsemble[index]=NPT;
          if(strcasecmp("NPH",keyword)==0) InitEnsemble[index]=NPH;
          if(strcasecmp("NPTPR",keyword)==0) InitEnsemble[index]=NPTPR;
          if(strcasecmp("NPHPR",keyword)==0) InitEnsemble[index]=NPHPR;
          if(strcasecmp("NPTPR_ISO",keyword)==0) {InitEnsemble[index]=NPTPR; NPTPRCellType[index]=ISOTROPIC;}
          if(strcasecmp("NPTPR_ANISO",keyword)==0) {InitEnsemble[index]=NPTPR; NPTPRCellType[index]=ANISOTROPIC;}
          if(strcasecmp("NPTPR_MONO",keyword)==0) {InitEnsemble[index]=NPTPR; NPTPRCellType[index]=MONOCLINIC;}
          if(strcasecmp("NPTPR_UPPER",keyword)==0) {InitEnsemble[index]=NPTPR; NPTPRCellType[index]=REGULAR_UPPER_TRIANGLE;}
          if(strcasecmp("NPTPR_MONO_UPPER",keyword)==0) {InitEnsemble[index]=NPTPR; NPTPRCellType[index]=MONOCLINIC_UPPER_TRIANGLE;}
          index++;
        }
      }while(index<NumberOfSystems);
    }
    if(strcasecmp("RunEnsemble",keyword)==0)
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          if(strcasecmp("NVE",keyword)==0) RunEnsemble[index]=NVE;
          if(strcasecmp("NVT",keyword)==0) RunEnsemble[index]=NVT;
          if(strcasecmp("NPT",keyword)==0) RunEnsemble[index]=NPT;
          if(strcasecmp("NPH",keyword)==0) RunEnsemble[index]=NPH;
          if(strcasecmp("NPTPR",keyword)==0) RunEnsemble[index]=NPTPR;
          if(strcasecmp("NPHPR",keyword)==0) RunEnsemble[index]=NPHPR;
          if(strcasecmp("NPTPR_ISO",keyword)==0) {RunEnsemble[index]=NPTPR; NPTPRCellType[index]=ISOTROPIC;}
          if(strcasecmp("NPTPR_ANISO",keyword)==0) {RunEnsemble[index]=NPTPR; NPTPRCellType[index]=ANISOTROPIC;}
          if(strcasecmp("NPTPR_MONO",keyword)==0) {RunEnsemble[index]=NPTPR; NPTPRCellType[index]=MONOCLINIC;}
          if(strcasecmp("NPTPR_UPPER",keyword)==0) {RunEnsemble[index]=NPTPR; NPTPRCellType[index]=REGULAR_UPPER_TRIANGLE;}
          if(strcasecmp("NPTPR_MONO_UPPER",keyword)==0) {RunEnsemble[index]=NPTPR; NPTPRCellType[index]=MONOCLINIC_UPPER_TRIANGLE;}
          index++;
        }
      }while(index<NumberOfSystems);
    } 
    if(strcasecmp("NPTPRCellType",keyword)==0)
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          if(strcasecmp("Regular",keyword)==0) NPTPRCellType[index]=REGULAR;
          if(strcasecmp("Monoclinic",keyword)==0) NPTPRCellType[index]=MONOCLINIC;
          if(strcasecmp("RegularUpperTriangle",keyword)==0) NPTPRCellType[index]=REGULAR_UPPER_TRIANGLE;
          if(strcasecmp("MonoclinicUpperTriangle",keyword)==0) NPTPRCellType[index]=MONOCLINIC_UPPER_TRIANGLE;
          if(strcasecmp("Isotropic",keyword)==0) NPTPRCellType[index]=ISOTROPIC;
          if(strcasecmp("Anisotropic",keyword)==0) NPTPRCellType[index]=ANISOTROPIC;
          index++;
        }
      }while(index<NumberOfSystems);
    }

    if(strcasecmp("MonoclinicAngleType",keyword)==0)
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          if(strcasecmp("Alpha",keyword)==0) MonoclinicAngleType[index]=MONOCLINIC_ALPHA_ANGLE;
          if(strcasecmp("Beta",keyword)==0) MonoclinicAngleType[index]=MONOCLINIC_BETA_ANGLE;
          if(strcasecmp("Gamma",keyword)==0) MonoclinicAngleType[index]=MONOCLINIC_GAMMA_ANGLE;
          index++;
        }
      }while(index<NumberOfSystems);
    }


    // read box-sizes and angles
    if(strcasecmp("Box",keyword)==0)
    {
      if(sscanf(arguments,"%d",&CurrentSystem))
      {
        if(CurrentSystem>=NumberOfSystems)
        {
          printf("Input error in 'Box'-command: the number of the system is incorrect (%d given, 0 <= system-number < %d)\n",CurrentSystem,NumberOfSystems);
          exit(0);
        }

        if(strlen(Framework[CurrentSystem].Name[0])>0)
        {
          printf("Input error in 'Box'-command: the number of the system %d is already taken\n",CurrentSystem);
          exit(0);
        }

        sprintf(Framework[CurrentSystem].Name[0],"Box");
        InitializeBox[CurrentSystem]=TRUE;
      }
    }
    if(strcasecmp("BoxMatrix",keyword)==0)
    {
      if(sscanf(arguments,"%d",&CurrentSystem))
      {
        fgets(arguments,1024,FilePtr);
        sscanf(arguments,"%lf %lf %lf\n",&Box[CurrentSystem].ax,&Box[CurrentSystem].bx,&Box[CurrentSystem].cx);
        fgets(arguments,1024,FilePtr);
        sscanf(arguments,"%lf %lf %lf\n",&Box[CurrentSystem].ay,&Box[CurrentSystem].by,&Box[CurrentSystem].cy);
        fgets(arguments,1024,FilePtr);
        sscanf(arguments,"%lf %lf %lf\n",&Box[CurrentSystem].az,&Box[CurrentSystem].bz,&Box[CurrentSystem].cz);

        NumberOfUnitCells[CurrentSystem].x=1;
        NumberOfUnitCells[CurrentSystem].y=1;
        NumberOfUnitCells[CurrentSystem].z=1;

        UnitCellBox[CurrentSystem]=Box[CurrentSystem];
        CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
        InverseUnitCellBox[CurrentSystem]=InverseBox[CurrentSystem];
        InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);
        CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);

        AlphaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bx);
        BetaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].by);
        GammaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bz);

        InitializeBox[CurrentSystem]=FALSE;

        sprintf(Framework[CurrentSystem].Name[CurrentFramework],"Box");
      }
    }
    if((strcasecmp("BoxLengths",keyword)==0)||(strcasecmp("CellLengths",keyword)==0))
    {
      UnitCellSize[CurrentSystem].x=UnitCellSize[CurrentSystem].y=UnitCellSize[CurrentSystem].z=0.0;
      switch(Dimension)
      {
        case 2:
          sscanf(arguments,"%lf %lf",&UnitCellSize[CurrentSystem].x,&UnitCellSize[CurrentSystem].y);
          break;
        case 3:
          sscanf(arguments,"%lf %lf %lf",&UnitCellSize[CurrentSystem].x,&UnitCellSize[CurrentSystem].y,&UnitCellSize[CurrentSystem].z);
          break;
      }
      InitializeBox[CurrentSystem]=TRUE;
    }
    if((strcasecmp("BoxAngles",keyword)==0)||(strcasecmp("CellAngles",keyword)==0))
    {
      AlphaAngle[CurrentSystem]=BetaAngle[CurrentSystem]=GammaAngle[CurrentSystem]=0.0;
      switch(Dimension)
      {
        case 2:
          sscanf(arguments,"%lf %lf",&AlphaAngle[CurrentSystem],&BetaAngle[CurrentSystem]);
          break;
        case 3:
          sscanf(arguments,"%lf %lf %lf",&AlphaAngle[CurrentSystem],&BetaAngle[CurrentSystem],&GammaAngle[CurrentSystem]);
          break;
      }
      InitializeBox[CurrentSystem]=TRUE;
    }



    // read framework options
    if(strcasecmp("Framework",keyword)==0) 
    {
      if(sscanf(arguments,"%d",&CurrentSystem))
      {
        if(strlen(Framework[CurrentSystem].Name[0])>0)
        {
          printf("Input error in 'Framework'-command: the number of the system %d is already taken\n",CurrentSystem);
          exit(0);
        }

        sprintf(Framework[CurrentSystem].Name[0],"Box");
        Framework[CurrentSystem].FrameworkModel=FULL;
        read_frameworks[CurrentSystem]=TRUE;
        CurrentFramework=-1;
        
      }
    }
    if(strcasecmp("FrameworkName",keyword)==0)
    {
      CurrentFramework++;
      sscanf(arguments,"%s",Framework[CurrentSystem].Name[CurrentFramework]);
    }
    if(strcasecmp("FrameworkDefinitions",keyword)==0) sscanf(arguments,"%s",Framework[CurrentSystem].FrameworkDefinitions);
    if(strcasecmp("HeliumVoidFraction",keyword)==0) sscanf(arguments,"%lf",&HeliumVoidFraction[CurrentSystem]);
    if(strcasecmp("ExcessVolume",keyword)==0) sscanf(arguments,"%lf",&ExcessVolume[CurrentSystem]);
    if(strcasecmp("IonsName",keyword)==0) sscanf(arguments,"%s",Framework[CurrentSystem].NameIons);

    if(strcasecmp("ContinuousFraction",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) UseContinuousFraction=TRUE;
      if(strcasecmp("no",firstargument)==0) UseContinuousFraction=FALSE;
    }

    if(strcasecmp("CalculateSpaceGroup",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].CalculateSpaceGroup[CurrentFramework]=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].CalculateSpaceGroup[CurrentFramework]=FALSE;
    }
    if(strcasecmp("RemoveAtomNumberCodeFromLabel",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].RemoveAtomNumberCodeFromLabel[CurrentFramework]=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].RemoveAtomNumberCodeFromLabel[CurrentFramework]=FALSE;
    }
    if(strcasecmp("AddAtomNumberCodeToLabel",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].AddAtomNumberCodeToLabel[CurrentFramework]=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].AddAtomNumberCodeToLabel[CurrentFramework]=FALSE;
    }
    if(strcasecmp("ForceSpaceGroupDetection",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].ForceSpaceGroupDetection=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].ForceSpaceGroupDetection=FALSE;
    }
    if(strcasecmp("ReadCIFAsCartesian",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].ReadCIFAsCartesian=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].ReadCIFAsCartesian=FALSE;
    }
    if(strcasecmp("RestrictFrameworkAtomsToBox",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].RestrictFrameworkAtomsToBox=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].RestrictFrameworkAtomsToBox=FALSE;
    }
    if(strcasecmp("InputFileType",keyword)==0)
    {
      if(strcasecmp("xyz",firstargument)==0) Framework[CurrentSystem].InputFileType[CurrentFramework]=XYZ_FORMAT;
      if(strcasecmp("mol",firstargument)==0) Framework[CurrentSystem].InputFileType[CurrentFramework]=MOL_FORMAT;
      if(strcasecmp("cssr",firstargument)==0) Framework[CurrentSystem].InputFileType[CurrentFramework]=CSSR_FORMAT;
      if(strcasecmp("pdb",firstargument)==0) Framework[CurrentSystem].InputFileType[CurrentFramework]=PDB_FORMAT;
      if(strcasecmp("cif",firstargument)==0) Framework[CurrentSystem].InputFileType[CurrentFramework]=CIF_FORMAT;
    }
    if(strcasecmp("UnitCells",keyword)==0) sscanf(arguments,"%d %d %d\n",&NumberOfUnitCells[CurrentSystem].x,
         &NumberOfUnitCells[CurrentSystem].y,&NumberOfUnitCells[CurrentSystem].z);
    if(strcasecmp("ReplicaUnitCells",keyword)==0) sscanf(arguments,"%d %d %d\n",&NumberOfReplicaCells[CurrentSystem].x,
         &NumberOfReplicaCells[CurrentSystem].y,&NumberOfReplicaCells[CurrentSystem].z);
    if(strcasecmp("ShiftUnitCells",keyword)==0) sscanf(arguments,"%lf %lf %lf",
         &Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].x,
         &Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].y,
         &Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].z);
    if(strcasecmp("FlexibleFramework",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) 
      {
        Framework[CurrentSystem].FrameworkModels[CurrentFramework]=FLEXIBLE;
        Framework[CurrentSystem].FrameworkModel=FLEXIBLE;
      }
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].FrameworkModels[CurrentFramework]=RIGID;
    }


    // read component and component options (MC-moves, starting bead, molfraction, etc)
    if(strcasecmp("Component",keyword)==0)
    {
      if(sscanf(arguments,"%d %*s %s",&CurrentComponent,string))
      {
        if(CurrentComponent>=NumberOfComponents)
        {
          printf("Input error in 'Component'-command: the number of the component is incorrect (%d given, 0 <= component-number < %d)\n",CurrentComponent,NumberOfComponents);
          exit(0);
        }

        if(strlen(Components[CurrentComponent].Name)>0)
        {
          printf("Input error in 'Component'-command: the number of the component %d is already taken\n",CurrentComponent);
          exit(0);
        }

        strcpy(Components[CurrentComponent].Name,string);
      }
    }
    if(strcasecmp("MoleculeDefinition",keyword)==0) sscanf(arguments,"%s",Components[CurrentComponent].MoleculeDefinition);
    if(strcasecmp("StartingBead",keyword)==0) sscanf(arguments,"%d",&Components[CurrentComponent].StartingBead);
    if(strcasecmp("Intra14VDWScalingValue",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].Intra14VDWScalingValue);
    if(strcasecmp("Intra14ChargeChargeScalingValue",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].Intra14ChargeChargeScalingValue);
    if(strcasecmp("BlockPockets",keyword)==0)
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          if(strcasecmp("yes",keyword)==0) Components[CurrentComponent].BlockPockets[index]=TRUE;
          if(strcasecmp("no",keyword)==0) Components[CurrentComponent].BlockPockets[index]=FALSE;
          index++;
        }
      }while(index<NumberOfSystems);
    }
    if(strcasecmp("BlockPocketsFileName",keyword)==0)
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          strcpy(Components[CurrentComponent].BlockPocketsFilename[index],keyword);
          index++;
        }
      }while(index<NumberOfSystems);
    }
    if(strcasecmp("ComputeFreeEnergyProfile",keyword)==0)
    {
      index=0;
      do
      {
        if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
        {
          if(strcasecmp("yes",keyword)==0) Components[CurrentComponent].ComputeFreeEnergyProfile[index]=ComputeFreeEnergyProfile[index]=TRUE;
          if(strcasecmp("no",keyword)==0) Components[CurrentComponent].ComputeFreeEnergyProfile[index]=FALSE;
          index++;
        }
      }while(index<NumberOfSystems);
    }
    if(strcasecmp("MolFraction",keyword)==0)
    {
      index=0;
      do
      {
        if(index==NumberOfSystems)
        {
          printf("Reading mol fraction for more systems than the maximum allowed: %d\n",NumberOfSystems);
          exit(0);
        }
      }while(sscanf(arguments,"%lf %[^\n]",&Components[CurrentComponent].MolFraction[index++],arguments)>1);
    }
    if(strcasecmp("FrameworkChangeMoveProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityFrameworkChangeMove);
    if(strcasecmp("FrameworkShiftMoveProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityFrameworkShiftMove);
    if(strcasecmp("VolumeChangeProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityVolumeChangeMove);
    if(strcasecmp("ParallelTemperingProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityParallelTemperingMove);
    if(strcasecmp("HyperParallelTemperingProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityHyperParallelTemperingMove);
    if(strcasecmp("ParallelMolFractionProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityParallelMolFractionMove);
    if(strcasecmp("ParallelMolFractionComponentA",keyword)==0) 
    {
      sscanf(arguments,"%d",&ParallelMolFractionComponentA);
      if((ParallelMolFractionComponentA<0)||(ParallelMolFractionComponentA>=NumberOfComponents))
      {
        printf("Input error in 'ParallelMolFractionComponentA'-command: the number of the component %d is incorrect (should be < %d)\n",ParallelMolFractionComponentA,NumberOfComponents);
        exit(0);
      }
    }
    if(strcasecmp("ParallelMolFractionComponentB",keyword)==0) 
    {
      sscanf(arguments,"%d",&ParallelMolFractionComponentB);
      if((ParallelMolFractionComponentB<0)||(ParallelMolFractionComponentB>=NumberOfComponents))
      {
        printf("Input error in 'ParallelMolFractionComponentB'-command: the number of the component %d is incorrect (should be < %d)\n",ParallelMolFractionComponentB,NumberOfComponents);
        exit(0);
      }
    }
    if(strcasecmp("ChiralInversionProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityChiralInversionMove);
    if(strcasecmp("MaximumVolumeChange",keyword)==0) 
      sscanf(arguments,"system %d: %lf",&CurrentSystem,&MaximumVolumeChange[CurrentSystem]);
    if(strcasecmp("BoxShapeChangeProbability",keyword)==0) 
    {
       sscanf(arguments,"%lf",&ProbabilityBoxShapeChangeMove);
       if(BoundaryCondition[CurrentSystem]==UNINITIALIZED_BOUNDARY_CONDITION)
         BoundaryCondition[CurrentSystem]=TRICLINIC;
    }
    if(strcasecmp("BoundaryCondition",keyword)==0)
    {
      if(strcasecmp("None",firstargument)==0) BoundaryCondition[CurrentSystem]=FINITE;
      if(strcasecmp("Reactangular",firstargument)==0) BoundaryCondition[CurrentSystem]=RECTANGULAR;
      if(strcasecmp("Triclinic",firstargument)==0) BoundaryCondition[CurrentSystem]=TRICLINIC;
    }
    if(strcasecmp("GibbsVolumeChangeProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityGibbsVolumeChangeMove);
    if(strcasecmp("GibbsSwapProbability",keyword)==0) 
    {
      sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityGibbsSwapChangeMove);
      if(Components[CurrentComponent].ProbabilityGibbsSwapChangeMove>0.0)
        Swapable=TRUE;
    }
    if(strcasecmp("GibbsIdentityChangeProbability",keyword)==0) 
    {
      sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove);
      if(Components[CurrentComponent].ProbabilityGibbsIdentityChangeMove>0.0)
        Swapable=TRUE;
    }
    if(strcasecmp("TranslationProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityTranslationMove);
    if(strcasecmp("TranslationDirection",keyword)==0)
    {
      if(strcasecmp("XYZ",firstargument)==0) Components[CurrentComponent].TranslationDirection=XYZ_DIR;
      if(strcasecmp("XY",firstargument)==0) Components[CurrentComponent].TranslationDirection=XY_DIR;
      if(strcasecmp("XZ",firstargument)==0) Components[CurrentComponent].TranslationDirection=XZ_DIR;
      if(strcasecmp("YZ",firstargument)==0) Components[CurrentComponent].TranslationDirection=YZ_DIR;
      if(strcasecmp("X",firstargument)==0) Components[CurrentComponent].TranslationDirection=X_DIR;
      if(strcasecmp("Y",firstargument)==0) Components[CurrentComponent].TranslationDirection=Y_DIR;
      if(strcasecmp("Z",firstargument)==0) Components[CurrentComponent].TranslationDirection=Z_DIR;

      if(strcasecmp("ABC",firstargument)==0) Components[CurrentComponent].TranslationDirection=ABC_DIR;
      if(strcasecmp("AB",firstargument)==0) Components[CurrentComponent].TranslationDirection=AB_DIR;
      if(strcasecmp("AC",firstargument)==0) Components[CurrentComponent].TranslationDirection=AC_DIR;
      if(strcasecmp("BC",firstargument)==0) Components[CurrentComponent].TranslationDirection=BC_DIR;

      if(strcasecmp("A",firstargument)==0) Components[CurrentComponent].TranslationDirection=A_DIR;
      if(strcasecmp("B",firstargument)==0) Components[CurrentComponent].TranslationDirection=B_DIR;
      if(strcasecmp("C",firstargument)==0) Components[CurrentComponent].TranslationDirection=C_DIR;

      if(strcasecmp("ORTHOGONAL_TO_AB_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_AB_DIR;
      if(strcasecmp("ORTHOGONAL_TO_AC_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_AC_DIR;
      if(strcasecmp("ORTHOGONAL_TO_BC_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_BC_DIR;

      if(strcasecmp("ORTHOGONAL_TO_O_AB_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_O_AB_DIR;
      if(strcasecmp("ORTHOGONAL_TO_O_AC_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_O_AC_DIR;
      if(strcasecmp("ORTHOGONAL_TO_O_BC_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_O_BC_DIR;

      if(strcasecmp("ORTHOGONAL_TO_A_BC_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_A_BC_DIR;
      if(strcasecmp("ORTHOGONAL_TO_B_AC_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_B_AC_DIR;
      if(strcasecmp("ORTHOGONAL_TO_C_AB_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_C_AB_DIR;
      if(strcasecmp("ORTHOGONAL_TO_O_ABC_DIR",firstargument)==0) Components[CurrentComponent].TranslationDirection=ORTHOGONAL_TO_O_ABC_DIR;
    }
    if(strcasecmp("RandomTranslationProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityRandomTranslationMove);
    if(strcasecmp("RotationProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityRotationMove);
    if(strcasecmp("CBMCProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityPartialReinsertionMove);
    if(strcasecmp("PartialReinsertionProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityPartialReinsertionMove);
    if(strcasecmp("ReinsertionProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityReinsertionMove);
    if(strcasecmp("RegrowProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityReinsertionMove);
    if(strcasecmp("ReinsertionInPlaceProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityReinsertionInPlaceMove);
    if(strcasecmp("RegrowInPlaceProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityReinsertionInPlaceMove);
    if(strcasecmp("HybridMCMDMoveProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityHybridNVEMove);
    if(strcasecmp("HybridNVEMoveProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityHybridNVEMove);
    if(strcasecmp("HybridNPHMoveProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityHybridNPHMove);
    if(strcasecmp("HybridNPHPRMoveProbability",keyword)==0) sscanf(arguments,"%lf",&ProbabilityHybridNPHPRMove);
    if(strcasecmp("NumberOfHybridNVESteps",keyword)==0) sscanf(arguments,"%d",&NumberOfHybridNVESteps);
    if(strcasecmp("NumberOfHybridNPHSteps",keyword)==0) sscanf(arguments,"%d",&NumberOfHybridNPHSteps);
    if(strcasecmp("NumberOfHybridNPHPRSteps",keyword)==0) sscanf(arguments,"%d",&NumberOfHybridNPHPRSteps);
    if(strcasecmp("ReinsertionInPlaneProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityReinsertionInPlaneMove);
    if(strcasecmp("IdentityChangeProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityIdentityChangeMove);
    if(strcasecmp("NumberOfIdentityChanges",keyword)==0) sscanf(arguments,"%d",&Components[CurrentComponent].NumberOfIdentityChanges);
    if(strcasecmp("IdentityChangesList",keyword)==0) 
      for(j=0;j<Components[CurrentComponent].NumberOfIdentityChanges;j++)
        sscanf(arguments,"%d %[^\n]",&Components[CurrentComponent].IdentityChanges[j],arguments);
    if(strcasecmp("NumberOfGibbsIdentityChanges",keyword)==0) sscanf(arguments,"%d",&Components[CurrentComponent].NumberOfGibbsIdentityChanges);
    if(strcasecmp("GibbsIdentityChangesList",keyword)==0) 
      for(j=0;j<Components[CurrentComponent].NumberOfGibbsIdentityChanges;j++)
        sscanf(arguments,"%d %[^\n]",&Components[CurrentComponent].GibbsIdentityChanges[j],arguments);
    if(strcasecmp("SwapProbability",keyword)==0) 
      if(sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilitySwapMove))
        if(Components[CurrentComponent].ProbabilitySwapMove>0.0)
        {
          Components[CurrentComponent].Swapable=TRUE;
          Swapable=TRUE;
        }
    if(strcasecmp("WidomProbability",keyword)==0) 
      if(sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilityWidomMove))
        if(Components[CurrentComponent].ProbabilityWidomMove>0.0)
          Components[CurrentComponent].Widom=TRUE;
    if(strcasecmp("WidomParticleInsertionComponent",keyword)==0) sscanf(arguments,"%d",&WidomParticleInsertionComponent);
    if(strcasecmp("SurfaceAreaProbability",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].ProbabilitySurfaceAreaMove);
    if(strcasecmp("RestrictEnantionface",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) Components[CurrentComponent].RestrictEnantionface=TRUE;
      if(strcasecmp("no",firstargument)==0) Components[CurrentComponent].RestrictEnantionface=FALSE;
    }
    if(strcasecmp("Enantioface",keyword)==0) 
    {
      if(strcasecmp("Re",firstargument)==0) Components[CurrentComponent].Enantioface=ENANTIOFACE_RE;
      if(strcasecmp("Si",firstargument)==0) Components[CurrentComponent].Enantioface=ENANTIOFACE_SI;
    }
    if(strcasecmp("EnantiofaceAtoms",keyword)==0) 
    {
      sscanf(arguments," %c%d%d %c%d%d %c%d%d %c%d%d %c%d%d %[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
            &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,&charinput5,&intinput9,&intinput10,arguments);

      if((charinput1=='F')||(charinput1=='f'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[0][0]=CATION;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[0][1]=intinput1;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[1][0]=CATION;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[1][1]=intinput3;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[2][0]=CATION;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[2][1]=intinput5;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[3][0]=CATION;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[3][1]=intinput7;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[3][2]=intinput8;

      if((charinput5=='F')||(charinput5=='f'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[4][0]=FRAMEWORK;
      else if((charinput5=='A')||(charinput5=='a'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[4][0]=ADSORBATE;
      else if((charinput5=='C')||(charinput5=='c'))
        Components[CurrentComponent].EnantiofaceAtomDefinitions[4][0]=CATION;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[4][1]=intinput9;
      Components[CurrentComponent].EnantiofaceAtomDefinitions[4][2]=intinput10;
    }

    if(strcasecmp("ExtraFrameworkMolecule",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) Components[CurrentComponent].ExtraFrameworkMolecule=TRUE;
      if(strcasecmp("no",firstargument)==0) Components[CurrentComponent].ExtraFrameworkMolecule=FALSE;
    }
    if(strcasecmp("CreateNumberOfMolecules",keyword)==0)
    {
      index=0;
      do
      {
        if(index==NumberOfSystems)
        {
          printf("Creating molecules for more systems than the maximum allowed: %d\n",NumberOfSystems);
          exit(0);
        }
      }while(sscanf(arguments,"%d %[^\n]",&Components[CurrentComponent].CreateNumberOfMolecules[index++],arguments)>1);
    }
    if(strcasecmp("FugacityCoefficient",keyword)==0)
    {
      index=0;
      do
      { 
        if(index==NumberOfSystems)
        {
          printf("Reading fugacity coefficients (%d) for more systems than the maximum allowed: %d\n",index,NumberOfSystems);
          exit(0);
        }
      }while(sscanf(arguments,"%lf %[^\n]",&Components[CurrentComponent].FugacityCoefficient[index++],arguments)>1);
      for(i=0;i<index;i++)
        if(Components[CurrentComponent].FugacityCoefficient[i]>=0.0) ComputeFugacityCoefficient[i][CurrentComponent]=FALSE;
    }
    if(strcasecmp("IdealGasRosenbluthWeight",keyword)==0)
    {
      index=0;
      do
      {
        if(index==NumberOfSystems)
        {
          printf("Reading ideal Rosenbluth weights for more systems than the maximum allowed: %d\n",NumberOfSystems);
          exit(0);
        }
      }while(sscanf(arguments,"%lf %[^\n]",&Components[CurrentComponent].IdealGasRosenbluthWeight[index++],arguments)>1);
    }
    if(strcasecmp("RestrictMovesToBox",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) 
      {
        Components[CurrentComponent].RestrictMovesToBox=TRUE;
        Components[CurrentComponent].RestrictMoves=TRUE;
      }
      if(strcasecmp("no",firstargument)==0) Components[CurrentComponent].RestrictMovesToBox=FALSE;
    }
    if(strcasecmp("BoxAxisABC_Min",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Min.x,
             &Components[CurrentComponent].BoxAxisABC_Min.y,&Components[CurrentComponent].BoxAxisABC_Min.z);
    if(strcasecmp("BoxAxisABC_Max",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Max.x,
             &Components[CurrentComponent].BoxAxisABC_Max.y,&Components[CurrentComponent].BoxAxisABC_Max.z);
    if(strcasecmp("BoxAxisABC_Min2",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Min2.x,
             &Components[CurrentComponent].BoxAxisABC_Min2.y,&Components[CurrentComponent].BoxAxisABC_Min2.z);
    if(strcasecmp("BoxAxisABC_Max2",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Max2.x,
             &Components[CurrentComponent].BoxAxisABC_Max2.y,&Components[CurrentComponent].BoxAxisABC_Max2.z);
    if(strcasecmp("BoxAxisABC_Min3",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Min3.x,
             &Components[CurrentComponent].BoxAxisABC_Min3.y,&Components[CurrentComponent].BoxAxisABC_Min3.z);
    if(strcasecmp("BoxAxisABC_Max3",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Max3.x,
             &Components[CurrentComponent].BoxAxisABC_Max3.y,&Components[CurrentComponent].BoxAxisABC_Max3.z);
    if(strcasecmp("BoxAxisABC_Min4",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Min4.x,
             &Components[CurrentComponent].BoxAxisABC_Min4.y,&Components[CurrentComponent].BoxAxisABC_Min4.z);
    if(strcasecmp("BoxAxisABC_Max4",keyword)==0) 
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].BoxAxisABC_Max4.x,
             &Components[CurrentComponent].BoxAxisABC_Max4.y,&Components[CurrentComponent].BoxAxisABC_Max4.z);

    if(strcasecmp("RestrictMovesToPrism",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0)
      {
        Components[CurrentComponent].RestrictMoves=TRUE;
        Components[CurrentComponent].RestrictMovesToPrisms=TRUE;
      }
      if(strcasecmp("no",firstargument)==0) Components[CurrentComponent].RestrictMovesToPrisms=FALSE;
    }
    if(strcasecmp("Prism",keyword)==0)
    {
      sscanf(arguments,"%d",&CurrentPrism);
      Components[CurrentComponent].RestrictMovesToPrism[CurrentPrism]=TRUE;
    }
    if(strcasecmp("PrismABC_Min",keyword)==0)
    {
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].RestrictPrismABC_Min[CurrentPrism].x,
             &Components[CurrentComponent].RestrictPrismABC_Min[CurrentPrism].y,&Components[CurrentComponent].RestrictPrismABC_Min[CurrentPrism].z);
    }
    if(strcasecmp("PrismABC_Max",keyword)==0)
    {
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].RestrictPrismABC_Max[CurrentPrism].x,
             &Components[CurrentComponent].RestrictPrismABC_Max[CurrentPrism].y,&Components[CurrentComponent].RestrictPrismABC_Max[CurrentPrism].z);
    }


    if(strcasecmp("RestrictMovesToCylinder",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) 
      {
        Components[CurrentComponent].RestrictMoves=TRUE;
        Components[CurrentComponent].RestrictMovesToCylinders=TRUE;
      }
      if(strcasecmp("no",firstargument)==0) Components[CurrentComponent].RestrictMovesToCylinders=FALSE;
    }
    if(strcasecmp("Cylinder",keyword)==0) 
    {
      sscanf(arguments,"%d",&CurrentCylinder);
      Components[CurrentComponent].RestrictMovesToCylinder[CurrentCylinder]=TRUE;
    }
    if(strcasecmp("CylinderABC_Min",keyword)==0)
    {
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].RestrictCylinderABC_Min[CurrentCylinder].x,
             &Components[CurrentComponent].RestrictCylinderABC_Min[CurrentCylinder].y,&Components[CurrentComponent].RestrictCylinderABC_Min[CurrentCylinder].z);
    }
    if(strcasecmp("CylinderABC_Max",keyword)==0)
    {
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].RestrictCylinderABC_Max[CurrentCylinder].x,
             &Components[CurrentComponent].RestrictCylinderABC_Max[CurrentCylinder].y,&Components[CurrentComponent].RestrictCylinderABC_Max[CurrentCylinder].z);
    }
    if(strcasecmp("CylinderCenter",keyword)==0)
    {
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].RestrictCylinderCenter[CurrentCylinder].x,
             &Components[CurrentComponent].RestrictCylinderCenter[CurrentCylinder].y,&Components[CurrentComponent].RestrictCylinderCenter[CurrentCylinder].z);
    }
    if(strcasecmp("CylinderRadius",keyword)==0)
      sscanf(arguments,"%lf",&Components[CurrentComponent].RestrictCylinderRadius[CurrentCylinder]);
    if(strcasecmp("CylinderDiameter",keyword)==0)
    {
      sscanf(arguments,"%lf",&Components[CurrentComponent].RestrictCylinderRadius[CurrentCylinder]);
      Components[CurrentComponent].RestrictCylinderRadius[CurrentCylinder]*=0.5;
    }
    if(strcasecmp("CylinderDirection",keyword)==0)
    {
      if(strcasecmp("X",firstargument)==0) Components[CurrentComponent].RestrictCylinderDirection[CurrentCylinder]=X_DIR;
      if(strcasecmp("Y",firstargument)==0) Components[CurrentComponent].RestrictCylinderDirection[CurrentCylinder]=Y_DIR;
      if(strcasecmp("Z",firstargument)==0) Components[CurrentComponent].RestrictCylinderDirection[CurrentCylinder]=Z_DIR;
    }

    if(strcasecmp("RestrictMovesToSphere",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) 
      {
        Components[CurrentComponent].RestrictMoves=TRUE;
        Components[CurrentComponent].RestrictMovesToSpheres=TRUE;
      }
      if(strcasecmp("no",firstargument)==0) Components[CurrentComponent].RestrictMovesToSpheres=FALSE;
    }
    if(strcasecmp("Sphere",keyword)==0) 
    {
      sscanf(arguments,"%d",&CurrentSphere);
      Components[CurrentComponent].RestrictMovesToSphere[CurrentSphere]=TRUE;
    }
    if(strcasecmp("SphereCenter",keyword)==0)
    {
      sscanf(arguments,"%lf %lf %lf",&Components[CurrentComponent].RestrictSphereCenter[CurrentSphere].x,
             &Components[CurrentComponent].RestrictSphereCenter[CurrentSphere].y,&Components[CurrentComponent].RestrictSphereCenter[CurrentSphere].z);
    }
    if(strcasecmp("SphereRadius",keyword)==0)
      sscanf(arguments,"%lf",&Components[CurrentComponent].RestrictSphereRadius[CurrentSphere]);
    if(strcasecmp("SphereDiameter",keyword)==0)
    {
      sscanf(arguments,"%lf",&Components[CurrentComponent].RestrictSphereRadius[CurrentSphere]);
      Components[CurrentComponent].RestrictSphereRadius[CurrentSphere]*=0.5;
    }

    if(strcasecmp("TransmissionCoefficients",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) Components[CurrentComponent].TransmissionCoefficient=TRUE;
      if(strcasecmp("no",firstargument)==0) Components[CurrentComponent].TransmissionCoefficient=FALSE;
    }
    if(strcasecmp("BiasingMethod",keyword)==0) 
    {
      if(strcasecmp("umbrella",firstargument)==0) Components[CurrentComponent].Biased=UMBRELLA;
      if(strcasecmp("ruizmontero",firstargument)==0) Components[CurrentComponent].Biased=RUIZ_MONTERO;
    }
    if(strcasecmp("RuizMonterofactor",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].RuizMonteroFactor);
    if(strcasecmp("UmbrellaFactor",keyword)==0) sscanf(arguments,"%lf",&Components[CurrentComponent].UmbrellaFactor);
    if(strcasecmp("BiasingProfile",keyword)==0) 
    { 
      Components[CurrentComponent].ReadBiasingFunction=TRUE;
      sscanf(arguments,"%s",Components[CurrentComponent].BiasingFunctionName);
    }
    if(strcasecmp("BiasingDirection",keyword)==0)
    {
      if(strcasecmp("A",firstargument)==0) Components[CurrentComponent].BiasingDirection=A_MAPPING;
      if(strcasecmp("B",firstargument)==0) Components[CurrentComponent].BiasingDirection=B_MAPPING;
      if(strcasecmp("C",firstargument)==0) Components[CurrentComponent].BiasingDirection=C_MAPPING;

      if(strcasecmp("AB_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_AB_DIAGONAL;
      if(strcasecmp("AC_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_AC_DIAGONAL;
      if(strcasecmp("BC_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_BC_DIAGONAL;

      if(strcasecmp("O_AB_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_O_AB_DIAGONAL;
      if(strcasecmp("O_AC_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_O_AC_DIAGONAL;
      if(strcasecmp("O_BC_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_O_BC_DIAGONAL;

      if(strcasecmp("A_BC_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_A_BC_DIAGONAL;
      if(strcasecmp("B_AC_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_B_AC_DIAGONAL;
      if(strcasecmp("C_AB_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_C_AB_DIAGONAL;
      if(strcasecmp("O_ABC_DIAGONAL",firstargument)==0) Components[CurrentComponent].BiasingDirection=MAP_O_ABC_DIAGONAL;
    }
    if(strcasecmp("Histogram",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputePositionHistogram[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputePositionHistogram[CurrentSystem]=FALSE;
    }
    if(strcasecmp("PositionHistogramMappingType",keyword)==0)
    {
      if(strcasecmp("A",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=A_MAPPING;
      if(strcasecmp("B",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=B_MAPPING;
      if(strcasecmp("C",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=C_MAPPING;
      if(strcasecmp("ABC",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=ABC_MAPPING;

      if(strcasecmp("AB_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_AB_DIAGONAL;
      if(strcasecmp("AC_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_AC_DIAGONAL;
      if(strcasecmp("BC_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_BC_DIAGONAL;

      if(strcasecmp("O_AB_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_O_AB_DIAGONAL;
      if(strcasecmp("O_AC_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_O_AC_DIAGONAL;
      if(strcasecmp("O_BC_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_O_BC_DIAGONAL;

      if(strcasecmp("A_BC_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_A_BC_DIAGONAL;
      if(strcasecmp("B_AC_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_B_AC_DIAGONAL;
      if(strcasecmp("C_AB_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_C_AB_DIAGONAL;
      if(strcasecmp("O_ABC_DIAGONAL",firstargument)==0) PositionHistogramMappingType[CurrentSystem]=MAP_O_ABC_DIAGONAL;
    }


    // reading options to measure properties
    // ==================================================================================================

    if(strcasecmp("Movies",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Movies[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) Movies[CurrentSystem]=FALSE;
    }
    if(strcasecmp("MovieScale",keyword)==0) sscanf(arguments,"%lf",&MovieScale);
    if(strcasecmp("WriteMoviesEvery",keyword)==0) sscanf(arguments,"%d",&WriteMoviesEvery[CurrentSystem]);


    // sampling the radial distribution function (RDF)
    if(strcasecmp("ComputeRDF",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeRDF[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeRDF[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteRDFEvery",keyword)==0) sscanf(arguments,"%d",&WriteRDFEvery[CurrentSystem]);
    if(strcasecmp("RDFHistogramSize",keyword)==0) sscanf(arguments,"%d",&RDFHistogramSize[CurrentSystem]);
    if(strcasecmp("RDFRange",keyword)==0) sscanf(arguments,"%lf",&RDFRange[CurrentSystem]);

    // sampling the number-of-molecules histogram
    if(strcasecmp("ComputeNumberOfMoleculesHistogram",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeNumberOfMoleculesHistogram[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeNumberOfMoleculesHistogram[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteNumberOfMoleculesHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteNumberOfMoleculesHistogramEvery[CurrentSystem]);
    if(strcasecmp("NumberOfMoleculesHistogramSize",keyword)==0) sscanf(arguments,"%d",&NumberOfMoleculesHistogramSize[CurrentSystem]);
    if(strcasecmp("NumberOfMoleculesRange",keyword)==0) sscanf(arguments,"%lf",&NumberOfMoleculesRange[CurrentSystem]);

    // sampling position histograms/free energies
    if(strcasecmp("ComputePositionHistogram",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputePositionHistogram[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputePositionHistogram[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WritePositionHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WritePositionHistogramEvery[CurrentSystem]);
    if(strcasecmp("PositionHistogramSize",keyword)==0) sscanf(arguments,"%d",&PositionHistogramSize[CurrentSystem]);

    // samples the free energy profiles in a,b,c directions
    if(strcasecmp("WriteFreeEnergyProfileEvery",keyword)==0) sscanf(arguments,"%d",&WriteFreeEnergyProfileEvery[CurrentSystem]);
    if(strcasecmp("FreeEnergyHistogramSize",keyword)==0) sscanf(arguments,"%d",&FreeEnergyHistogramSize[CurrentSystem]);

    // sampling the pore-size distribution (PSD)
    if((strcasecmp("ComputePSDHistogram",keyword)==0)||(strcasecmp("ComputePSD",keyword)==0))
    {
      if(strcasecmp("yes",firstargument)==0) ComputePSDHistogram[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputePSDHistogram[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WritePSDHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WritePSDHistogramEvery[CurrentSystem]);
    if(strcasecmp("PSDHistogramSize",keyword)==0) sscanf(arguments,"%d",&PSDHistogramSize[CurrentSystem]);
    if(strcasecmp("PSDRange",keyword)==0) sscanf(arguments,"%lf",&PSDRange[CurrentSystem]);
    if(strcasecmp("PSDProbeDistance",keyword)==0) 
    {
      if(strcasecmp("Minimum",firstargument)==0) Framework[CurrentSystem].PoreSizeDistributionProbeDistance=pow(2.0,1.0/6.0);
      if(strcasecmp("Sigma",firstargument)==0) Framework[CurrentSystem].PoreSizeDistributionProbeDistance=1.0;
    }

    // sampling the end-to-end histograms
    if(strcasecmp("ComputeEndToEndDistanceHistogram",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeEndToEndDistanceHistogram[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeEndToEndDistanceHistogram[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteEndToEndDistanceHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteEndToEndDistanceHistogramEvery[CurrentSystem]);
    if(strcasecmp("EndToEndHistogramSize",keyword)==0) sscanf(arguments,"%d",&EndToEndHistogramSize[CurrentSystem]);
    if(strcasecmp("EndToEndRange",keyword)==0) sscanf(arguments,"%lf",&EndToEndRange[CurrentSystem]);


    // sampling the energy histogram
    if(strcasecmp("ComputeEnergyHistogram",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeEnergyHistogram[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeEnergyHistogram[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteEnergyHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteEnergyHistogramEvery[CurrentSystem]);
    if(strcasecmp("EnergyHistogramSize",keyword)==0) sscanf(arguments,"%d",&EnergyHistogramSize[CurrentSystem]);
    if(strcasecmp("EnergyHistogramLowerLimit",keyword)==0) sscanf(arguments,"%lf",&EnergyHistogramLowerLimit[CurrentSystem]);
    if(strcasecmp("EnergyHistogramUpperLimit",keyword)==0) sscanf(arguments,"%lf",&EnergyHistogramUpperLimit[CurrentSystem]);

    // sampling the thermodynamic factor
    if(strcasecmp("ComputeThermoDynamicFactor",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeThermoDynamicFactor[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeThermoDynamicFactor[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteThermoDynamicFactorEvery",keyword)==0) sscanf(arguments,"%d",&WriteThermoDynamicFactorEvery[CurrentSystem]);

    // sampling the inter-framework spacing histogram
    if(strcasecmp("ComputeFrameworkSpacingHistogram",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeFrameworkSpacingHistogram[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeFrameworkSpacingHistogram[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteFrameworkSpacingHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteFrameworkSpacingHistogramEvery[CurrentSystem]);
    if(strcasecmp("FrameworkSpacingHistogramSize",keyword)==0) sscanf(arguments,"%d",&FrameworkSpacingHistogramSize[CurrentSystem]);
    if(strcasecmp("FrameworkSpacingRange",keyword)==0) sscanf(arguments,"%lf",&FrameworkSpacingRange[CurrentSystem]);

    // sampling histograms of the residence times
    if(strcasecmp("ComputeResidenceTimes",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeResidenceTimes[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeResidenceTimes[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteResidenceTimesEvery",keyword)==0) sscanf(arguments,"%d",&WriteResidenceTimesEvery[CurrentSystem]);
    if(strcasecmp("RangeResidenceTimes",keyword)==0) sscanf(arguments,"%lf",&RangeResidenceTimes[CurrentSystem]);
    if(strcasecmp("ResidenceTimesHistogramSize",keyword)==0) sscanf(arguments,"%d",&ResidenceTimesHistogramSize[CurrentSystem]);

    // sampling histograms of the distance between 2 selected atoms
    if(strcasecmp("ComputeDistanceHistograms",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeDistanceHistograms[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeDistanceHistograms[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteDistanceHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteDistanceHistogramsEvery[CurrentSystem]);
    if(strcasecmp("MaxRangeDistanceHistogram",keyword)==0) sscanf(arguments,"%lf",&MaxRangeDistanceHistogram);
    if(strcasecmp("NumberOfElementsDistanceHistogram",keyword)==0) sscanf(arguments,"%d",&NumberOfElementsDistanceHistogram);

    // sampling histograms of the bend angle between 3 selected atoms
    if(strcasecmp("ComputeBendAngleHistograms",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeBendAngleHistograms[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeBendAngleHistograms[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteBendAngleHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteBendAngleHistogramsEvery[CurrentSystem]);
    if(strcasecmp("MaxRangeBendAngleHistogram",keyword)==0) sscanf(arguments,"%lf",&MaxRangeBendAngleHistogram);
    if(strcasecmp("NumberOfElementsBendAngleHistogram",keyword)==0) sscanf(arguments,"%d",&NumberOfElementsBendAngleHistogram);

    // sampling histograms of the dihedral angle between 4 selected atoms
    if(strcasecmp("ComputeDihedralAngleHistograms",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeDihedralAngleHistograms[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeDihedralAngleHistograms[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteDihedralAngleHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteDihedralAngleHistogramsEvery[CurrentSystem]);
    if(strcasecmp("MaxRangeDihedralAngleHistogram",keyword)==0) sscanf(arguments,"%lf",&MaxRangeDihedralAngleHistogram);
    if(strcasecmp("NumberOfElementsDihedralAngleHistogram",keyword)==0) sscanf(arguments,"%d",&NumberOfElementsDihedralAngleHistogram);

    // sampling histograms of the angle between two planes (each formed by 3 chosen atoms)
    if(strcasecmp("ComputeAngleBetweenPlanesHistograms",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeAngleBetweenPlanesHistograms[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeAngleBetweenPlanesHistograms[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteAngleBetweenPlanesHistogramEvery",keyword)==0) sscanf(arguments,"%d",&WriteAngleBetweenPlanesHistogramsEvery[CurrentSystem]);
    if(strcasecmp("MaxRangeAngleBetweenPlanesHistogram",keyword)==0) sscanf(arguments,"%lf",&MaxRangeAngleBetweenPlanesHistogram);
    if(strcasecmp("NumberOfElementsAngleBetweenPlanesHistogram",keyword)==0) sscanf(arguments,"%d",&NumberOfElementsAngleBetweenPlanesHistogram);

    // sampling molecular properties (bond distance, bend angle, dihedral angle)
    if(strcasecmp("ComputeMoleculeProperties",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeMoleculeProperties[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeMoleculeProperties[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteMoleculePropertiesEvery",keyword)==0) sscanf(arguments,"%d",&WriteMoleculePropertiesEvery[CurrentSystem]);
    if(strcasecmp("BondLengthHistogramSize",keyword)==0) sscanf(arguments,"%d",&BondLengthHistogramSize[CurrentSystem]);
    if(strcasecmp("BendAngleHistogramSize",keyword)==0) sscanf(arguments,"%d",&BendAngleHistogramSize[CurrentSystem]);
    if(strcasecmp("DihedralHistogramSize",keyword)==0) sscanf(arguments,"%d",&DihedralHistogramSize[CurrentSystem]);
    if(strcasecmp("BondLengthRange",keyword)==0) sscanf(arguments,"%lf",&BondLengthRange[CurrentSystem]);
    if(strcasecmp("BendAngleRange",keyword)==0) sscanf(arguments,"%lf",&BendAngleRange[CurrentSystem]);
    if(strcasecmp("DihedralRange",keyword)==0) sscanf(arguments,"%lf",&DihedralRange[CurrentSystem]);

    // sampling the IR spectra (spacings: 2048, 4196, 8192, 16384, 32768 points)
    if(strcasecmp("ComputeInfraRedSpectra",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeInfraRedSpectra[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeInfraRedSpectra[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteInfraRedSpectraEvery",keyword)==0) sscanf(arguments,"%d",&WriteInfraRedSpectraEvery[CurrentSystem]);
    if(strcasecmp("SampleEveryInfraRed",keyword)==0) sscanf(arguments,"%d",&SampleEveryInfraRed);

    // sampling the mean-squared displacement using a modified order-N algorithm
    if(strcasecmp("ComputeMSD",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputeMSDOrderN[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeMSDOrderN[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteMSDEvery",keyword)==0) sscanf(arguments,"%d",&WriteMSDOrderNEvery[CurrentSystem]);
    if(strcasecmp("SampleMSDEvery",keyword)==0) sscanf(arguments,"%d",&SampleMSDOrderNEvery[CurrentSystem]);
    if(strcasecmp("NumberOfBlockElementsMSD",keyword)==0) sscanf(arguments,"%d",&NumberOfBlockElementsMSDOrderN);
    if(strcasecmp("NumberOfBlocksMSD",keyword)==0) sscanf(arguments,"%d",&MaxNumberOfBlocksMSDOrderN);
    if(strcasecmp("ComputeIndividualMSD",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputeIndividualMSDOrderN=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeIndividualMSDOrderN=FALSE;
    }
    if(strcasecmp("ComputeMSDPerPseudoAtom",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeMSDOrderNPerPseudoAtom=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeMSDOrderNPerPseudoAtom=FALSE;
    }

    // sampling the velocity autocorrelation function using a modified order-N algorithm
    if(strcasecmp("ComputeVACF",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeVACFOrderN[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeVACFOrderN[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteVACFEvery",keyword)==0) sscanf(arguments,"%d",&WriteVACFOrderNEvery[CurrentSystem]);
    if(strcasecmp("SampleVACFEvery",keyword)==0) sscanf(arguments,"%d",&SampleVACFOrderNEvery[CurrentSystem]);
    if(strcasecmp("NumberOfBlockElementsVACF",keyword)==0) sscanf(arguments,"%d",&NumberOfBlockElementsVACFOrderN);
    if(strcasecmp("NumberOfBlocksVACF",keyword)==0) sscanf(arguments,"%d",&MaxNumberOfBlocksVACFOrderN);
    if(strcasecmp("ComputeIndividualVACF",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeIndividualVACFOrderN=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeIndividualVACFOrderN=FALSE;
    }
    if(strcasecmp("ComputeVACFPerPseudoAtom",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeVACFOrderNPerPseudoAtom=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeVACFOrderNPerPseudoAtom=FALSE;
    }

    // sampling of the rotational velocity autocorrelation function using a modified order-N algorithm
    if(strcasecmp("ComputeRVACF",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeRVACFOrderN[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeRVACFOrderN[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteRVACFEvery",keyword)==0) sscanf(arguments,"%d",&WriteRVACFOrderNEvery[CurrentSystem]);
    if(strcasecmp("SampleRVACFEvery",keyword)==0) sscanf(arguments,"%d",&SampleRVACFOrderNEvery[CurrentSystem]);
    if(strcasecmp("NumberOfBlockElementsRVACF",keyword)==0) sscanf(arguments,"%d",&NumberOfBlockElementsRVACFOrderN);
    if(strcasecmp("NumberOfBlocksRVACF",keyword)==0) sscanf(arguments,"%d",&MaxNumberOfBlocksRVACFOrderN);
    if(strcasecmp("ComputeIndividualRVACF",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeIndividualRVACFOrderN=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeIndividualRVACFOrderN=FALSE;
    }
    if(strcasecmp("ComputeRVACFPerPseudoAtom",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeRVACFOrderNPerPseudoAtom=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeRVACFOrderNPerPseudoAtom=FALSE;
    }

    // sampling of the molecular orientation autocorrelation function using a modified order-N algorithm
    if(strcasecmp("ComputeMOACF",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeMolecularOrientationOrderN[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeMolecularOrientationOrderN[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteMOAACFEvery",keyword)==0) sscanf(arguments,"%d",&WriteMolecularOrientationOrderNEvery[CurrentSystem]);
    if(strcasecmp("SampleMOACFEvery",keyword)==0) sscanf(arguments,"%d",&SampleMolecularOrientationOrderNEvery[CurrentSystem]);
    if(strcasecmp("NumberOfBlockElementsMOACF",keyword)==0) sscanf(arguments,"%d",&NumberOfBlockElementsMolecularOrientationOrderN);
    if(strcasecmp("NumberOfBlocksMOACF",keyword)==0) sscanf(arguments,"%d",&MaxNumberOfBlocksMolecularOrientationOrderN);

    // sampling the mean-square displacement function using a conventional algorithm
    if(strcasecmp("ComputeMSDConventional",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeMSD[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeMSD[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteMSDConventionalEvery",keyword)==0) sscanf(arguments,"%d",&WriteMSDEvery[CurrentSystem]);
    if(strcasecmp("SampleMSDConventionalEvery",keyword)==0) sscanf(arguments,"%d",&SampleMSDEvery[CurrentSystem]);
    if(strcasecmp("NumberOfBuffersMSDConventional",keyword)==0) sscanf(arguments,"%d",&NumberOfBuffersMSD);
    if(strcasecmp("BufferLengthMSDConventional",keyword)==0) sscanf(arguments,"%d",&BufferLengthMSD);

    // sampling the velocity autocorrelation function using a conventional algorithm
    if(strcasecmp("ComputeVACFConventional",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeVACF[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeVACF[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteVACFConventionalEvery",keyword)==0) sscanf(arguments,"%d",&WriteVACFEvery[CurrentSystem]);
    if(strcasecmp("SampleVACFConventionalEvery",keyword)==0) sscanf(arguments,"%d",&SampleVACFEvery[CurrentSystem]);
    if(strcasecmp("NumberOfBuffersVACFConventional",keyword)==0) sscanf(arguments,"%d",&NumberOfBuffersVACF);
    if(strcasecmp("BufferLengthVACFConventional",keyword)==0) sscanf(arguments,"%d",&BufferLengthVACF);

    // sampling the 3D histograms of position (i.e. 3D free energy)
    if(strcasecmp("ComputeDensityProfile3DVTKGrid",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputeDensityProfile3DVTKGrid[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeDensityProfile3DVTKGrid[CurrentSystem]=FALSE;
    }
    if(strcasecmp("WriteDensityProfile3DVTKGridEvery",keyword)==0) sscanf(arguments,"%d",&WriteDensityProfile3DVTKGridEvery[CurrentSystem]);
    if(strcasecmp("DensityProfile3DVTKGridPoints",keyword)==0) 
      sscanf(arguments,"%d %d %d",&DensityProfile3DVTKGridPoints.x,&DensityProfile3DVTKGridPoints.y,&DensityProfile3DVTKGridPoints.z);

    // samples the cation sites and adsorption sites
    if(strcasecmp("ComputeCationAndAdsorptionSites",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeCationAndAdsorptionSites[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeCationAndAdsorptionSites[CurrentSystem]=FALSE;
    }

    if(strcasecmp("WriteCationAndAdsorptionSitesEvery",keyword)==0) sscanf(arguments,"%d",&WriteCationAndAdsorptionSitesEvery[CurrentSystem]);

    // samples initial configurations for the tranmission coefficient (dcTST)
    if(strcasecmp("WritedcTSTSnapShotsEvery",keyword)==0) sscanf(arguments,"%d",&WritedcTSTSnapShotsEvery[CurrentSystem]);
    if(strcasecmp("WritedcTSTSnapShotsToFile",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) WritedcTSTSnapShotsToFile[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) WritedcTSTSnapShotsToFile[CurrentSystem]=FALSE;
    }

    // samples the principle moment of inertia
    if(strcasecmp("ComputePrincipleMomentsOfInertia",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputePrincipleMomentsOfInertia=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputePrincipleMomentsOfInertia=FALSE;
    }

    if((strcasecmp("ComputeMolecularPressure",keyword)==0)||(strcasecmp("ComputePressure",keyword)==0))
    {
      if(strcasecmp("yes",firstargument)==0) ComputeMolecularPressure[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeMolecularPressure[CurrentSystem]=FALSE;
    }



    if(strcasecmp("SurfaceAreaProbeAtom",keyword)==0) sscanf(arguments,"%s",Framework[CurrentSystem].SurfaceAreaProbeAtom);
    if(strcasecmp("SurfaceAreaSamplingPointsPerSphere",keyword)==0) sscanf(arguments,"%d",&Framework[CurrentSystem].SurfaceAreaSamplingPointsPerShere);
    if(strcasecmp("SurfaceAreaProbeDistance",keyword)==0) 
    {
      if(strcasecmp("Minimum",firstargument)==0) Framework[CurrentSystem].SurfaceAreaProbeDistance=pow(2.0,1.0/6.0);
      if(strcasecmp("Sigma",firstargument)==0) Framework[CurrentSystem].SurfaceAreaProbeDistance=1.0;
    }

    if(strcasecmp("ComputeElasticConstants",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputeElasticConstants=ComputeBornTerm=ComputeCrossTerm=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeElasticConstants=ComputeBornTerm=ComputeCrossTerm=FALSE;
    }
    if(strcasecmp("ComputeElasticConstantsThirdOrder",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputeElasticConstantsThirdOrder=ComputeElasticConstants=ComputeBornTerm=ComputeCrossTerm=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeElasticConstantsThirdOrder=FALSE;
    }
    if(strcasecmp("ComputePowderDiffractionPattern",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputePowderDiffractionPattern=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputePowderDiffractionPattern=FALSE;
    }
    if(strcasecmp("DiffractionType",keyword)==0) 
    {
      if(strcasecmp("xray",firstargument)==0) Diffraction.Type=XRAY_DIFFRACTION;
      if(strcasecmp("neutron",firstargument)==0) Diffraction.Type=NEUTRON_DIFFRACTION;
      if(strcasecmp("electron",firstargument)==0) Diffraction.Type=ELECTRON_DIFFRACTION;
    }
    if(strcasecmp("DiffractionRadiationType",keyword)==0) 
    {
      if(strcasecmp("chromium",firstargument)==0) Diffraction.RadiationType=CHROMIUM_RADIATION;
      if(strcasecmp("iron",firstargument)==0) Diffraction.RadiationType=IRON_RADIATION;
      if(strcasecmp("copper",firstargument)==0) Diffraction.RadiationType=COPPER_RADIATION;
      if(strcasecmp("molybdenum",firstargument)==0) Diffraction.RadiationType=MOLYBDENUM_RADIATION;
      if(strcasecmp("silver",firstargument)==0) Diffraction.RadiationType=SILVER_RADIATION;
      if(strcasecmp("synchrotron",firstargument)==0) Diffraction.RadiationType=SYNCHROTRON_RADIATION;
    }
    if(strcasecmp("WaveLengthType",keyword)==0) 
    {
      if(strcasecmp("Single",firstargument)==0) Diffraction.lambda_type=DIFFRACTION_SINGLE;
      if(strcasecmp("Double",firstargument)==0) Diffraction.lambda_type=DIFFRACTION_DOUBLET;
    }
    if(strcasecmp("PeakShape",keyword)==0) 
    {
      if(strcasecmp("Gaussian",firstargument)==0) Diffraction.PeakShape=DIFFRACTION_GAUSSIAN;
      if(strcasecmp("Lorentzian",firstargument)==0) Diffraction.PeakShape=DIFFRACTION_LORENTZIAN;
      if(strcasecmp("PseudoVoigt",firstargument)==0) Diffraction.PeakShape=DIFFRACTION_PSEUDO_VOIGT;
    }
    if(strcasecmp("WaveLength",keyword)==0) sscanf(arguments,"%lf",&Wavelength);
    if(strcasecmp("TwoThetaMin",keyword)==0) sscanf(arguments,"%lf",&Diffraction.two_theta_min);
    if(strcasecmp("TwoThetaMax",keyword)==0) sscanf(arguments,"%lf",&Diffraction.two_theta_max);
    if(strcasecmp("TwoThetaStep",keyword)==0) sscanf(arguments,"%lf",&Diffraction.two_theta_step);
    if(strcasecmp("PeakWidthModifierU",keyword)==0) sscanf(arguments,"%lf",&Diffraction.u);
    if(strcasecmp("PeakWidthModifierV",keyword)==0) sscanf(arguments,"%lf",&Diffraction.v);
    if(strcasecmp("PeakWidthModifierW",keyword)==0) sscanf(arguments,"%lf",&Diffraction.w);


    if(strcasecmp("ComputeNormalModes",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) ComputeNormalModes=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeNormalModes=FALSE;
    }
    if(strcasecmp("CorrectNormalModesForConstraints",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) CorrectNormalModesForConstraints=TRUE;
      if(strcasecmp("no",firstargument)==0) CorrectNormalModesForConstraints=FALSE;
    }
    if(strcasecmp("MinimumMode",keyword)==0) sscanf(arguments,"%d",&MinimumMode);
    if(strcasecmp("MaximumMode",keyword)==0) sscanf(arguments,"%d",&MaximumMode);
    if(strcasecmp("ModeResolution",keyword)==0) sscanf(arguments,"%d",&ModeResolution);

    if(strcasecmp("WriteVTKGrids",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) WriteVTKGrids=TRUE;
      if(strcasecmp("no",firstargument)==0) WriteVTKGrids=FALSE;
    }

	  // read energy/force grid options
	  if(strcasecmp("UseTabularGrid",keyword)==0) 
		  if((strcasecmp(firstargument,"yes")==0)&&(SimulationType!=MAKE_GRID)) UseTabularGrid=TRUE;
	  if(strcasecmp("UseDynamicGrid",keyword)==0) 
		  if((strcasecmp(firstargument,"yes")==0)&&(SimulationType!=MAKE_GRID)) UseDynamicGrid=TRUE;
	  if(strcasecmp("UseDelaunayGrid",keyword)==0) 
		  if((strcasecmp(firstargument,"yes")==0)&&(SimulationType!=MAKE_GRID)) UseDelaunayGrid=TRUE;
	  if(strcasecmp("SpacingVDWGrid",keyword)==0) sscanf(arguments,"%lf",&SpacingVDWGrid);
	  if(strcasecmp("SpacingDelaunayGrid",keyword)==0) sscanf(arguments,"%lf",&SpacingDelaunayGrid);
	  if(strcasecmp("SpacingDynamicGrid",keyword)==0) sscanf(arguments,"%lf",&SpacingDynamicGrid);
	  if(strcasecmp("SpacingCoulombGrid",keyword)==0) sscanf(arguments,"%lf",&SpacingCoulombGrid);
	  if(strcasecmp("NumberOfGrids",keyword)==0)
    {
      sscanf(arguments,"%d",&NumberOfGrids);
      GridTypeListName=(char(*)[256])calloc(NumberOfGrids,sizeof(char[256]));
    }
    if(strcasecmp("GridTypes",keyword)==0) 
    {
      for(j=0;j<NumberOfGrids;j++)
        sscanf(arguments,"%s %[^\n]",GridTypeListName[j],arguments);
    }
    if(strcasecmp("BlockEnergyGrids",keyword)==0)
    {
      if(strcasecmp(firstargument,"yes")==0) BlockEnergyGrids=TRUE;
      if(strcasecmp(firstargument,"no")==0) BlockEnergyGrids=FALSE;
    }
    if(strcasecmp("BlockGridPockets",keyword)==0)
    {
      if(strcasecmp(firstargument,"yes")==0) BlockGridPockets=TRUE;
      if(strcasecmp(firstargument,"no")==0) BlockGridPockets=FALSE;
    }
    if(strcasecmp("BlockGridPores",keyword)==0)
    {
      if(strcasecmp(firstargument,"yes")==0) BlockGridPores=TRUE;
      if(strcasecmp(firstargument,"no")==0) BlockGridPores=FALSE;
    }
    if(strcasecmp("BlockEnergyGridOverlapCriteria",keyword)==0) sscanf(arguments,"%lf",&BlockEnergyGridOverlapCriteria);


    if(strcasecmp("GridSeed",keyword)==0)
    {
      GridSeeds=realloc(GridSeeds,(NumberOfGridSeeds+1)*sizeof(VECTOR));
      sscanf(arguments,"%lf %lf %lf",&GridSeeds[NumberOfGridSeeds].x,&GridSeeds[NumberOfGridSeeds].y,&GridSeeds[NumberOfGridSeeds].z);
		NumberOfGridSeeds++;
    }
	  
	  if(strcasecmp("UseDelaunayPolyhedraFile",keyword)==0)
	  {
		  if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
		  {
			  if(strcasecmp("yes",keyword)==0) DelaunayPolyhedraFile=TRUE;
			  if(strcasecmp("no",keyword)==0) DelaunayPolyhedraFile=FALSE;
		  }
	  }
	  
	  if(strcasecmp("DelaunayPolyhedraFileName",keyword)==0)
	  {
		  if(sscanf(arguments,"%s %[^\n]",keyword,arguments)>0)
		  {
			  strcpy(DelaunayPolyhedraFileName,keyword);
		  }
	  }
	  
	  // read method of minimization
    if(strcasecmp("MinimizationMethod",keyword)==0) 
    {
      if(strcasecmp("SteepestDescent",firstargument)==0) MinimizationMethod=STEEPEST_DESCENT_MINIMIZATION;
      if(strcasecmp("ConjugateGradient",firstargument)==0) MinimizationMethod=CONJUGATE_GRADIENT_MINIMIZATION;
      if(strcasecmp("bfgc",firstargument)==0) MinimizationMethod=BFGS_MINIMIZATION;
      if(strcasecmp("Snyman",firstargument)==0) MinimizationMethod=SNYMAN_MINIMIZATION;
      if((strcasecmp("Baker",firstargument)==0)||(strcasecmp("BakerMinimization",firstargument)==0)) MinimizationMethod=BAKER_MINIMIZATION;
      if((strcasecmp("NewtonRaphson",firstargument)==0)||(strcasecmp("NewtonRaphsonMinimization",firstargument)==0)) MinimizationMethod=NEWTON_RAPHSON_MINIMIZATION;
      if((strcasecmp("BakerSaddlePoint",firstargument)==0)||(strcasecmp("SaddlePoint",firstargument)==0)) MinimizationMethod=BAKER_SADDLE_POINT;
    }
    if(strcasecmp("MinimizationVariables",keyword)==0) 
    {
      if(strcasecmp("Fractional",firstargument)==0) MinimizationVariables=FRACTIONAL;
      if(strcasecmp("Cartesian",firstargument)==0) MinimizationVariables=CARTESIAN;
    }
    if(strcasecmp("UseSymmetryInMinimization",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) UseSymmetryInMinimization=TRUE;
      if(strcasecmp("no",firstargument)==0) UseSymmetryInMinimization=FALSE;
    }

    if(strcasecmp("MaximumNumberOfMinimizationSteps",keyword)==0) sscanf(arguments,"%d",&MaximumNumberOfMinimizationSteps);
    if(strcasecmp("RMSGradientTolerance",keyword)==0) sscanf(arguments,"%lf",&RMSGradientTolerance);
    if(strcasecmp("MaxGradientTolerance",keyword)==0) sscanf(arguments,"%lf",&MaxGradientTolerance);
    if(strcasecmp("MaximumStepLength",keyword)==0) 
    {
      sscanf(arguments,"%lf",&MaximumStepLengthInput);
      MaximumStepLength=MaximumStepLengthInput;
    }
    if(strcasecmp("MinimizationConvergenceFactor",keyword)==0) 
      sscanf(arguments,"%lf",&MinimizationConvergenceFactor);
    if(strcasecmp("UseGradientInLineMinimization",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) UseGradientInLineMinimization=TRUE;
      if(strcasecmp("no",firstargument)==0) UseGradientInLineMinimization=FALSE;
    }
    if(strcasecmp("TransformUnitCell",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) TransformUnitCell=TRUE;
      if(strcasecmp("no",firstargument)==0) TransformUnitCell=FALSE;
    }
    if(strcasecmp("UnitCellDeformation",keyword)==0) sscanf(arguments,"%lf",&UnitCellDeformation);
    if(strcasecmp("ElasticConstantEnergyVolume",keyword)==0) 
    {
      if(strcasecmp("CUBIC_C44",firstargument)==0) ElasticConstantEnergyVolume=ELASTIC_CONSTANT_CUBIC_C44;
      if(strcasecmp("CUBIC_CS",firstargument)==0) ElasticConstantEnergyVolume=ELASTIC_CONSTANT_CUBIC_CS;
      if(strcasecmp("CUBIC_C66",firstargument)==0) ElasticConstantEnergyVolume=ELASTIC_CONSTANT_CUBIC_C66;
      if(strcasecmp("HEXAGONAL_C44",firstargument)==0) ElasticConstantEnergyVolume=ELASTIC_CONSTANT_HEXAGONAL_C44;
      if(strcasecmp("BULK_MODULUS",firstargument)==0) ElasticConstantEnergyVolume=ELASTIC_CONSTANT_BULK_MODULUS;
    }
    if(strcasecmp("MinimizationPotentialMethod",keyword)==0) 
    {
      if(strcasecmp("Analytically",firstargument)==0) MinimizationPotentialMethod=ANALYTICALLY;
      if(strcasecmp("Numerically",firstargument)==0) MinimizationPotentialMethod=NUMERICALLY;
    }
    if(strcasecmp("ComputeLambda",keyword)==0) 
    {
      if(strcasecmp("LambdaMethod1",firstargument)==0) ComputeLambda=LAMBDA_METHOD_1;
      if(strcasecmp("LambdaMethod2",firstargument)==0) ComputeLambda=LAMBDA_METHOD_2;
    }
    if(strcasecmp("RemoveTranslationFromHessian",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) RemoveTranslationFromHessian=TRUE;
      if(strcasecmp("no",firstargument)==0) RemoveTranslationFromHessian=FALSE;
    }
    if(strcasecmp("RemoveRotationFromHessian",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) RemoveRotationFromHessian=TRUE;
      if(strcasecmp("no",firstargument)==0) RemoveRotationFromHessian=FALSE;
    }
    if(strcasecmp("FrameworkFixedInitialization",keyword)==0)
    {
      if(strcasecmp("free",firstargument)==0) FrameworkFixedInitialization[CurrentSystem][CurrentFramework]=FREE;
      if(strcasecmp("fixed",firstargument)==0) FrameworkFixedInitialization[CurrentSystem][CurrentFramework]=FIXED;
    } 
    if(strcasecmp("AdsorbateFixedInitialization",keyword)==0)
    {
      if(strcasecmp("free",firstargument)==0) AdsorbateFixedInitialization[CurrentSystem]=FREE;
      if(strcasecmp("fixed",firstargument)==0) AdsorbateFixedInitialization[CurrentSystem]=FIXED;
    }
    if(strcasecmp("CationFixedInitialization",keyword)==0)
    {
      if(strcasecmp("free",firstargument)==0) CationFixedInitialization[CurrentSystem]=FREE;
      if(strcasecmp("fixed",firstargument)==0) CationFixedInitialization[CurrentSystem]=FIXED;
    }
    if(strcasecmp("ActiveFrameworkAtom",keyword)==0)
    {
      ActiveFrameworkAtoms[CurrentSystem][CurrentFramework]=(int*)realloc(ActiveFrameworkAtoms[CurrentSystem][CurrentFramework],
                           (NumberOfActiveFrameworkAtoms[CurrentSystem][CurrentFramework]+1)*sizeof(int));
      sscanf(arguments,"%d",&ActiveFrameworkAtoms[CurrentSystem][CurrentFramework][NumberOfActiveFrameworkAtoms[CurrentSystem][CurrentFramework]]);
      NumberOfActiveFrameworkAtoms[CurrentSystem][CurrentFramework]++;
    }
    if(strcasecmp("ActiveFrameworkAtoms",keyword)==0)
    {
      sscanf(arguments,"%d%[^\n]",&temp_int,arguments);
      ActiveFrameworkAtoms[CurrentSystem][CurrentFramework]=(int*)realloc(ActiveFrameworkAtoms[CurrentSystem][CurrentFramework],
                           (NumberOfActiveFrameworkAtoms[CurrentSystem][CurrentFramework]+temp_int)*sizeof(int));
      for(i=0;i<temp_int;i++)
        sscanf(arguments,"%d%[^\n]",&ActiveFrameworkAtoms[CurrentSystem][CurrentFramework][NumberOfActiveFrameworkAtoms[CurrentSystem][CurrentFramework]++],arguments);
    }
    if(strcasecmp("FixedFrameworkAtom",keyword)==0)
    {
      FixedFrameworkAtoms[CurrentSystem][CurrentFramework]=(int*)realloc(FixedFrameworkAtoms[CurrentSystem][CurrentFramework],
                          (NumberOfFixedFrameworkAtoms[CurrentSystem][CurrentFramework]+1)*sizeof(int));
      sscanf(arguments,"%d",&FixedFrameworkAtoms[CurrentSystem][CurrentFramework][NumberOfFixedFrameworkAtoms[CurrentSystem][CurrentFramework]]);
      NumberOfFixedFrameworkAtoms[CurrentSystem][CurrentFramework]++;
    }
    if(strcasecmp("FixedFrameworkAtoms",keyword)==0)
    {
      sscanf(arguments,"%d%[^\n]",&temp_int,arguments);
      FixedFrameworkAtoms[CurrentSystem][CurrentFramework]=(int*)realloc(FixedFrameworkAtoms[CurrentSystem][CurrentFramework],
                          (NumberOfFixedFrameworkAtoms[CurrentSystem][CurrentFramework]+temp_int)*sizeof(int));
      for(i=0;i<temp_int;i++)
        sscanf(arguments,"%d%[^\n]",&FixedFrameworkAtoms[CurrentSystem][CurrentFramework][NumberOfFixedFrameworkAtoms[CurrentSystem][CurrentFramework]++],arguments);
    }
    if(strcasecmp("ActiveAdsorbateMolecule",keyword)==0)
    {
      ActiveAdsorbateMolecules[CurrentSystem]=(int*)realloc(ActiveAdsorbateMolecules[CurrentSystem],(NumberOfActiveAdsorbateMolecules[CurrentSystem]+1)*sizeof(int));
      sscanf(arguments,"%d",&ActiveAdsorbateMolecules[CurrentSystem][NumberOfActiveAdsorbateMolecules[CurrentSystem]]);
      NumberOfActiveAdsorbateMolecules[CurrentSystem]++;
    }
    if(strcasecmp("FixedAdsorbateMolecule",keyword)==0)
    {
      FixedAdsorbateMolecules[CurrentSystem]=(int*)realloc(FixedAdsorbateMolecules[CurrentSystem],(NumberOfFixedAdsorbateMolecules[CurrentSystem]+1)*sizeof(int));
      sscanf(arguments,"%d",&FixedAdsorbateMolecules[CurrentSystem][NumberOfFixedAdsorbateMolecules[CurrentSystem]]);
      NumberOfFixedAdsorbateMolecules[CurrentSystem]++;
    }
    if(strcasecmp("ActiveAdsorbateAtom",keyword)==0)
    {
      ActiveAdsorbateAtoms[CurrentSystem]=(PAIR*)realloc(ActiveAdsorbateAtoms[CurrentSystem],(NumberOfActiveAdsorbateAtoms[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveAdsorbateAtoms[CurrentSystem][NumberOfActiveAdsorbateAtoms[CurrentSystem]].A,
                       &ActiveAdsorbateAtoms[CurrentSystem][NumberOfActiveAdsorbateAtoms[CurrentSystem]].B);
      NumberOfActiveAdsorbateAtoms[CurrentSystem]++;
    }
    if(strcasecmp("FixedAdsorbateAtom",keyword)==0)
    {
      FixedAdsorbateAtoms[CurrentSystem]=(PAIR*)realloc(FixedAdsorbateAtoms[CurrentSystem],(NumberOfFixedAdsorbateAtoms[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedAdsorbateAtoms[CurrentSystem][NumberOfFixedAdsorbateAtoms[CurrentSystem]].A,
                       &FixedAdsorbateAtoms[CurrentSystem][NumberOfFixedAdsorbateAtoms[CurrentSystem]].B);
      NumberOfFixedAdsorbateAtoms[CurrentSystem]++;
    }
    if(strcasecmp("ActiveAdsorbateGroup",keyword)==0)
    {
      ActiveAdsorbateGroups[CurrentSystem]=(PAIR*)realloc(ActiveAdsorbateGroups[CurrentSystem],(NumberOfActiveAdsorbateGroups[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveAdsorbateGroups[CurrentSystem][NumberOfActiveAdsorbateGroups[CurrentSystem]].A,
                       &ActiveAdsorbateGroups[CurrentSystem][NumberOfActiveAdsorbateGroups[CurrentSystem]].B);
      NumberOfActiveAdsorbateGroups[CurrentSystem]++;
    }
    if(strcasecmp("FixedAdsorbateGroup",keyword)==0)
    {
      FixedAdsorbateGroups[CurrentSystem]=(PAIR*)realloc(FixedAdsorbateGroups[CurrentSystem],(NumberOfFixedAdsorbateGroups[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedAdsorbateGroups[CurrentSystem][NumberOfFixedAdsorbateGroups[CurrentSystem]].A,
                       &FixedAdsorbateGroups[CurrentSystem][NumberOfFixedAdsorbateGroups[CurrentSystem]].B);
      NumberOfFixedAdsorbateGroups[CurrentSystem]++;
    }
    if(strcasecmp("ActiveAdsorbateGroupCenterOfMass",keyword)==0)
    {
      ActiveAdsorbateGroupsCenterOfMass[CurrentSystem]=(PAIR*)realloc(ActiveAdsorbateGroupsCenterOfMass[CurrentSystem],
                       (NumberOfActiveAdsorbateGroupsCenterOfMass[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveAdsorbateGroupsCenterOfMass[CurrentSystem][NumberOfActiveAdsorbateGroupsCenterOfMass[CurrentSystem]].A,
                       &ActiveAdsorbateGroupsCenterOfMass[CurrentSystem][NumberOfActiveAdsorbateGroupsCenterOfMass[CurrentSystem]].B);
      NumberOfActiveAdsorbateGroupsCenterOfMass[CurrentSystem]++;
    }
    if(strcasecmp("FixedAdsorbateGroupCenterOfMass",keyword)==0)
    {
      FixedAdsorbateGroupsCenterOfMass[CurrentSystem]=(PAIR*)realloc(FixedAdsorbateGroupsCenterOfMass[CurrentSystem],
                                       (NumberOfFixedAdsorbateGroupsCenterOfMass[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedAdsorbateGroupsCenterOfMass[CurrentSystem][NumberOfFixedAdsorbateGroupsCenterOfMass[CurrentSystem]].A,
                       &FixedAdsorbateGroupsCenterOfMass[CurrentSystem][NumberOfFixedAdsorbateGroupsCenterOfMass[CurrentSystem]].B);
      NumberOfFixedAdsorbateGroupsCenterOfMass[CurrentSystem]++;
    }
    if(strcasecmp("ActiveAdsorbateGroupOrientation",keyword)==0)
    {
      ActiveAdsorbateGroupsOrientation[CurrentSystem]=(PAIR*)realloc(ActiveAdsorbateGroupsOrientation[CurrentSystem],
                       (NumberOfActiveAdsorbateGroupsOrientation[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveAdsorbateGroupsOrientation[CurrentSystem][NumberOfActiveAdsorbateGroupsOrientation[CurrentSystem]].A,
                       &ActiveAdsorbateGroupsOrientation[CurrentSystem][NumberOfActiveAdsorbateGroupsOrientation[CurrentSystem]].B);
      NumberOfActiveAdsorbateGroupsOrientation[CurrentSystem]++;
    }
    if(strcasecmp("FixedAdsorbateGroupOrientation",keyword)==0)
    {
      FixedAdsorbateGroupsOrientation[CurrentSystem]=(PAIR*)realloc(FixedAdsorbateGroupsOrientation[CurrentSystem],
                                       (NumberOfFixedAdsorbateGroupsOrientation[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedAdsorbateGroupsOrientation[CurrentSystem][NumberOfFixedAdsorbateGroupsOrientation[CurrentSystem]].A,
                       &FixedAdsorbateGroupsOrientation[CurrentSystem][NumberOfFixedAdsorbateGroupsOrientation[CurrentSystem]].B);
      NumberOfFixedAdsorbateGroupsOrientation[CurrentSystem]++;
    }
    if(strcasecmp("ActiveCationMolecule",keyword)==0)
    {
      ActiveCationMolecules[CurrentSystem]=(int*)realloc(ActiveCationMolecules[CurrentSystem],(NumberOfActiveCationMolecules[CurrentSystem]+1)*sizeof(int));
      sscanf(arguments,"%d",&ActiveCationMolecules[CurrentSystem][NumberOfActiveCationMolecules[CurrentSystem]]);
      NumberOfActiveCationMolecules[CurrentSystem]++;
    }
    if(strcasecmp("FixedAdsorbateMolecule",keyword)==0)
    {
      FixedAdsorbateMolecules[CurrentSystem]=(int*)realloc(FixedAdsorbateMolecules[CurrentSystem],(NumberOfFixedAdsorbateMolecules[CurrentSystem]+1)*sizeof(int));
      sscanf(arguments,"%d",&FixedAdsorbateMolecules[CurrentSystem][NumberOfFixedAdsorbateMolecules[CurrentSystem]]);
      NumberOfFixedAdsorbateMolecules[CurrentSystem]++;
    }
    if(strcasecmp("ActiveCationAtom",keyword)==0)
    {
      ActiveCationAtoms[CurrentSystem]=(PAIR*)realloc(ActiveCationAtoms[CurrentSystem],(NumberOfActiveCationAtoms[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveCationAtoms[CurrentSystem][NumberOfActiveCationAtoms[CurrentSystem]].A,
                       &ActiveCationAtoms[CurrentSystem][NumberOfActiveCationAtoms[CurrentSystem]].B);
      NumberOfActiveCationAtoms[CurrentSystem]++;
    }
    if(strcasecmp("FixedCationAtom",keyword)==0)
    {
      FixedCationAtoms[CurrentSystem]=(PAIR*)realloc(FixedCationAtoms[CurrentSystem],(NumberOfFixedCationAtoms[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedCationAtoms[CurrentSystem][NumberOfFixedCationAtoms[CurrentSystem]].A,
                       &FixedCationAtoms[CurrentSystem][NumberOfFixedCationAtoms[CurrentSystem]].B);
      NumberOfFixedCationAtoms[CurrentSystem]++;
    }
    if(strcasecmp("ActiveCationGroup",keyword)==0)
    {
      ActiveCationGroups[CurrentSystem]=(PAIR*)realloc(ActiveCationGroups[CurrentSystem],(NumberOfActiveCationGroups[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveCationGroups[CurrentSystem][NumberOfActiveCationGroups[CurrentSystem]].A,
                       &ActiveCationGroups[CurrentSystem][NumberOfActiveCationGroups[CurrentSystem]].B);
      NumberOfActiveCationGroups[CurrentSystem]++;
    }
    if(strcasecmp("FixedCationGroup",keyword)==0)
    {
      FixedCationGroups[CurrentSystem]=(PAIR*)realloc(FixedCationGroups[CurrentSystem],(NumberOfFixedCationGroups[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedCationGroups[CurrentSystem][NumberOfFixedCationGroups[CurrentSystem]].A,
                       &FixedCationGroups[CurrentSystem][NumberOfFixedCationGroups[CurrentSystem]].B);
      NumberOfFixedCationGroups[CurrentSystem]++;
    }
    if(strcasecmp("ActiveCationGroupCenterOfMass",keyword)==0)
    {
      ActiveCationGroupsCenterOfMass[CurrentSystem]=(PAIR*)realloc(ActiveCationGroupsCenterOfMass[CurrentSystem],
                       (NumberOfActiveCationGroupsCenterOfMass[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveCationGroupsCenterOfMass[CurrentSystem][NumberOfActiveCationGroupsCenterOfMass[CurrentSystem]].A,
                       &ActiveCationGroupsCenterOfMass[CurrentSystem][NumberOfActiveCationGroupsCenterOfMass[CurrentSystem]].B);
      NumberOfActiveCationGroupsCenterOfMass[CurrentSystem]++;
    }
    if(strcasecmp("FixedCationGroupCenterOfMass",keyword)==0)
    {
      FixedCationGroupsCenterOfMass[CurrentSystem]=(PAIR*)realloc(FixedCationGroupsCenterOfMass[CurrentSystem],
                                       (NumberOfFixedCationGroupsCenterOfMass[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedCationGroupsCenterOfMass[CurrentSystem][NumberOfFixedCationGroupsCenterOfMass[CurrentSystem]].A,
                       &FixedCationGroupsCenterOfMass[CurrentSystem][NumberOfFixedCationGroupsCenterOfMass[CurrentSystem]].B);
      NumberOfFixedCationGroupsCenterOfMass[CurrentSystem]++;
    }
    if(strcasecmp("ActiveCationGroupOrientation",keyword)==0)
    {
      ActiveCationGroupsOrientation[CurrentSystem]=(PAIR*)realloc(ActiveCationGroupsOrientation[CurrentSystem],
                       (NumberOfActiveCationGroupsOrientation[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&ActiveCationGroupsOrientation[CurrentSystem][NumberOfActiveCationGroupsOrientation[CurrentSystem]].A,
                       &ActiveCationGroupsOrientation[CurrentSystem][NumberOfActiveCationGroupsOrientation[CurrentSystem]].B);
      NumberOfActiveCationGroupsOrientation[CurrentSystem]++;
    }
    if(strcasecmp("FixedCationGroupOrientation",keyword)==0)
    {
      FixedCationGroupsOrientation[CurrentSystem]=(PAIR*)realloc(FixedCationGroupsOrientation[CurrentSystem],
                                       (NumberOfFixedCationGroupsOrientation[CurrentSystem]+1)*sizeof(PAIR));
      sscanf(arguments,"%d%d\n",&FixedCationGroupsOrientation[CurrentSystem][NumberOfFixedCationGroupsOrientation[CurrentSystem]].A,
                       &FixedCationGroupsOrientation[CurrentSystem][NumberOfFixedCationGroupsOrientation[CurrentSystem]].B);
      NumberOfFixedCationGroupsOrientation[CurrentSystem]++;
    }
    if(strcasecmp("ActiveAtomType",keyword)==0)
    {
      ActiveAtomTypes=(char(*)[256])realloc(ActiveAtomTypes,(NumberOfActiveAtomTypes+1)*sizeof(char[256]));
      sscanf(arguments,"%s",ActiveAtomTypes[NumberOfActiveAtomTypes]);
      NumberOfActiveAtomTypes++;
    }
    if(strcasecmp("ActiveAtomTypes",keyword)==0)
    {
      sscanf(arguments,"%d%[^\n]",&temp_int,arguments);
      ActiveAtomTypes=(char(*)[256])realloc(ActiveAtomTypes,(NumberOfActiveAtomTypes+temp_int)*sizeof(char[256]));
      for(i=0;i<temp_int;i++)
        sscanf(arguments,"%s%[^\n]",ActiveAtomTypes[NumberOfActiveAtomTypes++],arguments);
    }

    if(strcasecmp("FixedAtomType",keyword)==0)
    {
      FixedAtomTypes=(char(*)[256])realloc(FixedAtomTypes,(NumberOfFixedAtomTypes+1)*sizeof(char[256]));
      sscanf(arguments,"%s",FixedAtomTypes[NumberOfFixedAtomTypes]);
      NumberOfFixedAtomTypes++;
    }
    if(strcasecmp("FixedAtomTypes",keyword)==0)
    {
      sscanf(arguments,"%d%[^\n]",&temp_int,arguments);
      FixedAtomTypes=(char(*)[256])realloc(FixedAtomTypes,(NumberOfFixedAtomTypes+temp_int)*sizeof(char[256]));
      for(i=0;i<temp_int;i++)
        sscanf(arguments,"%s%[^\n]",FixedAtomTypes[NumberOfFixedAtomTypes++],arguments);
    }

    // read the input for the distance histograms between pairs of atoms
    if(strcasecmp("DistanceHistogramDefinition",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,arguments);
      DistanceHistogramDefinitions[CurrentSystem]=(int(*)[2][3])realloc(DistanceHistogramDefinitions[CurrentSystem],(NumberOfDistanceHistogramDefinitions[CurrentSystem]+1)*sizeof(int[2][3]));
      if((charinput1=='F')||(charinput1=='f'))
        DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][0][0]=CATION;
      DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][0][1]=intinput1;
      DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][0][2]=intinput2;
      if((charinput2=='F')||(charinput2=='f'))
        DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][1][0]=CATION;
      DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][1][1]=intinput3;
      DistanceHistogramDefinitions[CurrentSystem][NumberOfDistanceHistogramDefinitions[CurrentSystem]][1][2]=intinput4;

      NumberOfDistanceHistogramDefinitions[CurrentSystem]++;
    }

    // read the input for the bend-angle histograms between trimers of atoms
    if(strcasecmp("BendAngleHistogramDefinition",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d%*[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,&charinput3,&intinput5,&intinput6);
      BendAngleHistogramDefinitions[CurrentSystem]=(int(*)[3][3])realloc(BendAngleHistogramDefinitions[CurrentSystem],
              (NumberOfBendAngleHistogramDefinitions[CurrentSystem]+1)*sizeof(int[3][3]));

      if((charinput1=='F')||(charinput1=='f'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][0][0]=CATION;
      BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][0][1]=intinput1;
      BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][1][0]=CATION;
      BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][1][1]=intinput3;
      BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][2][0]=CATION;
      BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][2][1]=intinput5;
      BendAngleHistogramDefinitions[CurrentSystem][NumberOfBendAngleHistogramDefinitions[CurrentSystem]][2][2]=intinput6;

      NumberOfBendAngleHistogramDefinitions[CurrentSystem]++;
    }

    // read the input for the dihedral-angle histograms between quads of atoms
    if(strcasecmp("DihedralAngleHistogramDefinition",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d%*[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                         &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8);
      DihedralAngleHistogramDefinitions[CurrentSystem]=(int(*)[4][3])realloc(DihedralAngleHistogramDefinitions[CurrentSystem],
              (NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]+1)*sizeof(int[4][3]));

      if((charinput1=='F')||(charinput1=='f'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][0][0]=CATION;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][0][1]=intinput1;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][1][0]=CATION;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][1][1]=intinput3;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][2][0]=CATION;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][2][1]=intinput5;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][3][0]=CATION;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][3][1]=intinput7;
      DihedralAngleHistogramDefinitions[CurrentSystem][NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]][3][2]=intinput8;

      NumberOfDihedralAngleHistogramDefinitions[CurrentSystem]++;
    }

    // read the input for the angle-between-two-planes histograms
    if(strcasecmp("AngleBetweenPlanesHistogramDefinition",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d %c %d %d %c %d %d%*[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                                     &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,
                                                                     &charinput5,&intinput9,&intinput10,&charinput6,&intinput11,&intinput12);
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem]=(int(*)[6][3])realloc(AngleBetweenPlanesHistogramDefinitions[CurrentSystem],
              (NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]+1)*sizeof(int[6][3]));

      if((charinput1=='F')||(charinput1=='f'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][0][0]=CATION;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][0][1]=intinput1;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][1][0]=CATION;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][1][1]=intinput3;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][2][0]=CATION;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][2][1]=intinput5;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][3][0]=CATION;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][3][1]=intinput7;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][3][2]=intinput8;

      if((charinput5=='F')||(charinput5=='f'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][4][0]=FRAMEWORK;
      else if((charinput5=='A')||(charinput5=='a'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][4][0]=ADSORBATE;
      else if((charinput5=='C')||(charinput5=='c'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][4][0]=CATION;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][4][1]=intinput9;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][4][2]=intinput10;

      if((charinput6=='F')||(charinput6=='f'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][5][0]=FRAMEWORK;
      else if((charinput6=='A')||(charinput6=='a'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][5][0]=ADSORBATE;
      else if((charinput6=='C')||(charinput6=='c'))
        AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][5][0]=CATION;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][5][1]=intinput11;
      AngleBetweenPlanesHistogramDefinitions[CurrentSystem][NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]][5][2]=intinput12;

      NumberOfAngleBetweenPlanesHistogramDefinitions[CurrentSystem]++;
    } 

    if(strcasecmp("DistanceConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,arguments);
      sscanf(arguments,"%lf",&realinput1);
      DistanceDefinitions[CurrentSystem]=(int(*)[2][3])realloc(DistanceDefinitions[CurrentSystem],(NumberOfDistanceConstraints[CurrentSystem]+1)*sizeof(int[2][3]));
      if((charinput1=='F')||(charinput1=='f'))
        DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][0][0]=CATION;
      DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][0][1]=intinput1;
      DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][0][2]=intinput2;
      if((charinput2=='F')||(charinput2=='f'))
        DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][1][0]=CATION;
      DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][1][1]=intinput3;
      DistanceDefinitions[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]][1][2]=intinput4;

      // optional reference distance
      DistanceConstraintParameter[CurrentSystem]=(REAL*)realloc(DistanceConstraintParameter[CurrentSystem],
                                                 (NumberOfDistanceConstraints[CurrentSystem]+1)*sizeof(REAL));
      DistanceConstraintParameter[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]]=UNDEFINED;
      if(sscanf(arguments,"%lf",&realinput1)==1)
        DistanceConstraintParameter[CurrentSystem][NumberOfDistanceConstraints[CurrentSystem]]=realinput1;

      NumberOfDistanceConstraints[CurrentSystem]++;
    }

    if(strcasecmp("AngleConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,arguments);
      sscanf(arguments,"%lf",&realinput1);
      AngleDefinitions[CurrentSystem]=(int(*)[3][3])realloc(AngleDefinitions[CurrentSystem],
                                      (NumberOfAngleConstraints[CurrentSystem]+1)*sizeof(int[3][3]));
      if((charinput1=='F')||(charinput1=='f'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][0][0]=CATION;
      AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][0][1]=intinput1;
      AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][1][0]=CATION;
      AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][1][1]=intinput3;
      AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][2][0]=CATION;
      AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][2][1]=intinput5;
      AngleDefinitions[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]][2][2]=intinput6;

      // optional reference angle
      AngleConstraintParameter[CurrentSystem]=(REAL*)realloc(AngleConstraintParameter[CurrentSystem],
                                              (NumberOfAngleConstraints[CurrentSystem]+1)*sizeof(REAL));
      AngleConstraintParameter[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]]=UNDEFINED;
      if(sscanf(arguments,"%lf",&realinput1)==1)
        AngleConstraintParameter[CurrentSystem][NumberOfAngleConstraints[CurrentSystem]]=realinput1*DEG2RAD;

      NumberOfAngleConstraints[CurrentSystem]++;
    }

    if(strcasecmp("DihedralConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,arguments);
      sscanf(arguments,"%lf",&realinput1);
      DihedralDefinitions[CurrentSystem]=(int(*)[4][3])realloc(DihedralDefinitions[CurrentSystem],
                                         (NumberOfDihedralConstraints[CurrentSystem]+1)*sizeof(int[4][3]));
      if((charinput1=='F')||(charinput1=='f'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][0][0]=CATION;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][0][1]=intinput1;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][1][0]=CATION;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][1][1]=intinput3;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][2][0]=CATION;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][2][1]=intinput5;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][3][0]=CATION;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][3][1]=intinput7;
      DihedralDefinitions[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]][3][2]=intinput8;

      // optional reference dihedral angle
      DihedralConstraintParameter[CurrentSystem]=(REAL*)realloc(DihedralConstraintParameter[CurrentSystem],
                                                 (NumberOfDihedralConstraints[CurrentSystem]+1)*sizeof(REAL));
      DihedralConstraintParameter[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]]=UNDEFINED;
      if(sscanf(arguments,"%lf",&realinput1)==1)
        DihedralConstraintParameter[CurrentSystem][NumberOfDihedralConstraints[CurrentSystem]]=realinput1*DEG2RAD;
      
      NumberOfDihedralConstraints[CurrentSystem]++;
    }

    if(strcasecmp("ImproperDihedralConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,arguments);
      sscanf(arguments,"%lf",&realinput1);
      ImproperDihedralDefinitions[CurrentSystem]=(int(*)[4][3])realloc(ImproperDihedralDefinitions[CurrentSystem],
                                         (NumberOfImproperDihedralConstraints[CurrentSystem]+1)*sizeof(int[4][3]));
      if((charinput1=='F')||(charinput1=='f'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][0][0]=CATION;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][0][1]=intinput1;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][1][0]=CATION;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][1][1]=intinput3;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][2][0]=CATION;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][2][1]=intinput5;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][3][0]=CATION;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][3][1]=intinput7;
      ImproperDihedralDefinitions[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]][3][2]=intinput8;

      // optional reference dihedral angle
      ImproperDihedralConstraintParameter[CurrentSystem]=(REAL*)realloc(ImproperDihedralConstraintParameter[CurrentSystem],
                                                 (NumberOfImproperDihedralConstraints[CurrentSystem]+1)*sizeof(REAL));
      ImproperDihedralConstraintParameter[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]]=UNDEFINED;
      if(sscanf(arguments,"%lf",&realinput1)==1)
        ImproperDihedralConstraintParameter[CurrentSystem][NumberOfImproperDihedralConstraints[CurrentSystem]]=realinput1*DEG2RAD;

      NumberOfImproperDihedralConstraints[CurrentSystem]++;
    }

    if(strcasecmp("InversionBendConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,arguments);
      sscanf(arguments,"%lf",&realinput1);
      InversionBendDefinitions[CurrentSystem]=(int(*)[4][3])realloc(InversionBendDefinitions[CurrentSystem],
                                         (NumberOfInversionBendConstraints[CurrentSystem]+1)*sizeof(int[4][3]));
      if((charinput1=='F')||(charinput1=='f'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][0][0]=CATION;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][0][1]=intinput1;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][1][0]=CATION;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][1][1]=intinput3;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][2][0]=CATION;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][2][1]=intinput5;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][3][0]=CATION;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][3][1]=intinput7;
      InversionBendDefinitions[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]][3][2]=intinput8;

      // optional reference dihedral angle
      InversionBendConstraintParameter[CurrentSystem]=(REAL*)realloc(InversionBendConstraintParameter[CurrentSystem],
                                                 (NumberOfInversionBendConstraints[CurrentSystem]+1)*sizeof(REAL));
      InversionBendConstraintParameter[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]]=UNDEFINED;
      if(sscanf(arguments,"%lf",&realinput1)==1)
        InversionBendConstraintParameter[CurrentSystem][NumberOfInversionBendConstraints[CurrentSystem]]=realinput1*DEG2RAD;

      NumberOfInversionBendConstraints[CurrentSystem]++;
    }

    if(strcasecmp("OutOfPlaneConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,arguments);
      sscanf(arguments,"%lf",&realinput1);
      OutOfPlaneDistanceDefinitions[CurrentSystem]=(int(*)[4][3])realloc(OutOfPlaneDistanceDefinitions[CurrentSystem],
                                         (NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]+1)*sizeof(int[4][3]));
      if((charinput1=='F')||(charinput1=='f'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][0][0]=CATION;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][0][1]=intinput1;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][1][0]=CATION;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][1][1]=intinput3;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][2][0]=CATION;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][2][1]=intinput5;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][3][0]=CATION;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][3][1]=intinput7;
      OutOfPlaneDistanceDefinitions[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]][3][2]=intinput8;

      // optional reference dihedral angle
      OutOfPlaneDistanceConstraintParameter[CurrentSystem]=(REAL*)realloc(OutOfPlaneDistanceConstraintParameter[CurrentSystem],
                                                 (NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]+1)*sizeof(REAL));
      OutOfPlaneDistanceConstraintParameter[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]]=UNDEFINED;
      if(sscanf(arguments,"%lf",&realinput1)==1)
        OutOfPlaneDistanceConstraintParameter[CurrentSystem][NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]]=realinput1*DEG2RAD;

      NumberOfOutOfPlaneDistanceConstraints[CurrentSystem]++;
    }


    if(strcasecmp("HarmonicDistanceConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,arguments);

      HarmonicDistanceDefinitions[CurrentSystem]=(int(*)[2][3])realloc(HarmonicDistanceDefinitions[CurrentSystem],
                                                 (NumberOfHarmonicDistanceConstraints[CurrentSystem]+1)*sizeof(int[2][3]));
      if((charinput1=='F')||(charinput1=='f'))
        HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0][0]=CATION;
      HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0][1]=intinput1;
      HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0][2]=intinput2;
      if((charinput2=='F')||(charinput2=='f'))
        HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][1][0]=CATION;
      HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][1][1]=intinput3;
      HarmonicDistanceDefinitions[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][1][2]=intinput4;

      // optional bond-constant and reference distance
      HarmonicDistanceConstraintParameters[CurrentSystem]=(REAL(*)[2])realloc(HarmonicDistanceConstraintParameters[CurrentSystem],
                                                          (NumberOfHarmonicDistanceConstraints[CurrentSystem]+1)*sizeof(REAL[2]));
      HarmonicDistanceConstraintParameters[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0]=1e6/ENERGY_TO_KELVIN;
      HarmonicDistanceConstraintParameters[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][1]=UNDEFINED;
      if(sscanf(arguments,"%lf%lf",&realinput1,&realinput2)==2)
      {
        HarmonicDistanceConstraintParameters[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0]=realinput1/ENERGY_TO_KELVIN;
        HarmonicDistanceConstraintParameters[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][1]=realinput2;
      }
      else if(sscanf(arguments,"%lf",&realinput1)==1)
        HarmonicDistanceConstraintParameters[CurrentSystem][NumberOfHarmonicDistanceConstraints[CurrentSystem]][0]=realinput1/ENERGY_TO_KELVIN;


      NumberOfHarmonicDistanceConstraints[CurrentSystem]++;
    }

    if(strcasecmp("HarmonicAngleConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,arguments);
      HarmonicAngleDefinitions[CurrentSystem]=(int(*)[3][3])realloc(HarmonicAngleDefinitions[CurrentSystem],
                                              (NumberOfHarmonicAngleConstraints[CurrentSystem]+1)*sizeof(int[3][3]));
      if((charinput1=='F')||(charinput1=='f'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0][0]=CATION;
      HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0][1]=intinput1;
      HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][1][0]=CATION;
      HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][1][1]=intinput3;
      HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][2][0]=CATION;
      HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][2][1]=intinput5;
      HarmonicAngleDefinitions[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][2][2]=intinput6;

      // optional bend-constant and reference angle
      HarmonicAngleConstraintParameters[CurrentSystem]=(REAL(*)[2])realloc(HarmonicAngleConstraintParameters[CurrentSystem],
                                                       (NumberOfHarmonicAngleConstraints[CurrentSystem]+1)*sizeof(REAL[2]));
      HarmonicAngleConstraintParameters[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0]=1e5/ENERGY_TO_KELVIN;
      HarmonicAngleConstraintParameters[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][1]=UNDEFINED;
      if(sscanf(arguments,"%lf%lf",&realinput1,&realinput2)==2)
      {
        HarmonicAngleConstraintParameters[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0]=realinput1/ENERGY_TO_KELVIN;
        HarmonicAngleConstraintParameters[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][1]=realinput2*DEG2RAD;
      }
      else if(sscanf(arguments,"%lf",&realinput1)==1)
        HarmonicAngleConstraintParameters[CurrentSystem][NumberOfHarmonicAngleConstraints[CurrentSystem]][0]=realinput1/ENERGY_TO_KELVIN;


      NumberOfHarmonicAngleConstraints[CurrentSystem]++;
    }

    if(strcasecmp("HarmonicDihedralConstraint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,arguments);
      HarmonicDihedralDefinitions[CurrentSystem]=(int(*)[4][3])realloc(HarmonicDihedralDefinitions[CurrentSystem],
                                      (NumberOfHarmonicDihedralConstraints[CurrentSystem]+1)*sizeof(int[4][3]));
      if((charinput1=='F')||(charinput1=='f'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0][0]=CATION;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0][1]=intinput1;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][1][0]=CATION;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][1][1]=intinput3;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][2][0]=CATION;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][2][1]=intinput5;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][3][0]=CATION;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][3][1]=intinput7;
      HarmonicDihedralDefinitions[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][3][2]=intinput8;

      // optional bend-constant and reference angle
      HarmonicDihedralConstraintParameters[CurrentSystem]=(REAL(*)[2])realloc(HarmonicDihedralConstraintParameters[CurrentSystem],
                                                          (NumberOfHarmonicDihedralConstraints[CurrentSystem]+1)*sizeof(REAL[2]));
      HarmonicDihedralConstraintParameters[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0]=1e5/ENERGY_TO_KELVIN;
      HarmonicDihedralConstraintParameters[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][1]=UNDEFINED;
      if(sscanf(arguments,"%lf%lf",&realinput1,&realinput2)==2)
      {
        HarmonicDihedralConstraintParameters[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0]=realinput1/ENERGY_TO_KELVIN;
        HarmonicDihedralConstraintParameters[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][1]=realinput2*DEG2RAD;
      }
      else if(sscanf(arguments,"%lf",&realinput1)==1)
        HarmonicDihedralConstraintParameters[CurrentSystem][NumberOfHarmonicDihedralConstraints[CurrentSystem]][0]=realinput1/ENERGY_TO_KELVIN;
      
      NumberOfHarmonicDihedralConstraints[CurrentSystem]++;
    }

    if(strcasecmp("DistanceConstraintType",keyword)==0) 
    {
      if(strcasecmp("Distance",firstargument)==0) DistanceConstraintType=DISTANCE_R;
      if(strcasecmp("DistanceSquared",firstargument)==0) DistanceConstraintType=DISTANCE_R_SQUARED;
    }
    if(strcasecmp("BendConstraintType",keyword)==0)
    {
      if(strcasecmp("Theta",firstargument)==0) BendConstraintType=THETA;
      if(strcasecmp("CosTheta",firstargument)==0) BendConstraintType=COS_THETA;
      if(strcasecmp("CosThetaSquared",firstargument)==0) BendConstraintType=COS_THETA_SQUARED;
    }
    if(strcasecmp("DihedralConstraintType",keyword)==0)
    {
      if(strcasecmp("Phi",firstargument)==0) DihedralConstraintType=PHI;
      if(strcasecmp("CosPhi",firstargument)==0) DihedralConstraintType=COS_PHI;
      if(strcasecmp("CosPhiSquared",firstargument)==0) DihedralConstraintType=COS_PHI_SQUARED;
    }
    if(strcasecmp("ImproperDihedralConstraintType",keyword)==0)
    {
      if(strcasecmp("Phi",firstargument)==0) ImproperDihedralConstraintType=PHI;
      if(strcasecmp("CosPhi",firstargument)==0) ImproperDihedralConstraintType=COS_PHI;
      if(strcasecmp("CosPhiSquared",firstargument)==0) ImproperDihedralConstraintType=COS_PHI_SQUARED;
    }
    if(strcasecmp("InversionBendConstraintType",keyword)==0)
    {
      if(strcasecmp("Chi",firstargument)==0) InversionBendConstraintType=CHI;
      if(strcasecmp("SinChi",firstargument)==0) InversionBendConstraintType=SIN_CHI;
      if(strcasecmp("SinChiSquared",firstargument)==0) InversionBendConstraintType=SIN_CHI_SQUARED;
    }
    if(strcasecmp("OutOfPlaneConstraintType",keyword)==0)
    {
      if(strcasecmp("Distance",firstargument)==0) OutOfPlaneConstraintType=DISTANCE_R;
      if(strcasecmp("DistanceSquared",firstargument)==0) OutOfPlaneConstraintType=DISTANCE_R_SQUARED;
    }
    if(strcasecmp("ComputeRattleSteps",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) ComputeRattleSteps=TRUE;
      if(strcasecmp("no",firstargument)==0) ComputeRattleSteps=FALSE;
    }

    if(strcasecmp("MeasureDihedralWithMidpoint",keyword)==0)
    {
      sscanf(arguments," %c %d %d %c %d %d %c %d %d %c %d %d %c %d %d %c %d %d%[^\n]",&charinput1,&intinput1,&intinput2,&charinput2,&intinput3,&intinput4,
                                                           &charinput3,&intinput5,&intinput6,&charinput4,&intinput7,&intinput8,
                                                           &charinput5,&intinput9,&intinput10,&charinput6,&intinput11,&intinput12,arguments);
      sscanf(arguments,"%lf",&realinput1);
      TwoPointDihedralDefinitions[CurrentSystem]=(int(*)[6][3])realloc(TwoPointDihedralDefinitions[CurrentSystem],
                                         (NumberOfTwoPointDihedralDefinitions[CurrentSystem]+1)*sizeof(int[6][3]));
      if((charinput1=='F')||(charinput1=='f'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][0][0]=FRAMEWORK;
      else if((charinput1=='A')||(charinput1=='a'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][0][0]=ADSORBATE;
      else if((charinput1=='C')||(charinput1=='c'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][0][0]=CATION;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][0][1]=intinput1;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][0][2]=intinput2;

      if((charinput2=='F')||(charinput2=='f'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][1][0]=FRAMEWORK;
      else if((charinput2=='A')||(charinput2=='a'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][1][0]=ADSORBATE;
      else if((charinput2=='C')||(charinput2=='c'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][1][0]=CATION;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][1][1]=intinput3;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][1][2]=intinput4;

      if((charinput3=='F')||(charinput3=='f'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][2][0]=FRAMEWORK;
      else if((charinput3=='A')||(charinput3=='a'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][2][0]=ADSORBATE;
      else if((charinput3=='C')||(charinput3=='c'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][2][0]=CATION;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][2][1]=intinput5;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][2][2]=intinput6;

      if((charinput4=='F')||(charinput4=='f'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][3][0]=FRAMEWORK;
      else if((charinput4=='A')||(charinput4=='a'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][3][0]=ADSORBATE;
      else if((charinput4=='C')||(charinput4=='c'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][3][0]=CATION;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][3][1]=intinput7;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][3][2]=intinput8;

      if((charinput5=='F')||(charinput5=='f'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][4][0]=FRAMEWORK;
      else if((charinput5=='A')||(charinput5=='a'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][4][0]=ADSORBATE;
      else if((charinput5=='C')||(charinput5=='c'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][4][0]=CATION;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][4][1]=intinput9;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][4][2]=intinput10;

      if((charinput6=='F')||(charinput6=='f'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][5][0]=FRAMEWORK;
      else if((charinput6=='A')||(charinput6=='a'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][5][0]=ADSORBATE;
      else if((charinput6=='C')||(charinput6=='c'))
        TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][5][0]=CATION;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][5][1]=intinput11;
      TwoPointDihedralDefinitions[CurrentSystem][NumberOfTwoPointDihedralDefinitions[CurrentSystem]][5][2]=intinput12;

      NumberOfTwoPointDihedralDefinitions[CurrentSystem]++;
    }


    // read Monte Carlo settings
    if(strcasecmp("CBMCBiasingMethod",keyword)==0) 
    {
      if(strcasecmp("LJ_Biasing",firstargument)==0) BiasingMethod=LJ_BIASING;
      if(strcasecmp("LJ_And_Real_Biasing",firstargument)==0) BiasingMethod=LJ_AND_REAL_BIASING;
    }
    if(strcasecmp("OptimizeAcceptenceEvery",keyword)==0) sscanf(arguments,"%d",&OptimizeAcceptenceEvery);
    if(strcasecmp("MinimumInnerCycles",keyword)==0) sscanf(arguments,"%d",&MinimumInnerCycles);
    if(strcasecmp("NumberOfTrialPositions",keyword)==0) sscanf(arguments,"%d",&NumberOfTrialPositions);
    if(strcasecmp("NumberOfTrialPositionsTorsion",keyword)==0) sscanf(arguments,"%d",&NumberOfTrialPositionsTorsion);
    if(strcasecmp("NumberOfTrialPositionsForTheFirstBead",keyword)==0) 
       sscanf(arguments,"%d",&NumberOfTrialPositionsForTheFirstBead);
    if(strcasecmp("NumberOfTrialMovesPerOpenBead",keyword)==0) sscanf(arguments,"%d",&NumberOfTrialMovesPerOpenBead);
    if(strcasecmp("TargetAccRatioSmallMCScheme",keyword)==0) sscanf(arguments,"%lf",&TargetAccRatioSmallMCScheme);
    if(strcasecmp("TargetAccRatioTranslation",keyword)==0) sscanf(arguments,"%lf",&TargetAccRatioTranslation);
    if(strcasecmp("EnergyOverlapCriteria",keyword)==0) sscanf(arguments,"%lf",&EnergyOverlapCriteria);
    if(strcasecmp("MinimumRosenbluthFactor",keyword)==0) sscanf(arguments,"%lf",&MinimumRosenbluthFactor);
	if(strcasecmp("EBCBMC",keyword)==0)
	{
		if(strcasecmp("yes",firstargument)==0) EBCBMC=TRUE;
		if(strcasecmp("no",firstargument)==0) EBCBMC=FALSE;
	}
	if(strcasecmp("UseFirstBeadEnergy",keyword)==0)
	  {
		  if(strcasecmp("yes",firstargument)==0) UseEBCBMCEnergyFirstBead=TRUE;
		  if(strcasecmp("no",firstargument)==0) UseEBCBMCEnergyFirstBead=FALSE;
	  }	  
	if(strcasecmp("HermiteInterpolation",keyword)==0)
	{
		  if(strcasecmp("yes",firstargument)==0) HermiteInterpolation=TRUE;
		  if(strcasecmp("no",firstargument)==0) HermiteInterpolation=FALSE;
	}

    // transition state theory settings
    if(strcasecmp("FreeEnergyMappingType",keyword)==0) 
    {
      if(strcasecmp("A",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=A_MAPPING;
      if(strcasecmp("B",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=B_MAPPING;
      if(strcasecmp("C",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=C_MAPPING;
      if(strcasecmp("ABC",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=ABC_MAPPING;

      if(strcasecmp("AB_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_AB_DIAGONAL;
      if(strcasecmp("AC_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_AC_DIAGONAL;
      if(strcasecmp("BC_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_BC_DIAGONAL;

      if(strcasecmp("O_AB_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_O_AB_DIAGONAL;
      if(strcasecmp("O_AC_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_O_AC_DIAGONAL;
      if(strcasecmp("O_BC_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_O_BC_DIAGONAL;

      if(strcasecmp("A_BC_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_A_BC_DIAGONAL;
      if(strcasecmp("B_AC_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_B_AC_DIAGONAL;
      if(strcasecmp("C_AB_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_C_AB_DIAGONAL;
      if(strcasecmp("O_ABC_DIAGONAL",firstargument)==0) FreeEnergyMappingType[CurrentSystem]=MAP_O_ABC_DIAGONAL;
    }
    if(strcasecmp("PutMoleculeOnBarrier",keyword)==0) 
    {
      if(strcasecmp("yes",firstargument)==0) PutMoleculeOnBarrier[CurrentSystem]=TRUE;
      if(strcasecmp("no",firstargument)==0) PutMoleculeOnBarrier[CurrentSystem]=FALSE;
    }
    if(strcasecmp("BarrierPosition",keyword)==0) sscanf(arguments,"%lf %lf %lf",
       &BarrierPosition[CurrentSystem].x,&BarrierPosition[CurrentSystem].y,&BarrierPosition[CurrentSystem].z);
    if(strcasecmp("BarrierNormal",keyword)==0) sscanf(arguments,"%lf %lf %lf",
       &BarrierNormal[CurrentSystem].x,&BarrierNormal[CurrentSystem].y,&BarrierNormal[CurrentSystem].z);
    if(strcasecmp("BarrierAngle",keyword)==0) 
      if(sscanf(arguments,"%lf %lf %lf",&BarrierAngle[CurrentSystem].x,&BarrierAngle[CurrentSystem].y,&BarrierAngle[CurrentSystem].z))
         {BarrierAngle[CurrentSystem].x*=DEG2RAD; BarrierAngle[CurrentSystem].y*=DEG2RAD; BarrierAngle[CurrentSystem].z*=DEG2RAD;};
    if(strcasecmp("MaxBarrierDistance",keyword)==0) sscanf(arguments,"%lf",&MaxBarrierDistance[CurrentSystem]);
    if(strcasecmp("MaxBarrierTime",keyword)==0) sscanf(arguments,"%lf",&MaxBarrierTime[CurrentSystem]);
    if(strcasecmp("NumberOfVelocities",keyword)==0) sscanf(arguments,"%d",&NumberOfVelocities[CurrentSystem]);

    if(strcasecmp("PrintFrameworkBondStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkBondStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkBondStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkUreyBradleyStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkUreyBradleyStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkUreyBradleyStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkBendStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkInversionBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkInversionBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkInversionBendStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkImproperTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkImproperTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkImproperTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkBondBondStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkBondBondStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkBondBondStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkBondBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkBondBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkBondBendStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkBendBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkBendBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkBendBendStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkBendTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkBendTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkBendTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkIntraVDWStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkIntraVDWStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkIntraVDWStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkIntraChargeChargeStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkIntraChargeChargeStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkIntraChargeChargeStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkIntraChargeBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkIntraChargeBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkIntraChargeBondDipoleStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkIntraBondDipoleBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkIntraBondDipoleBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkIntraBondDipoleBondDipoleStatus=FALSE;
    }

    if(strcasecmp("PrintAdsorbateBondStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateBondStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateBondStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateUreyBradleyStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateUreyBradleyStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateUreyBradleyStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateBendStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateImproperTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateImproperTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateImproperTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateBondBondStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateBondBondStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateBondBondStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateBondBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateBondBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateBondBendStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateBendBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateBendBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateBendBendStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateBondTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateBondTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateBondTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateBendTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateBendTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateBendTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateIntraVDWStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateIntraVDWStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateIntraVDWStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateIntraChargeChargeStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateIntraChargeChargeStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateIntraChargeChargeStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateIntraChargeBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateIntraChargeBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateIntraChargeBondDipoleStatus=FALSE;
    }
    if(strcasecmp("PrintAdsorbateIntraBondDipoleBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintAdsorbateIntraBondDipoleBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintAdsorbateIntraBondDipoleBondDipoleStatus=FALSE;
    }

    if(strcasecmp("PrintCationBondStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationBondStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationBondStatus=FALSE;
    }
    if(strcasecmp("PrintCationUreyBradleyStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationUreyBradleyStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationUreyBradleyStatus=FALSE;
    }
    if(strcasecmp("PrintCationBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationBendStatus=FALSE;
    }
    if(strcasecmp("PrintCationTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintCationImproperTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationImproperTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationImproperTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintCationBondBondStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationBondBondStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationBondBondStatus=FALSE;
    }
    if(strcasecmp("PrintCationBondBendStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationBondBendStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationBondBendStatus=FALSE;
    }
    if(strcasecmp("PrintCationBondTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationBondTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationBondTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintCationBendTorsionStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationBendTorsionStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationBendTorsionStatus=FALSE;
    }
    if(strcasecmp("PrintCationIntraVDWStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationIntraVDWStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationIntraVDWStatus=FALSE;
    }
    if(strcasecmp("PrintCationIntraChargeChargeStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationIntraChargeChargeStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationIntraChargeChargeStatus=FALSE;
    }
    if(strcasecmp("PrintCationIntraChargeBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationIntraChargeBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationIntraChargeBondDipoleStatus=FALSE;
    }
    if(strcasecmp("PrintCationIntraBondDipoleBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintCationIntraBondDipoleBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintCationIntraBondDipoleBondDipoleStatus=FALSE;
    }

    if(strcasecmp("PrintInterVDWStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintInterVDWStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintInterVDWStatus=FALSE;
    }
    if(strcasecmp("PrintInterChargeChargeStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintInterChargeChargeStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintInterChargeChargeStatus=FALSE;
    }
    if(strcasecmp("PrintInterChargeBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintInterChargeBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintInterChargeBondDipoleStatus=FALSE;
    }
    if(strcasecmp("PrintInterBondDipoleBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintInterBondDipoleBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintInterBondDipoleBondDipoleStatus=FALSE;
    }

    if(strcasecmp("PrintFrameworkAdsorbateVDWStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkAdsorbateVDWStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkAdsorbateVDWStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkAdsorbateChargeChargeStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkAdsorbateChargeChargeStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkAdsorbateChargeChargeStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkAdsorbateChargeBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkAdsorbateChargeBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkAdsorbateChargeBondDipoleStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkAdsorbateBondDipoleBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkAdsorbateBondDipoleBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkAdsorbateBondDipoleBondDipoleStatus=FALSE;
    }

    if(strcasecmp("PrintFrameworkCationVDWStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkCationVDWStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkCationVDWStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkCationChargeChargeStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkCationChargeChargeStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkCationChargeChargeStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkCationChargeBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkCationChargeBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkCationChargeBondDipoleStatus=FALSE;
    }
    if(strcasecmp("PrintFrameworkCationBondDipoleBondDipoleStatus",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) PrintFrameworkCationBondDipoleBondDipoleStatus=TRUE;
      if(strcasecmp("no",firstargument)==0) PrintFrameworkCationBondDipoleBondDipoleStatus=FALSE;
    }



    if(strcasecmp("PrintFrameworkBondStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkBondStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkUreyBradleyStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkUreyBradleyStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkBendStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkInversionBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkInversionBendStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkImproperTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkImproperTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkBondBondStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkBondBondStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkBondBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkBondBendStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkBendBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkBendBendStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkBendTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkBendTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkIntraVDWStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkIntraVDWStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkIntraChargeChargeStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkIntraChargeChargeStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkIntraChargeBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkIntraChargeBondDipoleStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkIntraBondDipoleBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkIntraBondDipoleBondDipoleStatus=TRUE;
    }

    if(strcasecmp("PrintAdsorbateBondStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateBondStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateUreyBradleyStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateUreyBradleyStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateBendStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateImproperTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateImproperTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateBondBondStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateBondBondStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateBondBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateBondBendStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateBendBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateBendBendStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateBondTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateBondTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateBendTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateBendTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateIntraVDWStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateIntraVDWStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateIntraChargeChargeStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateIntraChargeChargeStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateIntraChargeBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateIntraChargeBondDipoleStatus=TRUE;
    }
    if(strcasecmp("PrintAdsorbateIntraBondDipoleBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintAdsorbateIntraBondDipoleBondDipoleStatus=TRUE;
    }

    if(strcasecmp("PrintCationBondStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationBondStatus=TRUE;
    }
    if(strcasecmp("PrintCationUreyBradleyStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationUreyBradleyStatus=TRUE;
    }
    if(strcasecmp("PrintCationBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationBendStatus=TRUE;
    }
    if(strcasecmp("PrintCationTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintCationImproperTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationImproperTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintCationBondBondStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationBondBondStatus=TRUE;
    }
    if(strcasecmp("PrintCationBondBendStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationBondBendStatus=TRUE;
    }
    if(strcasecmp("PrintCationBondTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationBondTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintCationBendTorsionStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationBendTorsionStatus=TRUE;
    }
    if(strcasecmp("PrintCationIntraVDWStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationIntraVDWStatus=TRUE;
    }
    if(strcasecmp("PrintCationIntraChargeChargeStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationIntraChargeChargeStatus=TRUE;
    }
    if(strcasecmp("PrintCationIntraChargeBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationIntraChargeBondDipoleStatus=TRUE;
    }
    if(strcasecmp("PrintCationIntraBondDipoleBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintCationIntraBondDipoleBondDipoleStatus=TRUE;
    }

    if(strcasecmp("PrintInterVDWStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintInterVDWStatus=TRUE;
    }
    if(strcasecmp("PrintInterChargeChargeStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintInterChargeChargeStatus=TRUE;
    }
    if(strcasecmp("PrintInterChargeBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintInterChargeBondDipoleStatus=TRUE;
    }
    if(strcasecmp("PrintInterBondDipoleBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintInterBondDipoleBondDipoleStatus=TRUE;
    }

    if(strcasecmp("PrintFrameworkAdsorbateVDWStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkAdsorbateVDWStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkAdsorbateChargeChargeStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkAdsorbateChargeChargeStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkAdsorbateChargeBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkAdsorbateChargeBondDipoleStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkAdsorbateBondDipoleBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkAdsorbateBondDipoleBondDipoleStatus=TRUE;
    }

    if(strcasecmp("PrintFrameworkCationVDWStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkCationVDWStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkCationChargeChargeStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkCationChargeChargeStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkCationChargeBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkCationChargeBondDipoleStatus=TRUE;
    }
    if(strcasecmp("PrintFrameworkCationBondDipoleBondDipoleStatusOnly",keyword)==0)
    {
      SetPrintStatusToFalse(); 
      PrintFrameworkCationBondDipoleBondDipoleStatus=TRUE;
    }


    // create frameworks
    if(strcasecmp("FrameworkProbability",keyword)==0) sscanf(arguments,"%lf",&Framework[CurrentSystem].FrameworkProbability[CurrentFramework]);
    if(strcasecmp("FrameworkExclusion",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].FrameworkExclusion=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].FrameworkExclusion=FALSE;
    }
    if(strcasecmp("RemoveHydrogenDisorder",keyword)==0)
    {
      if(strcasecmp("yes",firstargument)==0) Framework[CurrentSystem].RemoveHydrogenDisorder=TRUE;
      if(strcasecmp("no",firstargument)==0) Framework[CurrentSystem].RemoveHydrogenDisorder=FALSE;
    }
  }
  fclose(FilePtr);

  // if the frameworks can be moved apart, the system is set to 'flexible'
  // the individual frameworks can be either 'fleixble' or 'rigid'
  if(ProbabilityFrameworkShiftMove>0.0)
  {
    for(i=0;i<NumberOfSystems;i++)
      if(Framework[i].NumberOfFrameworks>1) Framework[i].FrameworkModel=FLEXIBLE;
  }


  // get the random number seed from the system-time if not specified
  if(seed<=0lu) seed=time(0lu);

  // Initialize random-number generator using the "random" seed
  InitializeRandomNumberGenerator(seed);

  if(BoundaryCondition[CurrentSystem]==FINITE)
  {
    CutOffVDW=DBL_MAX;
    CutOffVDWSquared=DBL_MAX;
    CutOffChargeCharge=DBL_MAX;
    CutOffChargeChargeSquared=DBL_MAX;
    CutOffChargeBondDipole=DBL_MAX;
    CutOffChargeBondDipoleSquared=DBL_MAX;
    CutOffBondDipoleBondDipole=DBL_MAX;
    CutOffBondDipoleBondDipoleSquared=DBL_MAX;
    CutOffVDWSwitch=DBL_MAX;
    CutOffVDWSwitchSquared=DBL_MAX;
    CutOffChargeChargeSwitch=DBL_MAX;
    CutOffChargeChargeSwitchSquared=DBL_MAX;
    CutOffChargeBondDipoleSwitch=DBL_MAX;
    CutOffChargeBondDipoleSwitchSquared=DBL_MAX;
    CutOffBondDipoleBondDipoleSwitch=DBL_MAX;
    CutOffBondDipoleBondDipoleSwitchSquared=DBL_MAX;
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
    {
      printf("ERROR: you can not use a finite system *and* the Ewald summation\n");
      exit(0);
    }
  }
  else
  {
    // Calculate cutoff-related values
    CutOffVDWSquared=SQR(CutOffVDW);
    CutOffChargeChargeSquared=SQR(CutOffChargeCharge);
    CutOffChargeBondDipoleSquared=SQR(CutOffChargeBondDipole);
    CutOffBondDipoleBondDipoleSquared=SQR(CutOffBondDipoleBondDipole);

    InverseCutOffVDW=1.0/CutOffVDW;
    InverseCutOffChargeCharge=1.0/CutOffChargeCharge;

    if(CutOffVDWSwitch<=0)
      CutOffVDWSwitch=0.9*CutOffVDW;
    CutOffVDWSwitchSquared=SQR(CutOffVDWSwitch);

    if(CutOffChargeChargeSwitch<=0.0)
      CutOffChargeChargeSwitch=0.65*CutOffChargeCharge;
    CutOffChargeChargeSwitchSquared=SQR(CutOffChargeChargeSwitch);

    if(CutOffChargeBondDipoleSwitch<=0.0)
      CutOffChargeBondDipoleSwitch=0.7*CutOffChargeBondDipole;
    CutOffChargeBondDipoleSwitchSquared=SQR(CutOffChargeBondDipoleSwitch);

    if(CutOffBondDipoleBondDipoleSwitch<=0.0)
      CutOffBondDipoleBondDipoleSwitch=0.75*CutOffBondDipoleBondDipole;
    CutOffBondDipoleBondDipoleSwitchSquared=SQR(CutOffBondDipoleBondDipoleSwitch);
  }

  ComputeSwitchingFactors();

  if(UseTabularGrid)
  {
    for(i=0;i<NumberOfSystems;i++)
      Framework[i].FrameworkModel=GRID,printf("Setting GRID\n");
  }

  // read pseudo-atoms definitions
  // Note these can be extended by atoms read from cif-files
  ReadPseudoAtomsDefinitions();
  //ReadForceFieldDefinitionsMixingRules();
  //ReadForceFieldDefinitions();

  // PrintEvery should be larger then zero
  if(PrintEvery<1) PrintEvery=5000;

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      if(read_frameworks[CurrentSystem])
      {
        for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
        {
          switch(Framework[CurrentSystem].InputFileType[CurrentFramework])
          {
            case CIF_FORMAT:
              ReadFrameworkDefinitionCIF();
              break;
            case CSSR_FORMAT:
              ReadFrameworkDefinitionCSSR();
              break;
            case MOL_FORMAT:
              ReadFrameworkDefinitionMOL();
              break;
            default:
              printf("Unknown structure-format\n");
              break;
          }
          if(Framework[CurrentSystem].CalculateSpaceGroup[CurrentFramework]) DetermineSpaceGroup();
        }

        InitializeReplicaBox();

        // check modification rules
        for(i=0;i<NumberOfModificationRules;i++)
        {
          // convert strings to ids
          for(k=0;k<6;k++)
          {
            Type=-1;
            for(j=0;j<NumberOfPseudoAtoms;j++)
              if(strcasecmp(PseudoAtoms[j].Name,ModifyFrameworkAtoms[i][k])==0) Type=j;
            ModifyFrameworkAtomTypes[i][k]=Type;
          }

          // check atom strings
          switch(ModificationRuleType[i])
          {
            case MODIFY_FRAMEWORKATOM_CONNECTED_TO:
              for(k=0;k<2;k++)
              {
                if(ModifyFrameworkAtomTypes[i][k]<0)
                {
                  printf("Unknown framework atom '%s' in modification rule %d atom %d\n",ModifyFrameworkAtoms[i][k],i,k);
                  exit(0);
                }
              }
              break;
            case MODIFY_FRAMEWORKATOM_DIMER:
              for(k=0;k<4;k++)
              {
                if(ModifyFrameworkAtomTypes[i][k]<0)
                {
                  printf("Unknown framework atom '%s' in modification rule %d atom %d\n",ModifyFrameworkAtoms[i][k],i,k);
                  exit(0);
                }
              }
              break;
            case MODIFY_FRAMEWORKATOM_TRIPLE:
              for(k=0;k<6;k++)
              {
                if(ModifyFrameworkAtomTypes[i][k]<0)
                {
                  printf("Unknown framework atom '%s' in modification rule %d atom %d\n",ModifyFrameworkAtoms[i][k],i,k);
                  exit(0);
                }
              }
              break;
          }
        }


        // check forbidden connectivity rules
        for(i=0;i<NumberOfForbiddenConnectivityRules;i++)
        {
          // convert strings to ids
          for(k=0;k<3;k++)
          {
            Type=-1;
            for(j=0;j<NumberOfPseudoAtoms;j++)
              if(strcasecmp(PseudoAtoms[j].Name,ForbiddenConnectivityAtoms[i][k])==0) Type=j;
            ForbiddenConnectivityTypes[i][k]=Type;
          }
          if(ForbiddenConnectivityTypes[i][0]<0)
          {
            printf("Unknown framework atom %s in forbidden connectivity rule %d\n",ForbiddenConnectivityAtoms[i][0],0);
            exit(0);
          }
          if(ForbiddenConnectivityTypes[i][1]<0)
          {
            printf("Unknown framework atom %s in forbidden connectivity rule %d\n",ForbiddenConnectivityAtoms[i][1],1);
            exit(0);
          }
          if(ForbiddenConnectivityTypes[i][2]<0)
          {
            printf("Unknown framework atom %s in forbidden connectivity rule %d\n",ForbiddenConnectivityAtoms[i][2],2);
            exit(0);
          }
        }


        // make a connectivity list to look for neighbors
        // the list is used to obtain bonds, bends, torsion etc 

        // read the flexible framework model definitions
        ReadFrameworkDefinition();

        MakeExclusionMatrix(CurrentSystem);
        MakeExcludedInteractionLists(CurrentSystem);

        ReadBlockingPockets();
      }

      // measure the original shift between the frameworks
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        CurrentFramework=f1;
        posf1=GetFrameworkCenterOfMassPosition();
        for(f2=f1+1;f2<Framework[CurrentSystem].NumberOfFrameworks;f2++)
        {
          CurrentFramework=f2;
          posf2=GetFrameworkCenterOfMassPosition();
          dr.x=posf1.x-posf2.x;
          dr.y=posf1.y-posf2.y;
          dr.z=posf1.z-posf2.z;
          dr=ApplyBoundaryCondition(dr);
          r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
          OriginalFrameworkShiftDirAvg[CurrentSystem][f1][f2]=r;
          OriginalFrameworkShift[CurrentSystem][f1][f2].x=dr.x;
          OriginalFrameworkShift[CurrentSystem][f1][f2].y=dr.y;
          OriginalFrameworkShift[CurrentSystem][f1][f2].z=dr.z;
        }
      }
      if(strlen(Framework[CurrentSystem].NameIons)>0)
        ReadIonSitingDefinition();
    }
  }

  // read forcefields after the structure has been read
  // (because CIF-files can add atom-types)
  ReadForceFieldDefinitionsMixingRules();
  ReadForceFieldDefinitions();
  ComputeDummyInteractions();
  ComputePotentialShifts();

  if(ChargeFromChargeEquilibration)
  {
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      ChargeEquilibration();
  }
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    CheckFrameworkCharges();

  for(i=0;i<NumberOfPseudoAtomsWithoutVDWInteraction;i++)
    PseudoAtoms[ReturnPseudoAtomNumber(PseudoAtomsWithoutVDWInteraction[i])].HasVDWInteraction=FALSE;


  // Read the definitions of the components
  MaxNumberOfCoulombicSites=0;
  LargestNumberOfCoulombicSites=0;
  MaxNumberOfBondDipoleSites=0;
  LargestNumberOfBondDipoleSites=0;
  for(i=0;i<NumberOfSystems;i++)
  {
    MaxNumberOfAdsorbateMolecules[i]=0;
    MaxNumberOfCationMolecules[i]=0;
    NumberOfAtomsPerSystem[i]=0;
    NumberOfChargesPerSystem[i]=0;
    NumberOfBondDipolesPerSystem[i]=0;
  }

  for(i=0;i<NumberOfComponents;i++)
  {
    ReadComponentDefinition(i);
    if(Components[i].StartingBead>=Components[i].NumberOfAtoms)
    {
      printf("Starting bead (%d) not within number of atoms: %d\n",
        Components[i].StartingBead,Components[i].NumberOfAtoms);
      exit(0);
    }
    for(j=0;j<NumberOfSystems;j++)
    {
      if(Components[i].ExtraFrameworkMolecule)
        MaxNumberOfCationMolecules[j]+=Components[i].CreateNumberOfMolecules[j];
      else
        MaxNumberOfAdsorbateMolecules[j]+=Components[i].CreateNumberOfMolecules[j];

      NumberOfAtomsPerSystem[j]+=Components[i].NumberOfAtoms*Components[i].CreateNumberOfMolecules[j];
      NumberOfChargesPerSystem[j]+=Components[i].NumberOfCharges*Components[i].CreateNumberOfMolecules[j];
      NumberOfBondDipolesPerSystem[j]+=Components[i].NumberOfBondDipoles*Components[i].CreateNumberOfMolecules[j];
    }
  }

  // update the maximum amount of atoms (maximum over all systems)
  LargestNumberOfCoulombicSites=0;
  for(i=0;i<NumberOfSystems;i++)
  {
    nr_sites=NumberOfChargesPerSystem[i];
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      nr_sites+=Framework[i].NumberOfCharges[j];
    if(nr_sites>LargestNumberOfCoulombicSites) LargestNumberOfCoulombicSites=nr_sites;

    NumberOfAtomsPerSystem[i]=0;
    NumberOfChargesPerSystem[i]=0;
  }

  LargestNumberOfBondDipoleSites=0;
  for(i=0;i<NumberOfSystems;i++)
  {
    nr_sites=NumberOfBondDipolesPerSystem[i];
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      nr_sites+=Framework[i].NumberOfBondDipoles[j];
    if(nr_sites>LargestNumberOfBondDipoleSites) LargestNumberOfBondDipoleSites=nr_sites;

    NumberOfBondDipolesPerSystem[i]=0;
  }


  // the components are read, and all limits are known
  // allocate the memory for use in the cbmc routines
  AllocateCBMCMemory();
  AllocateMCMovesMemory();
  AllocateGridMemory();

  for(i=0;i<NumberOfGrids;i++)
    GridTypeList[i]=ReturnPseudoAtomNumber(GridTypeListName[i]);

  for(i=0;i<NumberOfSystems;i++)
  {
    MaximumVolumeChange[i]=0.025;
    MaximumBoxShapeChange[i].ax=0.1;
    MaximumBoxShapeChange[i].ay=0.1;
    MaximumBoxShapeChange[i].az=0.1;
    MaximumBoxShapeChange[i].bx=0.1;
    MaximumBoxShapeChange[i].by=0.1;
    MaximumBoxShapeChange[i].bz=0.1;
    MaximumBoxShapeChange[i].cx=0.1;
    MaximumBoxShapeChange[i].cy=0.1;
    MaximumBoxShapeChange[i].cz=0.1;
    MaximumGibbsVolumeChange[i]=0.1;

    for(j=0;j<NumberOfComponents;j++)
    {
      switch(Components[j].TranslationDirection)
      {
        case XYZ_DIR:
        case XY_DIR:
        case XZ_DIR:
        case YZ_DIR:
        case X_DIR:
        case Y_DIR:
        case Z_DIR:
          MaximumTranslation[i][j].x=1.3;
          MaximumTranslation[i][j].y=1.3;
          MaximumTranslation[i][j].z=1.3;
          break;
        case ABC_DIR:
        case AB_DIR:
        case AC_DIR:
        case BC_DIR:
        case A_DIR:
        case B_DIR:
        case C_DIR:
        default:
          MaximumTranslation[i][j].x=0.03;
          MaximumTranslation[i][j].y=0.03;
          MaximumTranslation[i][j].z=0.03;
          break;
      }
      MaximumTranslationInPlane[i][j].x=0.3;
      MaximumTranslationInPlane[i][j].y=0.3;
      MaximumTranslationInPlane[i][j].z=0.3;
    }
  }


  // components have been read and the maximum amount of molecules has been determined for NVT
  if(!Swapable)
  {
    // the +1 molecule for the Widom particle insertion particle
    MaxNumberOfCoulombicSites=LargestNumberOfCoulombicSites+2*MaxNumberOfBeads;
    MaxNumberOfBondDipoleSites=LargestNumberOfBondDipoleSites+2*MaxNumberOfBeads;

    Adsorbates=(ADSORBATE_MOLECULE**)calloc(NumberOfSystems,sizeof(ADSORBATE_MOLECULE*));
    for(i=0;i<NumberOfSystems;i++)
    {
      if(MaxNumberOfAdsorbateMolecules[i]>0)
        Adsorbates[i]=(ADSORBATE_MOLECULE*)calloc(MaxNumberOfAdsorbateMolecules[i],sizeof(ADSORBATE_MOLECULE));
    }

    Cations=(CATION_MOLECULE**)calloc(NumberOfSystems,sizeof(CATION_MOLECULE*));
    for(i=0;i<NumberOfSystems;i++)
    {
      if(MaxNumberOfCationMolecules[i]>0)
        Cations[i]=(CATION_MOLECULE*)calloc(MaxNumberOfCationMolecules[i],sizeof(CATION_MOLECULE));
    }
  }
  else
  {
    // make sure these are off for a grand-canonical system
    for(i=0;i<NumberOfSystems;i++)
    {
      ComputeInfraRedSpectra[i]=FALSE;
      ComputeMSDOrderN[i]=FALSE;
      ComputeVACFOrderN[i]=FALSE;
      ComputeRVACFOrderN[i]=FALSE;
      ComputeMSD[i]=FALSE;
      ComputeVACF[i]=FALSE;
    }

    // default: start with 512 atoms in addition to a possible framework
    MaxNumberOfCoulombicSites=LargestNumberOfCoulombicSites+MAX2(MaxNumberOfBeads,512);
    MaxNumberOfBondDipoleSites=LargestNumberOfBondDipoleSites+MAX2(MaxNumberOfBeads,512);

    // default: start with 256 adsorbates
    Adsorbates=(ADSORBATE_MOLECULE**)calloc(NumberOfSystems,sizeof(ADSORBATE_MOLECULE*));
    for(i=0;i<NumberOfSystems;i++)
    {
      MaxNumberOfAdsorbateMolecules[i]=256;
      Adsorbates[i]=(ADSORBATE_MOLECULE*)calloc(MaxNumberOfAdsorbateMolecules[i],sizeof(ADSORBATE_MOLECULE));
    }

    // default: start with 1 cation (zero can not be reallocated later)
    Cations=(CATION_MOLECULE**)calloc(NumberOfSystems,sizeof(CATION_MOLECULE*));
    for(i=0;i<NumberOfSystems;i++)
    {
      MaxNumberOfCationMolecules[i]=1;
      Cations[i]=(CATION_MOLECULE*)calloc(MaxNumberOfCationMolecules[i],sizeof(CATION_MOLECULE));
    }
  }

  // HERE

  // should be improved
  BondTypes[RIGID_BOND].nr_args=1;

  // initialize the box for "Box" and "Angles", both are specified independently and in order to allow the
  // statements to be specified in any order the construction is done after both are known
  for(i=0;i<NumberOfSystems;i++)
  {
    if(InitializeBox[i])
    {

      switch(Dimension)
      {
        case 2:
          // Construct transformation matrix Box to go from abc-coordinates to xyz-coordinates
          A=UnitCellSize[i].x;
          B=UnitCellSize[i].y;
          C=0.0;
          UnitCellBox[i].ax=A;   UnitCellBox[i].bx=B*cos(GammaAngle[i]); UnitCellBox[i].cx=0.0;
          UnitCellBox[i].ay=0.0; UnitCellBox[i].by=B*sin(GammaAngle[i]); UnitCellBox[i].cy=0.0;
          UnitCellBox[i].az=0.0; UnitCellBox[i].bz=0.0;                  UnitCellBox[i].cz=0.0;

          tempd=1.0/(UnitCellBox[i].ax*UnitCellBox[i].by-B*cos(GammaAngle[i]));
          InverseUnitCellBox[i].ax=tempd*UnitCellBox[i].by;  InverseUnitCellBox[i].bx=-tempd*UnitCellBox[i].bx; InverseUnitCellBox[i].cx=0.0;
          InverseUnitCellBox[i].ay=-tempd*UnitCellBox[i].ay; InverseUnitCellBox[i].by=tempd*UnitCellBox[i].ax;  InverseUnitCellBox[i].cy=0.0;
          InverseUnitCellBox[i].az=0.0;                      InverseUnitCellBox[i].cy=0.0;                      InverseUnitCellBox[i].cz=0.0;

          A=(REAL)NumberOfUnitCells[i].x*UnitCellSize[i].x;
          B=(REAL)NumberOfUnitCells[i].y*UnitCellSize[i].y;
          C=0.0;
          Box[i].ax=A;   Box[i].bx=B*cos(GammaAngle[i]); Box[i].cx=0.0;
          Box[i].ay=0.0; Box[i].by=B*sin(GammaAngle[i]); Box[i].cy=0.0;
          Box[i].az=0.0; Box[i].bz=0.0;                  Box[i].cz=0.0;

          tempd=1.0/(Box[i].ax*Box[i].by-B*cos(GammaAngle[i]));
          InverseBox[i].ax=tempd*Box[i].by;  InverseBox[i].bx=-tempd*Box[i].bx; InverseBox[i].cx=0.0;
          InverseBox[i].ay=-tempd*Box[i].ay; InverseBox[i].by=tempd*Box[i].ax;  InverseBox[i].cy=0.0;
          InverseBox[i].az=0.0;              InverseBox[i].cy=0.0;              InverseBox[i].cz=0.0;

          // Calculate box-properties
          CellProperties(&Box[i],&BoxProperties[i],&Volume[i]);

          // Compute the invers box and properties of the inverse box
          CellProperties(&InverseBox[i],&InverseBoxProperties[i],&det);
          break;
        case 3:
          // Construct transformation matrix Box to go from abc-coordinates to xyz-coordinates
          A=UnitCellSize[i].x;
          B=UnitCellSize[i].y;
          C=UnitCellSize[i].z;
          tempd=(cos(AlphaAngle[i])-cos(GammaAngle[i])*cos(BetaAngle[i]))/sin(GammaAngle[i]);
          UnitCellBox[i].ax=A;   UnitCellBox[i].bx=B*cos(GammaAngle[i]); UnitCellBox[i].cx=C*cos(BetaAngle[i]);
          UnitCellBox[i].ay=0.0; UnitCellBox[i].by=B*sin(GammaAngle[i]); UnitCellBox[i].cy=C*tempd;
          UnitCellBox[i].az=0.0; UnitCellBox[i].bz=0.0;               
          UnitCellBox[i].cz=C*sqrt(1.0-SQR(cos(BetaAngle[i]))-SQR(tempd));
          Invert3x3Matrix(&UnitCellBox[i],&InverseUnitCellBox[i],&det);

          A=(REAL)NumberOfUnitCells[i].x*UnitCellSize[i].x;
          B=(REAL)NumberOfUnitCells[i].y*UnitCellSize[i].y;
          C=(REAL)NumberOfUnitCells[i].z*UnitCellSize[i].z;
          tempd=(cos(AlphaAngle[i])-cos(GammaAngle[i])*cos(BetaAngle[i]))/sin(GammaAngle[i]);
          Box[i].ax=A;   Box[i].bx=B*cos(GammaAngle[i]); Box[i].cx=C*cos(BetaAngle[i]);
          Box[i].ay=0.0; Box[i].by=B*sin(GammaAngle[i]); Box[i].cy=C*tempd;
          Box[i].az=0.0; Box[i].bz=0.0;                  Box[i].cz=C*sqrt(1.0-SQR(cos(BetaAngle[i]))-SQR(tempd));

          // Calculate box-properties
          CellProperties(&Box[i],&BoxProperties[i],&Volume[i]);

          // Compute the invers box and properties of the inverse box
          Invert3x3Matrix(&Box[i],&InverseBox[i],&det);
          CellProperties(&InverseBox[i],&InverseBoxProperties[i],&det);
          break;
      }
    }

    // determine boundary conditions from angles
    if(BoundaryCondition[i]==UNINITIALIZED_BOUNDARY_CONDITION)
    {
      if((fabs(AlphaAngle[i]-90.0)>0.001)||(fabs(BetaAngle[i]-90.0)>0.001)||(fabs(GammaAngle[i]-90.0)>0.001))
        BoundaryCondition[i]=TRICLINIC;
      else BoundaryCondition[i]=RECTANGULAR;
    }
  }

  // the probabilities are just relative with respect to each other
  // so rescale and clean up all the specified probabilities after they are known
  RescaleComponentProbabilities();
  RescaleMolarFractions();

  // calculate box-properties
  for(i=0;i<NumberOfSystems;i++)
  {
    CellProperties(&Box[i],&BoxProperties[i],&Volume[i]);
    InverseBoxMatrix(&Box[i],&InverseBox[i]);
    CellProperties(&InverseBox[i],&InverseBoxProperties[i],&det);

    CurrentSystem=i;
    InitializeReplicaBox();
  }

  CurrentComponent=0;

  // set the wavelength
  SetDiffractionWaveLength(Wavelength);

  // get the indices to the scattering factors for all (pseudo-)atoms
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    PseudoAtoms[i].ScatteringType=GetScatteringNumber(PseudoAtoms[i].ChemicalElement);
    PseudoAtoms[i].AnomalousScatteringType=GetAnomalousScatteringNumber(PseudoAtoms[i].PrintToPDBName);
  }

  InitializeEwald(EwaldPrecision,EwaldAutomatic);
	AllocateEwaldMemory();
	
	if (UseDynamicGrid) 
		MakeDynamicGrid();
	
	if(Framework[0].FrameworkModel==GRID)
	{
		CurrentSystem=0;
		ReadVDWGrid();
		if(UseDelaunayGrid)
			ReadDelaunayGrid();  
		
		if(ChargeMethod!=NONE)
			ReadCoulombGrid();
		
		if(BlockEnergyGrids)
			BlockingVDWGrid();
		
		if (EBCBMC) {
			CreateEBCBMCProbBins();
		}
		
	}

  if(Restart)
  {
    if(RestartStyle==RASPA_STYLE)
      ReadRestartFile();
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
      InitializeReplicaBox();
  }

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {

    if(PutMoleculeOnBarrier[CurrentSystem])
    {
      s=BarrierPosition[CurrentSystem];

      do
      {
        NumberOfBeadsAlreadyPlaced=0;
        CurrentAdsorbateMolecule=0;
        CurrentComponent=0;
        StartingBead=Components[CurrentComponent].StartingBead;
        NewPosition[CurrentSystem][StartingBead].x=UnitCellBox[CurrentSystem].ax*s.x+UnitCellBox[CurrentSystem].bx*s.y+UnitCellBox[CurrentSystem].cx*s.z;
        NewPosition[CurrentSystem][StartingBead].y=UnitCellBox[CurrentSystem].ay*s.x+UnitCellBox[CurrentSystem].by*s.y+UnitCellBox[CurrentSystem].cy*s.z;
        NewPosition[CurrentSystem][StartingBead].z=UnitCellBox[CurrentSystem].az*s.x+UnitCellBox[CurrentSystem].bz*s.y+UnitCellBox[CurrentSystem].cz*s.z;
        GrowMolecule(2);
      }
      while(OVERLAP==TRUE);
      InsertAdsorbateMolecule();
      Components[CurrentComponent].CreateNumberOfMolecules[CurrentSystem]--;
    }
    SetBarierNormal();
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    for(j=0;j<NumberOfComponents;j++)
    {
      if(Components[j].CreateNumberOfMolecules[i]>0)
      {
        CurrentSystem=i;
        CalculateAnisotropicSites();
        if(Components[j].ExtraFrameworkMolecule)
          MakeInitialCations(Components[j].CreateNumberOfMolecules[i],j);
        else
          MakeInitialAdsorbates(Components[j].CreateNumberOfMolecules[i],j);
      }
    }
  }


  // compute the Inertia-tensors and quaternions for all rigid molecules
  ComputeQuaternions();

  CheckConfigMoves();

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    PrintRestartFile();

    // compute the intial center-of-mass position
    Framework[CurrentSystem].IntialCenterOfMassPosition=GetFrameworkCenterOfMass();

    IntialCenterOfMassPosition[CurrentSystem]=GetCenterOfMassCurrentSystem();

    if(Ensemble[CurrentSystem]==NPTPR) BoundaryCondition[CurrentSystem]=TRICLINIC;
    if(ProbabilityBoxShapeChangeMove>0.0) BoundaryCondition[CurrentSystem]=TRICLINIC;
  }

  // if the boundary-condition is still not set, let it default to 'triclinic'
  for(i=0;i<NumberOfSystems;i++)
  {
    if(BoundaryCondition[i]==UNINITIALIZED_BOUNDARY_CONDITION)
      BoundaryCondition[i]=TRICLINIC;
  }

  for(i=0;i<NumberOfSystems;i++)
    WriteVTK(i);

  for(i=0;i<NumberOfSystems;i++)
  {
    if(Framework[i].FrameworkModel==FLEXIBLE)
    {
      // fill in the atom that are kept fixed in e.g. minimization
      // set overall fixed/free setting for framework
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        if(Framework[i].FrameworkModels[j]==FLEXIBLE)
        {
          for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
            Framework[i].Atoms[j][k].Fixed=FALSE;

          switch(FrameworkFixedInitialization[i][j])
          {
            case FIXED:
              for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
                Framework[i].Atoms[j][k].Fixed=TRUE;
              break;
            case FREE:
              for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
                Framework[i].Atoms[j][k].Fixed=FALSE;
              break;
            default:
              break;
          }
        }
        else
        {
          for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
            Framework[i].Atoms[j][k].Fixed=TRUE;

          switch(FrameworkFixedInitialization[i][j])
          {
            case FIXED:
              for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
                Framework[i].Atoms[j][k].Fixed=TRUE;
              break;
            case FREE:
              for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
                Framework[i].Atoms[j][k].Fixed=FALSE;
              break;
            default:
              break;
          }
        }
      }
    }
    else 
    {
      // set overall fixed/free setting for framework
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
          Framework[i].Atoms[j][k].Fixed=TRUE;
      }
    }
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    // set overall fixed/free setting for adsorbates
    for(m=0;m<NumberOfAdsorbateMolecules[i];m++)
    {
      Type=Adsorbates[i][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          if(AdsorbateFixedInitialization[i]==FIXED)
          {
            Adsorbates[i][m].Groups[l].FixedCenterOfMass=TRUE;
            Adsorbates[i][m].Groups[l].FixedOrientation=TRUE;
          }
          else
          {
            Adsorbates[i][m].Groups[l].FixedCenterOfMass=FALSE;
            Adsorbates[i][m].Groups[l].FixedOrientation=FALSE;
          }
        }
        else
        {
          for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
          {
            A=Components[Type].Groups[l].Atoms[k];
            if(AdsorbateFixedInitialization[i]==FIXED)
              Adsorbates[i][m].Atoms[k].Fixed=TRUE;
            else
              Adsorbates[i][m].Atoms[k].Fixed=FALSE;
          }
        }
      }
    }

    // set overall fixed/free setting for cations
    for(m=0;m<NumberOfCationMolecules[i];m++)
    {
      Type=Cations[i][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          if(CationFixedInitialization[i]==FIXED)
          {
            Cations[i][m].Groups[l].FixedCenterOfMass=TRUE;
            Cations[i][m].Groups[l].FixedOrientation=TRUE;
          }
          else
          {
            Cations[i][m].Groups[l].FixedCenterOfMass=FALSE;
            Cations[i][m].Groups[l].FixedOrientation=FALSE;
          }
        }
        else
        {
          for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
          {
            A=Components[Type].Groups[l].Atoms[k];
            if(CationFixedInitialization[i]==FIXED)
              Cations[i][m].Atoms[k].Fixed=TRUE;
            else
              Cations[i][m].Atoms[k].Fixed=FALSE;
          }
        }
      }
    }

    // overwrite based on atom type

    for(f1=0;f1<Framework[i].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[i].NumberOfAtoms[f1];j++)
      {
        if(IsActiveAtomType(Framework[i].Atoms[f1][j].Type))
          Framework[CurrentSystem].Atoms[f1][j].Fixed=FALSE;
        if(IsFixedAtomType(Framework[i].Atoms[f1][j].Type))
          Framework[CurrentSystem].Atoms[f1][j].Fixed=TRUE;
      }
    }

    for(m=0;m<NumberOfAdsorbateMolecules[i];m++)
    {
      Type=Adsorbates[i][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A1=Components[Type].Groups[l].Atoms[k];

          // handle fixed atoms
          if(IsActiveAtomType(Adsorbates[i][m].Atoms[A1].Type))
          {
            if(Components[Type].Groups[l].Rigid) // rigid unit, set all atoms as 'active'
            {
              for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
              {
                 B1=Components[Type].Groups[l].Atoms[j];
                 Adsorbates[i][m].Atoms[B1].Fixed=FALSE;
              }
            }
            else Adsorbates[i][m].Atoms[A1].Fixed=FALSE;
          }
          // handle fixed atoms
          if(IsFixedAtomType(Adsorbates[i][m].Atoms[A1].Type))
          {
            if(Components[Type].Groups[l].Rigid) // rigid unit, set all atoms as 'fixed'
            {
              for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
              {
                 B1=Components[Type].Groups[l].Atoms[j];
                 Adsorbates[i][m].Atoms[B1].Fixed=TRUE;
              }
            }
            else Adsorbates[i][m].Atoms[A1].Fixed=TRUE;
          }
        }
      }
    }

    for(m=0;m<NumberOfCationMolecules[i];m++)
    {
      Type=Cations[i][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A1=Components[Type].Groups[l].Atoms[k];

          // handle fixed atoms
          if(IsActiveAtomType(Cations[i][m].Atoms[A1].Type))
          {
            if(Components[Type].Groups[l].Rigid) // rigid unit, set all atoms as 'active'
            {
              for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
              {
                 B1=Components[Type].Groups[l].Atoms[j];
                 Cations[i][m].Atoms[B1].Fixed=FALSE;
              }
            }
            else Cations[i][m].Atoms[A1].Fixed=FALSE;
          }
          // handle fixed atoms
          if(IsFixedAtomType(Cations[i][m].Atoms[A1].Type))
          {
            if(Components[Type].Groups[l].Rigid) // rigid unit, set all atoms as 'fixed'
            {
              for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
              {
                 B1=Components[Type].Groups[l].Atoms[j];
                 Cations[i][m].Atoms[B1].Fixed=TRUE;
              }
            }
            else Cations[i][m].Atoms[A1].Fixed=TRUE;
          }
        }
      }
    }

    // overwrite individual settings for the framework
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
    {
      for(k=0;k<NumberOfActiveFrameworkAtoms[i][j];k++)
        Framework[i].Atoms[j][ActiveFrameworkAtoms[i][j][k]].Fixed=FALSE;

      for(k=0;k<NumberOfFixedFrameworkAtoms[i][j];k++)
        Framework[i].Atoms[j][FixedFrameworkAtoms[i][j][k]].Fixed=TRUE;
    }

    // overwrite individual settings for the adsorbates
    for(k=0;k<NumberOfActiveAdsorbateMolecules[i];k++)
    {
      A1=ActiveAdsorbateMolecules[i][k];

      Type=Adsorbates[i][A1].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          Adsorbates[i][A1].Groups[l].FixedCenterOfMass=FALSE;
          Adsorbates[i][A1].Groups[l].FixedOrientation=FALSE;
        }
        else
        {
          for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
            Adsorbates[i][A1].Atoms[Components[Type].Groups[l].Atoms[j]].Fixed=FALSE;
        }
      }
    }
    for(k=0;k<NumberOfFixedAdsorbateMolecules[i];k++)
    {
      A1=FixedAdsorbateMolecules[i][k];

      Type=Adsorbates[i][A1].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          Adsorbates[i][A1].Groups[l].FixedCenterOfMass=TRUE;
          Adsorbates[i][A1].Groups[l].FixedOrientation=TRUE;
        }
        else
        {
          for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
            Adsorbates[i][A1].Atoms[Components[Type].Groups[l].Atoms[j]].Fixed=TRUE;
        }
      }
    }
    for(k=0;k<NumberOfActiveAdsorbateAtoms[i];k++)
    {
      A1=ActiveAdsorbateAtoms[i][k].A;
      B1=ActiveAdsorbateAtoms[i][k].B;
      Adsorbates[i][A1].Atoms[B1].Fixed=FALSE;
    }
    for(k=0;k<NumberOfFixedAdsorbateAtoms[i];k++)
    {
      A1=FixedAdsorbateAtoms[i][k].A;
      B1=FixedAdsorbateAtoms[i][k].B;
      Adsorbates[i][A1].Atoms[B1].Fixed=TRUE;
    }
    for(k=0;k<NumberOfActiveAdsorbateGroups[i];k++)
    {
      A1=ActiveAdsorbateGroups[i][k].A;
      B1=ActiveAdsorbateGroups[i][k].B;
      Adsorbates[i][A1].Groups[B1].FixedCenterOfMass=FALSE;
      Adsorbates[i][A1].Groups[B1].FixedOrientation=FALSE;
    }
    for(k=0;k<NumberOfFixedAdsorbateGroups[i];k++)
    {
      A1=FixedAdsorbateGroups[i][k].A;
      B1=FixedAdsorbateGroups[i][k].B;
      Adsorbates[i][A1].Groups[B1].FixedCenterOfMass=TRUE;
      Adsorbates[i][A1].Groups[B1].FixedOrientation=TRUE;
    }
    for(k=0;k<NumberOfActiveAdsorbateGroupsCenterOfMass[i];k++)
    {
      A1=ActiveAdsorbateGroupsCenterOfMass[i][k].A;
      B1=ActiveAdsorbateGroupsCenterOfMass[i][k].B;
      Adsorbates[i][A1].Groups[B1].FixedCenterOfMass=FALSE;
    }
    for(k=0;k<NumberOfFixedAdsorbateGroupsCenterOfMass[i];k++)
    {
      A1=FixedAdsorbateGroupsCenterOfMass[i][k].A;
      B1=FixedAdsorbateGroupsCenterOfMass[i][k].B;
      Adsorbates[i][A1].Groups[B1].FixedCenterOfMass=TRUE;
    }
    for(k=0;k<NumberOfActiveAdsorbateGroupsOrientation[i];k++)
    {
      A1=ActiveAdsorbateGroupsOrientation[i][k].A;
      B1=ActiveAdsorbateGroupsOrientation[i][k].B;
      Adsorbates[i][A1].Groups[B1].FixedOrientation=FALSE;
    }
    for(k=0;k<NumberOfFixedAdsorbateGroupsOrientation[i];k++)
    {
      A1=FixedAdsorbateGroupsOrientation[i][k].A;
      B1=FixedAdsorbateGroupsOrientation[i][k].B;
      Adsorbates[i][A1].Groups[B1].FixedOrientation=TRUE;
    }

    // overwrite individual settings for the cations
    for(k=0;k<NumberOfActiveCationMolecules[i];k++)
    {
      A1=ActiveCationMolecules[i][k];

      Type=Cations[i][A1].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          Cations[i][A1].Groups[l].FixedCenterOfMass=FALSE;
          Cations[i][A1].Groups[l].FixedOrientation=FALSE;
        }
        else
        {
          for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
            Cations[i][A1].Atoms[Components[Type].Groups[l].Atoms[j]].Fixed=FALSE;
        }
      }
    }
    for(k=0;k<NumberOfFixedCationMolecules[i];k++)
    {
      A1=FixedCationMolecules[i][k];

      Type=Cations[i][A1].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          Cations[i][A1].Groups[l].FixedCenterOfMass=TRUE;
          Cations[i][A1].Groups[l].FixedOrientation=TRUE;
        }
        else
        {
          for(j=0;j<Components[Type].Groups[l].NumberOfGroupAtoms;j++)
            Cations[i][A1].Atoms[Components[Type].Groups[l].Atoms[j]].Fixed=TRUE;
        }
      }
    }
    for(k=0;k<NumberOfActiveCationAtoms[i];k++)
    {
      A1=ActiveCationAtoms[i][k].A;
      B1=ActiveCationAtoms[i][k].B;
      Cations[i][A1].Atoms[B1].Fixed=FALSE;
    }
    for(k=0;k<NumberOfFixedCationAtoms[i];k++)
    {
      A1=FixedCationAtoms[i][k].A;
      B1=FixedCationAtoms[i][k].B;
      Cations[i][A1].Atoms[B1].Fixed=TRUE;
    }
    for(k=0;k<NumberOfActiveCationGroups[i];k++)
    {
      A1=ActiveCationGroups[i][k].A;
      B1=ActiveCationGroups[i][k].B;
      Cations[i][A1].Groups[B1].FixedCenterOfMass=FALSE;
      Cations[i][A1].Groups[B1].FixedOrientation=FALSE;
    }
    for(k=0;k<NumberOfFixedCationGroups[i];k++)
    {
      A1=FixedCationGroups[i][k].A;
      B1=FixedCationGroups[i][k].B;
      Cations[i][A1].Groups[B1].FixedCenterOfMass=TRUE;
      Cations[i][A1].Groups[B1].FixedOrientation=TRUE;
    }
    for(k=0;k<NumberOfActiveCationGroupsCenterOfMass[i];k++)
    {
      A1=ActiveCationGroupsCenterOfMass[i][k].A;
      B1=ActiveCationGroupsCenterOfMass[i][k].B;
      Cations[i][A1].Groups[B1].FixedCenterOfMass=FALSE;
    }
    for(k=0;k<NumberOfFixedCationGroupsCenterOfMass[i];k++)
    {
      A1=FixedCationGroupsCenterOfMass[i][k].A;
      B1=FixedCationGroupsCenterOfMass[i][k].B;
      Cations[i][A1].Groups[B1].FixedCenterOfMass=TRUE;
    }
    for(k=0;k<NumberOfActiveCationGroupsOrientation[i];k++)
    {
      A1=ActiveCationGroupsOrientation[i][k].A;
      B1=ActiveCationGroupsOrientation[i][k].B;
      Cations[i][A1].Groups[B1].FixedOrientation=FALSE;
    }
    for(k=0;k<NumberOfFixedCationGroupsOrientation[i];k++)
    {
      A1=FixedCationGroupsOrientation[i][k].A;
      B1=FixedCationGroupsOrientation[i][k].B;
      Cations[i][A1].Groups[B1].FixedOrientation=TRUE;
    }
  }

  // distance pairs for histograms
  for(i=0;i<NumberOfSystems;i++)
  {
    if(NumberOfDistanceHistogramDefinitions[i]>0)
    {
      DistanceHistogramPairs[i]=(ATOM*(*)[2])calloc(NumberOfDistanceHistogramDefinitions[i],sizeof(ATOM*[2]));
      for(j=0;j<NumberOfDistanceHistogramDefinitions[i];j++)
      {
        intinput1=DistanceHistogramDefinitions[i][j][0][1];
        intinput2=DistanceHistogramDefinitions[i][j][0][2];
        switch(DistanceHistogramDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            DistanceHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            DistanceHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            DistanceHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=DistanceHistogramDefinitions[i][j][1][1];
        intinput4=DistanceHistogramDefinitions[i][j][1][2];
        switch(DistanceHistogramDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            DistanceHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            DistanceHistogramPairs[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            DistanceHistogramPairs[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
      }
    }
  }

  // bend-angle trimers for histograms
  for(i=0;i<NumberOfSystems;i++)
  {
    if(NumberOfBendAngleHistogramDefinitions[i]>0)
    {
      BendAngleHistogramPairs[i]=(ATOM*(*)[3])calloc(NumberOfBendAngleHistogramDefinitions[i],sizeof(ATOM*[3]));
      for(j=0;j<NumberOfBendAngleHistogramDefinitions[i];j++)
      {
        intinput1=BendAngleHistogramDefinitions[i][j][0][1];
        intinput2=BendAngleHistogramDefinitions[i][j][0][2];
        switch(BendAngleHistogramDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=BendAngleHistogramDefinitions[i][j][1][1];
        intinput2=BendAngleHistogramDefinitions[i][j][1][2];
        switch(BendAngleHistogramDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][1]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][1]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=BendAngleHistogramDefinitions[i][j][2][1];
        intinput2=BendAngleHistogramDefinitions[i][j][2][2];
        switch(BendAngleHistogramDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][2]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][2]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            BendAngleHistogramPairs[i][j][2]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
      }
    }
  }

  // Dihedral-angle quads for histograms
  for(i=0;i<NumberOfSystems;i++)
  {
    if(NumberOfDihedralAngleHistogramDefinitions[i]>0)
    {
      DihedralAngleHistogramPairs[i]=(ATOM*(*)[4])calloc(NumberOfDihedralAngleHistogramDefinitions[i],sizeof(ATOM*[4]));
      for(j=0;j<NumberOfDihedralAngleHistogramDefinitions[i];j++)
      {
        intinput1=DihedralAngleHistogramDefinitions[i][j][0][1];
        intinput2=DihedralAngleHistogramDefinitions[i][j][0][2];
        switch(DihedralAngleHistogramDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=DihedralAngleHistogramDefinitions[i][j][1][1];
        intinput2=DihedralAngleHistogramDefinitions[i][j][1][2];
        switch(DihedralAngleHistogramDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][1]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][1]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=DihedralAngleHistogramDefinitions[i][j][2][1];
        intinput2=DihedralAngleHistogramDefinitions[i][j][2][2];
        switch(DihedralAngleHistogramDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][2]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][2]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][2]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=DihedralAngleHistogramDefinitions[i][j][3][1];
        intinput2=DihedralAngleHistogramDefinitions[i][j][3][2];
        switch(DihedralAngleHistogramDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][3]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][3]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            DihedralAngleHistogramPairs[i][j][3]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
      }
    }
  }

  // Dihedral-angle quads for histograms
  for(i=0;i<NumberOfSystems;i++)
  {
    if(NumberOfAngleBetweenPlanesHistogramDefinitions[i]>0)
    {
      AngleBetweenPlanesHistogramPairs[i]=(ATOM*(*)[6])calloc(NumberOfAngleBetweenPlanesHistogramDefinitions[i],sizeof(ATOM*[6]));
      for(j=0;j<NumberOfAngleBetweenPlanesHistogramDefinitions[i];j++)
      {
        intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][0][1];
        intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][0][2];
        switch(AngleBetweenPlanesHistogramDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][1][1];
        intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][1][2];
        switch(AngleBetweenPlanesHistogramDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][1]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][1]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][1]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][2][1];
        intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][2][2];
        switch(AngleBetweenPlanesHistogramDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][2]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][2]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][2]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][3][1];
        intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][3][2];
        switch(AngleBetweenPlanesHistogramDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][3]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][3]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][3]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][4][1];
        intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][4][2];
        switch(AngleBetweenPlanesHistogramDefinitions[i][j][4][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][4]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][4]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][4]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }

        intinput1=AngleBetweenPlanesHistogramDefinitions[i][j][5][1];
        intinput2=AngleBetweenPlanesHistogramDefinitions[i][j][5][2];
        switch(AngleBetweenPlanesHistogramDefinitions[i][j][5][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][5]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][5]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            AngleBetweenPlanesHistogramPairs[i][j][5]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
      }
    }
  }



  for(i=0;i<NumberOfSystems;i++)
  {
    if(NumberOfDistanceConstraints[i]>0)
    {
      DistanceConstraints[i]=(ATOM*(*)[2])calloc(NumberOfDistanceConstraints[i],sizeof(ATOM*[2]));
      DistanceConstraintsDerivatives[i]=(VECTOR(*)[2])calloc(NumberOfDistanceConstraints[i],sizeof(VECTOR[2]));
      for(j=0;j<NumberOfDistanceConstraints[i];j++)
      {
        intinput1=DistanceDefinitions[i][j][0][1];
        intinput2=DistanceDefinitions[i][j][0][2];
        switch(DistanceDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            DistanceConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            DistanceConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            DistanceConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=DistanceDefinitions[i][j][1][1];
        intinput4=DistanceDefinitions[i][j][1][2];
        switch(DistanceDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            DistanceConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            DistanceConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            DistanceConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
      }
    }

    if(NumberOfAngleConstraints[i]>0)
    {
      AngleConstraints[i]=(ATOM*(*)[3])calloc(NumberOfAngleConstraints[i],sizeof(ATOM*[3]));
      AngleConstraintsDerivatives[i]=(VECTOR(*)[3])calloc(NumberOfAngleConstraints[i],sizeof(VECTOR[3]));
      for(j=0;j<NumberOfAngleConstraints[i];j++)
      {
        intinput1=AngleDefinitions[i][j][0][1];
        intinput2=AngleDefinitions[i][j][0][2];
        switch(AngleDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            AngleConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            AngleConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            AngleConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=AngleDefinitions[i][j][1][1];
        intinput4=AngleDefinitions[i][j][1][2];
        switch(AngleDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            AngleConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            AngleConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            AngleConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=AngleDefinitions[i][j][2][1];
        intinput6=AngleDefinitions[i][j][2][2];
        switch(AngleDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            AngleConstraints[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            AngleConstraints[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            AngleConstraints[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
      }
    }

    if(NumberOfDihedralConstraints[i]>0)
    {
      DihedralConstraints[i]=(ATOM*(*)[4])calloc(NumberOfDihedralConstraints[i],sizeof(ATOM*[4]));
      DihedralConstraintsDerivatives[i]=(VECTOR(*)[4])calloc(NumberOfDihedralConstraints[i],sizeof(VECTOR[4]));
      for(j=0;j<NumberOfDihedralConstraints[i];j++)
      {
        intinput1=DihedralDefinitions[i][j][0][1];
        intinput2=DihedralDefinitions[i][j][0][2];
        switch(DihedralDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            DihedralConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            DihedralConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            DihedralConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=DihedralDefinitions[i][j][1][1];
        intinput4=DihedralDefinitions[i][j][1][2];
        switch(DihedralDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            DihedralConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            DihedralConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            DihedralConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=DihedralDefinitions[i][j][2][1];
        intinput6=DihedralDefinitions[i][j][2][2];
        switch(DihedralDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            DihedralConstraints[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            DihedralConstraints[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            DihedralConstraints[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
        intinput7=DihedralDefinitions[i][j][3][1];
        intinput8=DihedralDefinitions[i][j][3][2];
        switch(DihedralDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput7,intinput8);
            DihedralConstraints[i][j][3]=&Framework[i].Atoms[intinput7][intinput8];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput7,intinput8);
            DihedralConstraints[i][j][3]=&Adsorbates[i][intinput7].Atoms[intinput8];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput7,intinput8);
            DihedralConstraints[i][j][3]=&Cations[i][intinput7].Atoms[intinput8];
            break;
        }
      }
    }

    if(NumberOfImproperDihedralConstraints[i]>0)
    {
      ImproperDihedralConstraints[i]=(ATOM*(*)[4])calloc(NumberOfImproperDihedralConstraints[i],sizeof(ATOM*[4]));
      ImproperDihedralConstraintsDerivatives[i]=(VECTOR(*)[4])calloc(NumberOfImproperDihedralConstraints[i],sizeof(VECTOR[4]));
      for(j=0;j<NumberOfImproperDihedralConstraints[i];j++)
      {
        intinput1=ImproperDihedralDefinitions[i][j][0][1];
        intinput2=ImproperDihedralDefinitions[i][j][0][2];
        switch(ImproperDihedralDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            ImproperDihedralConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            ImproperDihedralConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            ImproperDihedralConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=ImproperDihedralDefinitions[i][j][1][1];
        intinput4=ImproperDihedralDefinitions[i][j][1][2];
        switch(ImproperDihedralDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            ImproperDihedralConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            ImproperDihedralConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            ImproperDihedralConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=ImproperDihedralDefinitions[i][j][2][1];
        intinput6=ImproperDihedralDefinitions[i][j][2][2];
        switch(ImproperDihedralDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            ImproperDihedralConstraints[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            ImproperDihedralConstraints[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            ImproperDihedralConstraints[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
        intinput7=ImproperDihedralDefinitions[i][j][3][1];
        intinput8=ImproperDihedralDefinitions[i][j][3][2];
        switch(ImproperDihedralDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput7,intinput8);
            ImproperDihedralConstraints[i][j][3]=&Framework[i].Atoms[intinput7][intinput8];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput7,intinput8);
            ImproperDihedralConstraints[i][j][3]=&Adsorbates[i][intinput7].Atoms[intinput8];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput7,intinput8);
            ImproperDihedralConstraints[i][j][3]=&Cations[i][intinput7].Atoms[intinput8];
            break;
        }
      }
    }

    if(NumberOfInversionBendConstraints[i]>0)
    {
      InversionBendConstraints[i]=(ATOM*(*)[4])calloc(NumberOfInversionBendConstraints[i],sizeof(ATOM*[4]));
      InversionBendConstraintsDerivatives[i]=(VECTOR(*)[4])calloc(NumberOfInversionBendConstraints[i],sizeof(VECTOR[4]));
      for(j=0;j<NumberOfInversionBendConstraints[i];j++)
      {
        intinput1=InversionBendDefinitions[i][j][0][1];
        intinput2=InversionBendDefinitions[i][j][0][2];
        switch(InversionBendDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            InversionBendConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            InversionBendConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            InversionBendConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=InversionBendDefinitions[i][j][1][1];
        intinput4=InversionBendDefinitions[i][j][1][2];
        switch(InversionBendDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            InversionBendConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            InversionBendConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            InversionBendConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=InversionBendDefinitions[i][j][2][1];
        intinput6=InversionBendDefinitions[i][j][2][2];
        switch(InversionBendDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            InversionBendConstraints[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            InversionBendConstraints[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            InversionBendConstraints[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
        intinput7=InversionBendDefinitions[i][j][3][1];
        intinput8=InversionBendDefinitions[i][j][3][2];
        switch(InversionBendDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput7,intinput8);
            InversionBendConstraints[i][j][3]=&Framework[i].Atoms[intinput7][intinput8];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput7,intinput8);
            InversionBendConstraints[i][j][3]=&Adsorbates[i][intinput7].Atoms[intinput8];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput7,intinput8);
            InversionBendConstraints[i][j][3]=&Cations[i][intinput7].Atoms[intinput8];
            break;
        }
      }
    }


    if(NumberOfOutOfPlaneDistanceConstraints[i]>0)
    {
      OutOfPlaneDistanceConstraints[i]=(ATOM*(*)[4])calloc(NumberOfOutOfPlaneDistanceConstraints[i],sizeof(ATOM*[4]));
      OutOfPlaneDistanceConstraintsDerivatives[i]=(VECTOR(*)[4])calloc(NumberOfOutOfPlaneDistanceConstraints[i],sizeof(VECTOR[4]));
      for(j=0;j<NumberOfOutOfPlaneDistanceConstraints[i];j++)
      {
        intinput1=OutOfPlaneDistanceDefinitions[i][j][0][1];
        intinput2=OutOfPlaneDistanceDefinitions[i][j][0][2];
        switch(OutOfPlaneDistanceDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            OutOfPlaneDistanceConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            OutOfPlaneDistanceConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            OutOfPlaneDistanceConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=OutOfPlaneDistanceDefinitions[i][j][1][1];
        intinput4=OutOfPlaneDistanceDefinitions[i][j][1][2];
        switch(OutOfPlaneDistanceDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            OutOfPlaneDistanceConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            OutOfPlaneDistanceConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            OutOfPlaneDistanceConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=OutOfPlaneDistanceDefinitions[i][j][2][1];
        intinput6=OutOfPlaneDistanceDefinitions[i][j][2][2];
        switch(OutOfPlaneDistanceDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            OutOfPlaneDistanceConstraints[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            OutOfPlaneDistanceConstraints[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            OutOfPlaneDistanceConstraints[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
        intinput7=OutOfPlaneDistanceDefinitions[i][j][3][1];
        intinput8=OutOfPlaneDistanceDefinitions[i][j][3][2];
        switch(OutOfPlaneDistanceDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput7,intinput8);
            OutOfPlaneDistanceConstraints[i][j][3]=&Framework[i].Atoms[intinput7][intinput8];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput7,intinput8);
            OutOfPlaneDistanceConstraints[i][j][3]=&Adsorbates[i][intinput7].Atoms[intinput8];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput7,intinput8);
            OutOfPlaneDistanceConstraints[i][j][3]=&Cations[i][intinput7].Atoms[intinput8];
            break;
        }
      }
    }


    if(NumberOfHarmonicDistanceConstraints[i]>0)
    {
      HarmonicDistanceConstraints[i]=(ATOM*(*)[2])calloc(NumberOfHarmonicDistanceConstraints[i],sizeof(ATOM*[2]));

      for(j=0;j<NumberOfHarmonicDistanceConstraints[i];j++)
      {
        intinput1=HarmonicDistanceDefinitions[i][j][0][1];
        intinput2=HarmonicDistanceDefinitions[i][j][0][2];
        switch(HarmonicDistanceDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            HarmonicDistanceConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            HarmonicDistanceConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            HarmonicDistanceConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=HarmonicDistanceDefinitions[i][j][1][1];
        intinput4=HarmonicDistanceDefinitions[i][j][1][2];
        switch(HarmonicDistanceDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            HarmonicDistanceConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            HarmonicDistanceConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            HarmonicDistanceConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
      }
    }

    if(NumberOfHarmonicAngleConstraints[i]>0)
    {
      HarmonicAngleConstraints[i]=(ATOM*(*)[3])calloc(NumberOfHarmonicAngleConstraints[i],sizeof(ATOM*[3]));
      for(j=0;j<NumberOfHarmonicAngleConstraints[i];j++)
      {
        intinput1=HarmonicAngleDefinitions[i][j][0][1];
        intinput2=HarmonicAngleDefinitions[i][j][0][2];
        switch(HarmonicAngleDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            HarmonicAngleConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            HarmonicAngleConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            HarmonicAngleConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=HarmonicAngleDefinitions[i][j][1][1];
        intinput4=HarmonicAngleDefinitions[i][j][1][2];
        switch(HarmonicAngleDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            HarmonicAngleConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            HarmonicAngleConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            HarmonicAngleConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=HarmonicAngleDefinitions[i][j][2][1];
        intinput6=HarmonicAngleDefinitions[i][j][2][2];
        switch(HarmonicAngleDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            HarmonicAngleConstraints[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            HarmonicAngleConstraints[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            HarmonicAngleConstraints[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
      }
    }

    if(NumberOfHarmonicDihedralConstraints[i]>0)
    {
      HarmonicDihedralConstraints[i]=(ATOM*(*)[4])calloc(NumberOfHarmonicDihedralConstraints[i],sizeof(ATOM*[4]));
      for(j=0;j<NumberOfHarmonicDihedralConstraints[i];j++)
      {
        intinput1=HarmonicDihedralDefinitions[i][j][0][1];
        intinput2=HarmonicDihedralDefinitions[i][j][0][2];
        switch(HarmonicDihedralDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            HarmonicDihedralConstraints[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            HarmonicDihedralConstraints[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            HarmonicDihedralConstraints[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=HarmonicDihedralDefinitions[i][j][1][1];
        intinput4=HarmonicDihedralDefinitions[i][j][1][2];
        switch(HarmonicDihedralDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            HarmonicDihedralConstraints[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            HarmonicDihedralConstraints[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            HarmonicDihedralConstraints[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=HarmonicDihedralDefinitions[i][j][2][1];
        intinput6=HarmonicDihedralDefinitions[i][j][2][2];
        switch(HarmonicDihedralDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            HarmonicDihedralConstraints[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            HarmonicDihedralConstraints[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            HarmonicDihedralConstraints[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
        intinput7=HarmonicDihedralDefinitions[i][j][3][1];
        intinput8=HarmonicDihedralDefinitions[i][j][3][2];
        switch(HarmonicDihedralDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput7,intinput8);
            HarmonicDihedralConstraints[i][j][3]=&Framework[i].Atoms[intinput7][intinput8];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput7,intinput8);
            HarmonicDihedralConstraints[i][j][3]=&Adsorbates[i][intinput7].Atoms[intinput8];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput7,intinput8);
            HarmonicDihedralConstraints[i][j][3]=&Cations[i][intinput7].Atoms[intinput8];
            break;
        }
      }
    }


    if(NumberOfTwoPointDihedralDefinitions[i]>0)
    {
      TwoPointDihedrals[i]=(ATOM*(*)[6])calloc(NumberOfTwoPointDihedralDefinitions[i],sizeof(ATOM*[6]));
      for(j=0;j<NumberOfTwoPointDihedralDefinitions[i];j++)
      {
        intinput1=TwoPointDihedralDefinitions[i][j][0][1];
        intinput2=TwoPointDihedralDefinitions[i][j][0][2];
        switch(TwoPointDihedralDefinitions[i][j][0][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput1,intinput2);
            TwoPointDihedrals[i][j][0]=&Framework[i].Atoms[intinput1][intinput2];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
            TwoPointDihedrals[i][j][0]=&Adsorbates[i][intinput1].Atoms[intinput2];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput1,intinput2);
            TwoPointDihedrals[i][j][0]=&Cations[i][intinput1].Atoms[intinput2];
            break;
        }
        intinput3=TwoPointDihedralDefinitions[i][j][1][1];
        intinput4=TwoPointDihedralDefinitions[i][j][1][2];
        switch(TwoPointDihedralDefinitions[i][j][1][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput3,intinput4);
            TwoPointDihedrals[i][j][1]=&Framework[i].Atoms[intinput3][intinput4];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput3,intinput4);
            TwoPointDihedrals[i][j][1]=&Adsorbates[i][intinput3].Atoms[intinput4];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput3,intinput4);
            TwoPointDihedrals[i][j][1]=&Cations[i][intinput3].Atoms[intinput4];
            break;
        }
        intinput5=TwoPointDihedralDefinitions[i][j][2][1];
        intinput6=TwoPointDihedralDefinitions[i][j][2][2];
        switch(TwoPointDihedralDefinitions[i][j][2][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput5,intinput6);
            TwoPointDihedrals[i][j][2]=&Framework[i].Atoms[intinput5][intinput6];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput5,intinput6);
            TwoPointDihedrals[i][j][2]=&Adsorbates[i][intinput5].Atoms[intinput6];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput5,intinput6);
            TwoPointDihedrals[i][j][2]=&Cations[i][intinput5].Atoms[intinput6];
            break;
        }
        intinput7=TwoPointDihedralDefinitions[i][j][3][1];
        intinput8=TwoPointDihedralDefinitions[i][j][3][2];
        switch(TwoPointDihedralDefinitions[i][j][3][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput7,intinput8);
            TwoPointDihedrals[i][j][3]=&Framework[i].Atoms[intinput7][intinput8];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput7,intinput8);
            TwoPointDihedrals[i][j][3]=&Adsorbates[i][intinput7].Atoms[intinput8];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput7,intinput8);
            TwoPointDihedrals[i][j][3]=&Cations[i][intinput7].Atoms[intinput8];
            break;
        }
        intinput9=TwoPointDihedralDefinitions[i][j][4][1];
        intinput10=TwoPointDihedralDefinitions[i][j][4][2];
        switch(TwoPointDihedralDefinitions[i][j][4][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput9,intinput10);
            TwoPointDihedrals[i][j][4]=&Framework[i].Atoms[intinput9][intinput10];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput9,intinput10);
            TwoPointDihedrals[i][j][4]=&Adsorbates[i][intinput9].Atoms[intinput10];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput9,intinput10);
            TwoPointDihedrals[i][j][4]=&Cations[i][intinput9].Atoms[intinput10];
            break;
        }
        intinput11=TwoPointDihedralDefinitions[i][j][5][1];
        intinput12=TwoPointDihedralDefinitions[i][j][5][2];
        switch(TwoPointDihedralDefinitions[i][j][5][0])
        {
          case FRAMEWORK:
            CheckConstraintInputFramework(Framework,intinput11,intinput12);
            TwoPointDihedrals[i][j][5]=&Framework[i].Atoms[intinput11][intinput12];
            break;
          case ADSORBATE:
            CheckConstraintInputAdsorbate(Adsorbates,intinput11,intinput12);
            TwoPointDihedrals[i][j][5]=&Adsorbates[i][intinput11].Atoms[intinput12];
            break;
          case CATION:
            CheckConstraintInputCation(Cations,intinput11,intinput12);
            TwoPointDihedrals[i][j][5]=&Cations[i][intinput11].Atoms[intinput12];
            break;
        }
      }
    }

    for(j=0;j<NumberOfDistanceConstraints[i];j++)
      if(DistanceConstraintParameter[i][j]<0.0)
        DistanceConstraintParameter[i][j]=ReturnBondDistance(DistanceConstraints[i][j][0]->Position,
              DistanceConstraints[i][j][1]->Position);

    for(j=0;j<NumberOfAngleConstraints[i];j++)
      if(AngleConstraintParameter[i][j]<-180.0)
        AngleConstraintParameter[i][j]=ReturnBendAngle(AngleConstraints[i][j][0]->Position,
              AngleConstraints[i][j][1]->Position,AngleConstraints[i][j][2]->Position);

    for(j=0;j<NumberOfDihedralConstraints[i];j++)
      if(DihedralConstraintParameter[i][j]<-180.0)
        DihedralConstraintParameter[i][j]=ReturnDihedralAngle(DihedralConstraints[i][j][0]->Position,
              DihedralConstraints[i][j][1]->Position,DihedralConstraints[i][j][2]->Position,
              DihedralConstraints[i][j][3]->Position);

    for(j=0;j<NumberOfImproperDihedralConstraints[i];j++)
      if(ImproperDihedralConstraintParameter[i][j]<-180.0)
        ImproperDihedralConstraintParameter[i][j]=ReturnDihedralAngle(ImproperDihedralConstraints[i][j][0]->Position,
              ImproperDihedralConstraints[i][j][1]->Position,ImproperDihedralConstraints[i][j][2]->Position,
              ImproperDihedralConstraints[i][j][3]->Position);

    for(j=0;j<NumberOfInversionBendConstraints[i];j++)
      if(InversionBendConstraintParameter[i][j]<-180.0)
        InversionBendConstraintParameter[i][j]=ReturnInversionBendAngle(InversionBendConstraints[i][j][0]->Position,
              InversionBendConstraints[i][j][1]->Position,InversionBendConstraints[i][j][2]->Position,
              InversionBendConstraints[i][j][3]->Position);



    for(j=0;j<NumberOfHarmonicDistanceConstraints[i];j++)
      if(HarmonicDistanceConstraintParameters[i][j][1]<0.0)
        HarmonicDistanceConstraintParameters[i][j][1]=ReturnBondDistance(
              HarmonicDistanceConstraints[i][j][0]->Position,
              HarmonicDistanceConstraints[i][j][1]->Position);

    for(j=0;j<NumberOfHarmonicAngleConstraints[i];j++)
      if(HarmonicAngleConstraintParameters[i][j][1]<-180.0)
        HarmonicAngleConstraintParameters[i][j][1]=
              ReturnBendAngle(HarmonicAngleConstraints[i][j][0]->Position,
              HarmonicAngleConstraints[i][j][1]->Position,HarmonicAngleConstraints[i][j][2]->Position);

    for(j=0;j<NumberOfHarmonicDihedralConstraints[i];j++)
      if(HarmonicDihedralConstraintParameters[i][j][1]<-180.0)
        HarmonicDihedralConstraintParameters[i][j][1]=
              ReturnDihedralAngle(HarmonicDihedralConstraints[i][j][0]->Position,
              HarmonicDihedralConstraints[i][j][1]->Position,
              HarmonicDihedralConstraints[i][j][2]->Position,
              HarmonicDihedralConstraints[i][j][3]->Position);
  }

  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].RestrictEnantionface)
    {
      intinput1=Components[i].EnantiofaceAtomDefinitions[0][1];
      intinput2=Components[i].EnantiofaceAtomDefinitions[0][2];
      switch(Components[i].EnantiofaceAtomDefinitions[0][0])
      {
        case FRAMEWORK:
          CheckConstraintInputFramework(Framework,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[0]=&Framework[i].Atoms[intinput1][intinput2];
          break;
        case ADSORBATE:
          CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[0]=&Adsorbates[i][intinput1].Atoms[intinput2];
          break;
        case CATION:
          CheckConstraintInputCation(Cations,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[0]=&Cations[i][intinput1].Atoms[intinput2];
          break;
      }

      intinput1=Components[i].EnantiofaceAtomDefinitions[1][1];
      intinput2=Components[i].EnantiofaceAtomDefinitions[1][2];
      switch(Components[i].EnantiofaceAtomDefinitions[1][0])
      {
        case FRAMEWORK:
          CheckConstraintInputFramework(Framework,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[1]=&Framework[i].Atoms[intinput1][intinput2];
          break;
        case ADSORBATE:
          CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[1]=&Adsorbates[i][intinput1].Atoms[intinput2];
          break;
        case CATION:
          CheckConstraintInputCation(Cations,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[1]=&Cations[i][intinput1].Atoms[intinput2];
          break;
      }

      intinput1=Components[i].EnantiofaceAtomDefinitions[2][1];
      intinput2=Components[i].EnantiofaceAtomDefinitions[2][2];
      switch(Components[i].EnantiofaceAtomDefinitions[2][0])
      {
        case FRAMEWORK:
          CheckConstraintInputFramework(Framework,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[2]=&Framework[i].Atoms[intinput1][intinput2];
          break;
        case ADSORBATE:
          CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[2]=&Adsorbates[i][intinput1].Atoms[intinput2];
          break;
        case CATION:
          CheckConstraintInputCation(Cations,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[2]=&Cations[i][intinput1].Atoms[intinput2];
          break;
      }

      intinput1=Components[i].EnantiofaceAtomDefinitions[3][1];
      intinput2=Components[i].EnantiofaceAtomDefinitions[3][2];
      switch(Components[i].EnantiofaceAtomDefinitions[3][0])
      {
        case FRAMEWORK:
          CheckConstraintInputFramework(Framework,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[3]=&Framework[i].Atoms[intinput1][intinput2];
          break;
        case ADSORBATE:
          CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[3]=&Adsorbates[i][intinput1].Atoms[intinput2];
          break;
        case CATION:
          CheckConstraintInputCation(Cations,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[3]=&Cations[i][intinput1].Atoms[intinput2];
          break;
      }

      intinput1=Components[i].EnantiofaceAtomDefinitions[4][1];
      intinput2=Components[i].EnantiofaceAtomDefinitions[4][2];
      switch(Components[i].EnantiofaceAtomDefinitions[4][0])
      {
        case FRAMEWORK:
          CheckConstraintInputFramework(Framework,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[4]=&Framework[i].Atoms[intinput1][intinput2];
          break;
        case ADSORBATE:
          CheckConstraintInputAdsorbate(Adsorbates,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[4]=&Adsorbates[i][intinput1].Atoms[intinput2];
          break;
        case CATION:
          CheckConstraintInputCation(Cations,intinput1,intinput2);
          Components[i].EnantiofaceAtoms[4]=&Cations[i][intinput1].Atoms[intinput2];
          break;
      }
    }
  }

  for(i=0;i<NumberOfSystems;i++)
  {
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
    {
      Framework[i].NumberOfFreeAtoms[j]=0;
      Framework[i].NumberOfFixedAtoms[j]=0;
      for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
      {
        if(Framework[i].Atoms[j][k].Fixed)
          Framework[i].NumberOfFixedAtoms[j]++;
        else
          Framework[i].NumberOfFreeAtoms[j]++;
      }
    }
  }

  // mark also atoms as swapable if the molecule is swapable
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].Swapable)
    {
      for(j=0;j<Components[i].NumberOfAtoms;j++)
        PseudoAtoms[Components[i].Type[j]].Swapable=TRUE;
    }
  }

  // compute the degrees of freedom
  ComputeDegreesOfFreedom();

  // generate warnings
  CreateWarnings();

  // precompute contributions to the Ewald sums from fixed atoms
  // setup k-vectors in advance
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    SetupKVectors();
    EwaldEnergyIon();
    PrecomputeFixedEwaldContributions();
    PrecomputeTotalEwaldContributions();
  }

  if(UseContinuousFraction)
  {
    for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    {
      for(CurrentComponent=0;CurrentComponent<NumberOfComponents;CurrentComponent++)
      {
        Components[CurrentComponent].FractionalMolecule[CurrentSystem]=(int)(RandomNumber()*Components[CurrentComponent].NumberOfMolecules[CurrentSystem]);
        if(Components[CurrentComponent].ExtraFrameworkMolecule)
          Cations[CurrentSystem][Components[CurrentComponent].FractionalMolecule[CurrentSystem]].Fractional=TRUE;
        else
          Adsorbates[CurrentSystem][Components[CurrentComponent].FractionalMolecule[CurrentSystem]].Fractional=TRUE;

        Components[CurrentComponent].Lambda=0.5;
      }
    }
  }

  CurrentSystem=0;
  CurrentComponent=0;
  CurrentFramework=0;
  
  printf("Done reading '%s'\n",inputfilename);
  return 0;
}

void ReadRestartFile(void)
{
  int i;
  int NumberOfComponentsRead;
  int extra_framework_boolean;
  FILE *FilePtrIn;
  char buffer[256];
  char line[1024],keyword[1024],arguments[1024];
  int CurrentComponentRead;
  char ExtraFrameworkMoleculeRead[256];
  int NumberOfMoleculesRead,AtomNumberRead,MoleculeNumberRead,FrameworkAtomNumberRead;
  char ComponentNameRead[256];
  VECTOR tmp;
  double temp1,temp2,temp3;

  extra_framework_boolean=FALSE;
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    sprintf(buffer,"RestartInitial/System_%d",CurrentSystem);
    mkdir(buffer,S_IRWXU);
  }

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    sprintf(buffer,"RestartInitial/System_%d/restart_%s_%d.%d.%d_%lf_%lf",
            CurrentSystem,
            Framework[CurrentSystem].Name[0],
            NumberOfUnitCells[CurrentSystem].x,
            NumberOfUnitCells[CurrentSystem].y,
            NumberOfUnitCells[CurrentSystem].z,
            (double)therm_baro_stats.ExternalTemperature[CurrentSystem],
            (double)therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR);

    if(!(FilePtrIn=fopen(buffer,"r")))
    {
      printf("Could NOT open file: %s\n",buffer);
      exit(0);
    }

    while(fgets(line,1024,FilePtrIn))
    {
      // extract first word
      strcpy(keyword,"keyword");
      sscanf(line,"%s %[^\n]",keyword,arguments);

      if(strcasecmp(keyword,"Box:")==0)
      {
        fgets(line,1024,FilePtrIn);
        fgets(line,1024,FilePtrIn);
        sscanf(line,"%lf %lf %lf\n",&temp1,&temp2,&temp3);
        Box[CurrentSystem].ax=(REAL)temp1;
        Box[CurrentSystem].bx=(REAL)temp2;
        Box[CurrentSystem].cx=(REAL)temp3;
        fgets(line,1024,FilePtrIn);
        sscanf(line,"%lf %lf %lf\n",&temp1,&temp2,&temp3);
        Box[CurrentSystem].ay=(REAL)temp1;
        Box[CurrentSystem].by=(REAL)temp2;
        Box[CurrentSystem].cy=(REAL)temp3;
        fgets(line,1024,FilePtrIn);
        sscanf(line,"%lf %lf %lf\n",&temp1,&temp2,&temp3);
        Box[CurrentSystem].az=(REAL)temp1;
        Box[CurrentSystem].bz=(REAL)temp2;
        Box[CurrentSystem].cz=(REAL)temp3;

        printf("BOX: %g %g %g\n",Box[CurrentSystem].ax,Box[CurrentSystem].by,Box[CurrentSystem].cz);

        InverseBoxMatrix(&Box[CurrentSystem],&InverseBox[CurrentSystem]);

        //UnitCellBox[CurrentSystem]=Box[CurrentSystem];
        //InverseUnitCellBox[CurrentSystem]=InverseBox[CurrentSystem];

        CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

        AlphaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bx);
        BetaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].by);
        GammaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bz);

        // determine boundary conditions from angles
        if((BoundaryCondition[CurrentSystem]!=UNINITIALIZED_BOUNDARY_CONDITION)&&(BoundaryCondition[CurrentSystem]!=FINITE))
        {
          if((fabs(AlphaAngle[CurrentSystem]-90.0)>0.001)||
             (fabs(BetaAngle[CurrentSystem]-90.0)>0.001)||
             (fabs(GammaAngle[CurrentSystem]-90.0)>0.001))
            BoundaryCondition[CurrentSystem]=TRICLINIC;
          else BoundaryCondition[CurrentSystem]=RECTANGULAR;
        }

        if(Framework[CurrentSystem].FrameworkModel==NONE)
          sprintf(Framework[CurrentSystem].Name[0],"Box");
      }

      if(strcasecmp(keyword,"Components:")==0)
      {
        sscanf(arguments,"%d%*[^\n]",&NumberOfComponentsRead);
        if(NumberOfComponentsRead!=NumberOfComponents)
          printf("Warning: NumberOfComponents does not match restart-file !\n");
      }

      if(strcasecmp(keyword,"Volume")==0)
      {
        sscanf(arguments,"change, maximum: %lf\n",&temp1);
        MaximumVolumeChange[CurrentSystem]=(REAL)temp1;
      }
      if(strcasecmp(keyword,"Translation")==0)
      {
        sscanf(arguments,"change component %d [ %s ], maximum: %lf,%lf,%lf\n",
          &CurrentComponent,
          buffer,&temp1,&temp2,&temp3);
        tmp.x=(REAL)temp1;
        tmp.y=(REAL)temp2;
        tmp.z=(REAL)temp3;
        MaximumTranslation[CurrentSystem][CurrentComponent]=tmp;
      }

      // parse "Component: 0     Cation       96 molecules of sodium"
      if(strcasecmp(keyword,"Component:")==0)
      {
        sscanf(arguments,"%d %s %d molecules of %s%*[^\n]",
           &CurrentComponentRead,
           ExtraFrameworkMoleculeRead,
           &NumberOfMoleculesRead,
           ComponentNameRead);
        CurrentComponent=CurrentComponentRead;

        if(!strncasecmp(ExtraFrameworkMoleculeRead,"Cation",MAX2(strlen(ExtraFrameworkMoleculeRead),strlen("Cation"))))
          extra_framework_boolean=TRUE;
        else
          extra_framework_boolean=FALSE;
      }

      if(strcasecmp(keyword,"Framework")==0)
      {
        FrameworkAtomNumberRead=-1;
        sscanf(arguments,"Atom: %d %d %[^\n]",
           &FrameworkAtomNumberRead,
           &CurrentFramework,
           arguments);
        sscanf(arguments,"Position: %lf %lf %lf %[^\n]",&temp1,&temp2,&temp3,arguments);
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Position.x=(REAL)temp1;
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Position.y=(REAL)temp2;
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Position.z=(REAL)temp3;

        sscanf(arguments,"Velocity: %lf %lf %lf %[^\n]",&temp1,&temp2,&temp3,arguments);
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Velocity.x=(REAL)temp1;
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Velocity.y=(REAL)temp2;
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Velocity.z=(REAL)temp3;

        sscanf(arguments,"Force: %lf %lf %lf %[^\n]",&temp1,&temp2,&temp3,arguments);
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Force.x=(REAL)temp1;
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Force.y=(REAL)temp2;
        Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Force.z=(REAL)temp3;

        if(sscanf(arguments,"Charge: %lf%*[^\n]",&temp1))
          Framework[CurrentSystem].Atoms[CurrentFramework][FrameworkAtomNumberRead].Charge=(REAL)temp1;
      }

      if(strcasecmp(keyword,"Molecule:")==0)
      {
        if(NumberOfComponents==0)
        {
          printf("\nReadRestartFile: component present, but not defined in inputfile\n");
          exit(0);
        }

        MoleculeNumberRead=-1;
        sscanf(arguments,"%d %[^\n]",
           &MoleculeNumberRead,
           arguments);

        AtomNumberRead=-1;
        sscanf(arguments,"Atom: %d %[^\n]",
           &AtomNumberRead,
           arguments);

        sscanf(arguments,"Position: %lf %lf %lf %[^\n]",&temp1,&temp2,&temp3,arguments);
        NewPosition[CurrentSystem][AtomNumberRead].x=(REAL)temp1;
        NewPosition[CurrentSystem][AtomNumberRead].y=(REAL)temp2;
        NewPosition[CurrentSystem][AtomNumberRead].z=(REAL)temp3;

        sscanf(arguments,"Velocity: %lf %lf %lf %[^\n]",&temp1,&temp2,&temp3,arguments);
        NewVelocity[CurrentSystem][AtomNumberRead].x=(REAL)temp1;
        NewVelocity[CurrentSystem][AtomNumberRead].y=(REAL)temp2;
        NewVelocity[CurrentSystem][AtomNumberRead].z=(REAL)temp3;

        sscanf(arguments,"Force: %lf %lf %lf %*[^\n]",&temp1,&temp2,&temp3);
        NewForce[CurrentSystem][AtomNumberRead].x=(REAL)temp1;
        NewForce[CurrentSystem][AtomNumberRead].y=(REAL)temp2;
        NewForce[CurrentSystem][AtomNumberRead].z=(REAL)temp3;

        if(AtomNumberRead==Components[CurrentComponent].NumberOfAtoms-1)
        {
          if(extra_framework_boolean)
            InsertCationMolecule();
          else
            InsertAdsorbateMolecule();
          AtomNumberRead=0;
        }
      }
    }
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassAdsorbate(i);
    for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
      UpdateGroupCenterOfMassCation(i);
  }
}

void ReadBinaryRestartFiles(void)
{
  FILE *FilePtr;
  char buffer[1024];
  REAL Check;

  printf("read binary file\n");

  sprintf(buffer,"CrashRestart/binary_restart.dat");
  if((FilePtr=fopen(buffer,"r")))
  {
    ReadRestartConstants(FilePtr);
    ReadRestartSimulation(FilePtr);
    ReadRestartWarnings(FilePtr);
    ReadRestartPseudoAtoms(FilePtr);
    ReadRestartComponent(FilePtr);
    ReadRestartMolecules(FilePtr);
    ReadRestartFramework(FilePtr);
    ReadRestartCBMC(FilePtr);
    ReadRestartEwald(FilePtr);
    ReadRestartStatistics(FilePtr);
    ReadRestartMcMoves(FilePtr);
    ReadRestartSample(FilePtr);
    ReadRestartThermoBarostats(FilePtr);
    ReadRestartEquationOfState(FilePtr);
    ReadRestartGrids(FilePtr);
    ReadRestartMinimization(FilePtr);
    ReadRestartUtils(FilePtr);
    ReadRestartMovies(FilePtr);
    ReadRestartOutput(FilePtr);

    fread(&Check,1,sizeof(REAL),FilePtr);
    printf("Check number: %18.10f\n",Check);
    if(fabs(Check-123456789.0)>1e-10)
    {
      printf("Error in binary restart-file\n");
      exit(0);
    }

    fclose(FilePtr);
  }
  else
  {
    ContinueAfterCrash=FALSE;
    printf("Crash set to false\n");
  }
}

