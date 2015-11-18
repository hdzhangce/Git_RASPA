/*****************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2012 David Dubbeldam, Sofia Calero, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'grids.h' is part of RASPA.

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

#ifndef GRIDS_H
#define GRIDS_H

extern int UseTabularGrid;
extern int UseDynamicGrid;
extern int UseDelaunayGrid;
extern int NumberOfGrids;
extern int *GridTypeList;
extern char (*GridTypeListName)[256];
extern int EBCBMC;
extern int HermiteInterpolation;

extern float *****VDWGrid;
extern int *****DelaunayGrid;
extern int ****DynamicGrid;
extern float ****CoulombGrid;

extern INT_VECTOR3 NumberOfVDWGridPoints;
extern INT_VECTOR3 NumberOfDynamicGridPoints;
extern INT_VECTOR3 NumberOfDelaunayGridPoints;

extern REAL SpacingVDWGrid;
extern REAL SpacingCoulombGrid;
extern REAL SpacingDelaunayGrid;
extern REAL SpacingDynamicGrid;

extern int BlockEnergyGrids;
extern int BlockGridPockets;
extern int BlockGridPores;
extern REAL BlockEnergyGridOverlapCriteria;
extern int NumberOfGridSeeds;
extern VECTOR *GridSeeds;

extern int DelaunayPolyhedraFile;
extern char DelaunayPolyhedraFileName[256];
extern VECTOR *DelaunayUnitCellSize;

extern float ***EBCBMCValue;
extern INT_VECTOR3 ***EBCBMCGridLocation;
extern int **SizeOfEBCBMCBins;
extern REAL *SumVDWGridRosenbluth; 
extern int NumberOfEBCBMCBins;
extern int EBCBMCVolume;
extern int UseEBCBMCEnergyFirstBead;

void MakeASCIGrid(void);


VECTOR MapToUnitCell(VECTOR pos);
VECTOR MapZToBox(VECTOR pos);
void MakeRigidFrameworkList(void);
void RigidFrameworkGrid(VECTOR pos,int typeA,REAL *Uvdw,REAL *Ucoul);

void MakeGrid(void);
void MakeDynamicGrid(void);
int WriteVDWGrid(int l);
void ReadVDWGrid(void);
int WriteCoulombGrid(void);
void ReadCoulombGrid(void);
REAL InterpolateVDWGrid(int typeA,VECTOR pos);
int InterpolateDelaunayGrid(int typeA,VECTOR pos);
int InterpolateDynamicGrid(VECTOR pos);
REAL InterpolateVDWForceGrid(int typeA,VECTOR pos,VECTOR *Force);
REAL InterpolateCoulombGrid(int typeA,VECTOR pos);
REAL InterpolateCoulombForceGrid(int typeA,VECTOR pos,VECTOR *Force);
void TestGrid(FILE *FilePtr);
void TestForceGrid(FILE *FilePtr);
INT_VECTOR3 ConvertXYZPositionToGridIndex(VECTOR pos);
VECTOR ConvertGridIndexToXYZIndex(INT_VECTOR3 GridIndex);
void BlockingVDWGrid(void);
void WriteRestartGrids(FILE *FilePtr);
void AllocateGridMemory(void);
void ReadRestartGrids(FILE *FilePtr);
void AdjustDynamicGrid(int GridValue, int typeA, VECTOR pos);
void PrintDynamicGrid(void);

//Delaunay Functions
int WriteDelaunayGrid(int l);
void ReadDelaunayGrid(void);
void ReadPolyhedraData(void);
int PointInPolygon(VECTOR pos,VECTOR Ray,int i, int j);
VECTOR PointProject(VECTOR PointToProject,int i, int j);
VECTOR CoordinateSystemShift(VECTOR PointToShift, int i, int j);
void PlaneFinder(int i, int j);
VECTOR PBCCorrection(VECTOR pos,int i, int j);
VECTOR RandomRay(void);
int CollinearVertex(VECTOR Vertex,VECTOR Ray, VECTOR pos);
int EdgeIntersection(VECTOR VertexOne, VECTOR VertexTwo,VECTOR Ray, VECTOR pos);
int LinePlaneIntersect(int i,int j, VECTOR pos, VECTOR Ray);
int PointInPolyhedra(int i,int Inaccessible,VECTOR pos);
void CreateEBCBMCProbBins(void);
void ReadEBCBMCBins(void);
int WriteEBCBMCBins(int l);
VECTOR ConvertCBMCGridPosToXYZ(VECTOR pos);




#endif
