/******************************************************************************
 *  input_param
 *   Input various parameters
 *
 * Jiannan Tu 9/13/2013, 2/4/2014
 ******************************************************************************/
#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>

using namespace std;

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"

#define tag1 1

int input_param(char *workdir,AppCtx *user)
{
  char fname[150];
  int  prid, numpr, file_free=0;
  MPI_Status status;

  strncpy(fname, workdir, 150);
  strcat(fname, "/iditm2d.dat");

  MPI_Comm_rank(PETSC_COMM_WORLD,&prid);
  MPI_Comm_size(PETSC_COMM_WORLD,&numpr);

  if (prid==0) file_free=1;
  else {
    MPI_Recv(&file_free,1,MPI_INT,prid-1,tag1,PETSC_COMM_WORLD,&status);
  }

  if (file_free==1) {
    fstream infstr(fname, fstream::in);
    if(!infstr) {
      cout << "Can't open file " << fname << endl;
      return -1;
    }
    infstr >> smod;
    infstr.ignore(200, '\n');
    infstr >> iyr;
    infstr.ignore(200, '\n');
    infstr >> mon;
    infstr.ignore(200, '\n');
    infstr >> idate;
    infstr.ignore(200, '\n');
    infstr >> sec;
    infstr.ignore(200, '\n');
    infstr >> LT;
    infstr.ignore(200, '\n');
    infstr >> f107;
    infstr.ignore(200, '\n');
    infstr >> f107a;
    infstr.ignore(200, '\n');
    infstr >> Ap;
    infstr.ignore(200, '\n');
    infstr >> r0;
    infstr.ignore(200, '\n');
    infstr >> n0;
    infstr.ignore(200, '\n');
    infstr >> B0;
    infstr.ignore(200, '\n');
    infstr >> ntot;
    infstr.ignore(200, '\n');
    infstr >> nout;
    infstr.ignore(200, '\n');
    infstr >> user->nosub;
    infstr.ignore(200, '\n');
    infstr >> npre;
    infstr.ignore(200, '\n');
    infstr >> dt;
    infstr.ignore(200, '\n');
    infstr >> rb;
    infstr.ignore(200, '\n');
    infstr >> ru;
    infstr.ignore(200, '\n');
    infstr >> a1;
    infstr.ignore(200, '\n');
    infstr >> a2;
    infstr.ignore(200, '\n');
    infstr >> sl;
    infstr.ignore(200, '\n');
    infstr >> sm;
    infstr.ignore(200, '\n');
    infstr >> user->prefln;
    infstr.ignore(200, '\n');
    infstr >> nvt;
    infstr.ignore(200, '\n');
    infstr >> bimp;
    infstr.ignore(200, '\n');
    infstr >> bimpf;
    infstr.ignore(200, '\n');
    infstr >> tsl;
    infstr.ignore(200, '\n');
    infstr >> TL;
    infstr.ignore(200, '\n');
    infstr >> lat[0] >> lat[1];
    infstr.ignore(200, '\n');
    infstr >> lat[2] >> lat[3];
    infstr.ignore(200, '\n');
    infstr >> d;
    if (smod==1) {
        infstr.ignore(200, '\n');
        infstr >> user->preflm;
        infstr.ignore(200, '\n');
        infstr >> user->preflo;
    }

    infstr.close();

    if(a1 % 5 != 0 || a2 % 5 !=0) {
      cout << "Number of grids along x1 and x2 must be multiples of 5" << endl;
      return -1;
    }

    user->ti=0.0;
    user->nt=0;
  }

  if(prid != numpr-1) MPI_Send(&file_free,1,MPI_INT,prid+1,tag1,PETSC_COMM_WORLD);

  slm=sl+sm;

  //# of state variables on a gird
  a3=5*sl+4;

  //# of grids along x1 curvilinear coordinate
  N1=a1-1;

  //# of grids along x2 curvilinear coordinate
  N2=a2-1;

  aa=a3*a2;

  dth=360.0/double(a2);

  //convert boundary heights to radial distances (meters))
  rb=1.0e3*rb+Re;
  ru=1.0e3*ru+Re;

  return 0;
}
