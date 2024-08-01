/************************************************************************
 Routine: euvflux_seg
 Purpose: get solar EUV flux, neutrla absorption and ionization cross
          sections

 Input:  workdir - working directory
         f107, f107a - solar 10.7 cm emission flux & its 81-day average
 output: euvflux[37]    - EUV flux at 37 wavelength rages
         segabs[37][5]  - EUV absorption cross sections of N, O, He, N2, O2
         segabs[37][5]  - Ionization cross sections of N, O, He, N2, O2

 By Jiannan Tu
 * 2/10/2011, 12/9/2013
************************************************************************/
#include <iostream>
#include <fstream>
#include <string.h>
//#include <iomanip>
using namespace std;

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"

#define tag1 1

int euvflux_seg(DM da,char *workdir)
{
    double f74, ai, pp, xflux;
    int    j, i, prid, numpr, freef=0;
    char   fname[150];
    MPI_Comm comm;
    MPI_Status status;

    fstream euvfstr, absfstr, ionfstr;

    PetscObjectGetComm((PetscObject)da,&comm);
    MPI_Comm_rank(comm,&prid);
    MPI_Comm_size(comm,&numpr);

    if (prid==0) freef=1;
    else {
      MPI_Recv(&freef,1,MPI_INT,prid-1,tag1,PETSC_COMM_WORLD,&status);
    }

    if(freef==1) {
      strncpy(fname, workdir, 150);
      strcat(fname, "/inp/euvflux.inp");

      euvfstr.open(fname, fstream::in);
      if(!euvfstr) {
        cout << "Can't open file " <<fname << endl;
        return -1;
      }

      //EUV flux from EUVAC model
      pp=0.50*(f107+f107a);

      for (j = 0; j < 37; j++) {
        euvfstr >> f74 >> ai >> pene[j];
        euvfstr.ignore(100, '\n');

        xflux=1.0+ai*(pp-80.0);
        if (xflux < 0.8) xflux=0.8;

        euvflux[j]=f74*xflux*1.0e13;  //euvflux in unit of photons m^-2 s^-1
      }
      euvfstr.close();

      //solar radiation scattered through atmosphere into night sector
      euvfluxn[0]=1.0e10;   //He II     303.78 A
      euvfluxn[1]=0.3e10;   //He I      584.33 A
      euvfluxn[2]=3.0e10;   //Ly-beta  1025.70 A
      euvfluxn[3]=1.5e13;   //Ly-alpha 1215.67 A

      //----photo-absorption cross sections
      strncpy(fname, workdir, 150);
      strcat(fname, "/inp/phabsdt.inp");

      absfstr.open(fname, fstream::in);
      if(!absfstr) {
        cout << "Can't open file " <<fname << endl;
        return -1;
      }

/* ----- absorption cross sections of O, O2, N2, N & He*/
      for (j = 0; j < 37; j++) {
        absfstr >> i >> segabs[j][0] >> segabs[j][1] >> segabs[j][2]
                >> segabs[j][3] >> segabs[j][4];
        absfstr.ignore(70, '\n');

        //in unit of m^2
        for (i = 0; i < 5; i++) segabs[j][i]=segabs[j][i]*1.0e-22;

        //test
        /*cout << segabs[j][0] << "  " << segabs[j][1] << "  " << segabs[j][2]
             << "  " << segabs[j][3] << "  " << segabs[j][4] << endl; */
      }

      //nighttime absorption cross sections for O, O2, N2, NO
      absfstr.ignore(120, '\n');
      for (j = 0; j < 4; j++) {
        absfstr >> i >> segabsn[j][0] >> segabsn[j][1] >> segabsn[j][2]
                >> segabsn[j][3];
        absfstr.ignore(70, '\n');

        //in unit of m^2
        for (i = 0; i < 4; i++) segabsn[j][i]=segabsn[j][i]*1.0e-22;
      }
      absfstr.close();

      //----photoionization cross sections
      strncpy(fname, workdir, 150);
      strcat(fname, "/inp/phiondt.inp");

      ionfstr.open(fname, fstream::in);
      if(!ionfstr) {
        cout << "Can't open file " <<fname << endl;
        return -1;
      }

/* ionization cross sections of following 4 photoionization processes
 *   (0) O + hv --> O+ + e, (1) O2 + hv --> O+ + O + e, (2) O2 + hv --> O2+ + e, 
 *   (3) N2 + hv --> N2+ + e, (4) N2+hv --> N+ + N + e, (5) N + hv --> N+  + e
     (6) He + hv --> He+ + e */
      for (j = 0; j < 37; j++) {
        ionfstr >> i >> segion[j][0] >> segion[j][1] >> segion[j][2]
                >> segion[j][3] >> segion[j][4] >> segion[j][5] >> segion[j][6];
        ionfstr.ignore(70, '\n');

        //in unit of m^2
        for (i = 0; i < 7; i++) segion[j][i]=segion[j][i]*1.0e-22;

        //test
        /*cout << setw(8) << setprecision(4) << segion[j][0]/1.0e-22
             << setw(8) << setprecision(4) << segion[j][1]/1.0e-22
             << setw(8) << setprecision(4) << segion[j][2]/1.0e-22
             << setw(8) << setprecision(4) << segion[j][3]/1.0e-22
             << setw(8) << setprecision(4) << segion[j][4]/1.0e-22
             << setw(8) << setprecision(4) << segion[j][5]/1.0e-22
             << setw(8) << setprecision(4) << segion[j][6]/1.0e-22 << endl;*/
      }

      //nighttime ionization cross sections for O, O2, N2, NO
      ionfstr.ignore(120, '\n');
      for (j = 0; j < 4; j++) {
        ionfstr >> i >> segionn[j][0] >> segionn[j][1] >> segionn[j][2]
                >> segionn[j][3];
        ionfstr.ignore(120, '\n');

        //in unit of m^2
        for (i = 0; i < 4; i++) segionn[j][i]=segionn[j][i]*1.0e-22;
      }
      ionfstr.close();
    }

    if(prid != numpr-1) MPI_Send(&freef,1,MPI_INT,prid+1,tag1,PETSC_COMM_WORLD);

    return 0;
}