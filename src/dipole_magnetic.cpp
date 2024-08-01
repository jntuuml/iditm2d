/***************************************************************************** 
 *   dipole_magnetic
 *   Calculate background magnetic field components with magnetic meridian in
 *   curvilinear coordinates
 *
 *   Jiannan Tu
 *   12/9/2013, 1/10/2014
 *****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <cmath>
using namespace std;

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"
#include "funcdef.h"

#define tag1  1

int dipole_magnetic(DM da,Vec U,char *workdir)
{
    int    i, j, prid, numpr, xs, xm, ys, ym, file_free=0, ierr;
    double ar3, BB;
    const double Beq=3.12e-5; //magnetic induction at equator on the earth's surface
    char fname[150];
    fstream divbfs;
    MPI_Comm comm;
    MPI_Status status;
    Fieldu **uu;

    PetscObjectGetComm((PetscObject)da,&comm);
    MPI_Comm_rank(comm,&prid);
    MPI_Comm_size(comm,&numpr);

    /*
     Get pointers to vector data
    */
    ierr = DMDAVecGetArray(da,U,&uu);CHKERRQ(ierr);

    DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);

    BB=Beq/B0;
    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            /* normalized dipole magnetic field in curvilinear coordinates.
             * here theta is polar angle (0 - 2*pi, not 0 - pi). It is relative
             * to geomagnetic north pole (assuming as x-axis). */
            ar3=pow(Re/(rr[i]*r0), 3.0)*BB;
            uu[j][i].B0[0]=-2.0*ar3*cos(theta[j]);
            uu[j][i].B0[1]=-ar3*sin(theta[j]);
        }
    }

    strncpy(fname, workdir, 150);
    strcat(fname, "/output/B0divB0.dat");

    if (prid==0) {
      file_free=1;

      divbfs.open(fname, fstream::out);
      if(!divbfs) {
        cout << "Can't open file " << fname << endl;
        MPI_Abort(comm,-1);
      }
    }
    else {
      MPI_Recv(&file_free,1,MPI_INT,prid-1,tag1,comm,&status);

      /* open the file to append contents*/
      divbfs.open(fname);
      if(!divbfs) {
        cout << "Can't open file " << fname << endl;
        MPI_Abort(comm,-1);
      }
      divbfs.seekp(0,fstream::end);
    }

    if(file_free==1) {
      for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
          divbfs<<scientific<<setw(12)<<setprecision(4)<<(rr[i]*r0-Re)/1.0e3
                <<scientific<<setw(14)<<setprecision(5)<<uu[j][i].B0[0]*B0
                <<scientific<<setw(14)<<setprecision(5)<<uu[j][i].B0[1]*B0<<endl;
        }
      }
      divbfs.close();
    }

    if(prid != numpr-1) MPI_Send(&file_free,1,MPI_INT,prid+1,tag1,comm);

    ierr = DMDAVecRestoreArray(da,U,&uu);CHKERRQ(ierr);

    return (0);
}
