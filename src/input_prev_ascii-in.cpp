/*****************************************************************************
 *  input_prev
 *   input parameters & simulation results from a previous simulation. 
 *
 *   Jiannan Tu
 *   12/11/2013, 1/29/2014
 ****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <cmath>

using namespace std;

#include "param.h"
#include "funcdef.h"

#define tag1 1

int input_prev(DM da,Field **xx,Fieldu **uu,char *workdir,char *prefln)
{
  int        i,j,l,s0,d1=80,ip=99;
  double     f20[55];
  char       fname[150];
  fstream    infstr;
  int        ierr,prid,numpr,xs,xm,ys,ym,file_free=0;
  MPI_Comm   comm;
  MPI_Status status;
  long       pos;
  double     *zz,**ab,sc;

  zz=new double[d1];
  ab=new double*[a3+22];
  for (l=0; l<a3+22; l++) ab[l]=new double[d1];

  for (i=0; i<d1; i++) {
      if (zh[i]<=2000.0e3) ip=i;
  }

  ierr=PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);

  MPI_Comm_rank(comm,&prid);
  MPI_Comm_size(comm,&numpr);

  ierr=DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  if (prid==0) file_free=1;
  else {
    MPI_Recv(&file_free,1,MPI_INT,prid-1,tag1,comm,&status);
    MPI_Recv(&pos,1,MPI_LONG,prid-1,tag1,comm,&status);
  }

  if (file_free==1) {
    strncpy(fname, workdir, 150);
    strcat(fname, "/inp/");
    strcat(fname, prefln);

    infstr.open(fname, fstream::in);
    if(!infstr) {
          cout << "Can't open file " << fname << endl;
          return -1;
    }

    if (prid==0) {
        infstr.ignore(60,'\n');
    }
    else infstr.seekp(pos+1); //move to the position after just read 

    for (j = ys; j < ys+ym; j++) {
        infstr.ignore(60, '\n');

        for (i = xs; i < d1; i++) {
            for (l = 0; l < 55; l++) infstr >> f20[l];
            infstr.ignore(1,'\n');

            zz[i]=f20[0]*1.0e3;

            for (l = 0; l < 3; l++) {
                //Normalized perturbation magnetic field components
                ab[l][i]=f20[l+1]/B0;
            }

            //electron temperature (un-normalized)
            ab[3][i]=f20[7]/T0;

            //normalized density, velocity, and temperature
            for (l = 0; l < slm; l++) {
                s0=5*l;
                ab[4+s0][i]=log(f20[8+s0]/n0);
                ab[5+s0][i]=f20[9+s0]/v0;
                ab[6+s0][i]=f20[10+s0]/v0;
                ab[7+s0][i]=f20[11+s0]/v0;
                ab[8+s0][i]=f20[12+s0]/T0;
            }

            //natural logarithm of normalized neutral NO and N densities
            ab[49][i]=log(f20[8+5*(slm)]/n0);
            ab[50][i]=log(f20[9+5*(slm)]/n0);
        }

        for (i=0; i<=ip; i++) {
            for (l=0; l<3; l++) {
                xx[j][i].fx[l]=lgrg(zz,ab[l],d1,zh[i]);
            }
            xx[j][i].fx[3]=lgrg(zz,ab[3],d1,zh[i]);
            for (l=0; l<sl; l++) {
                s0=5*l;
                xx[j][i].fx[4+s0]=exp(lgrg(zz,ab[4+s0],d1,zh[i]));                
                xx[j][i].fx[5+s0]=lgrg(zz,ab[5+s0],d1,zh[i]);
                xx[j][i].fx[6+s0]=lgrg(zz,ab[6+s0],d1,zh[i]);
                xx[j][i].fx[7+s0]=lgrg(zz,ab[7+s0],d1,zh[i]);
                xx[j][i].fx[8+s0]=lgrg(zz,ab[8+s0],d1,zh[i]);
            }

            //natural logarithm of normalized neutral density, velocity, and temperature
            for (l = 0; l < sm; l++) {
                s0=5*(l+sl);
                uu[j][i].fu[5*l]=  lgrg(zz,ab[4+s0],d1,zh[i]);
                uu[j][i].fu[1+5*l]=lgrg(zz,ab[5+s0],d1,zh[i]);
                uu[j][i].fu[2+5*l]=lgrg(zz,ab[6+s0],d1,zh[i]);
                uu[j][i].fu[3+5*l]=lgrg(zz,ab[7+s0],d1,zh[i]);
                uu[j][i].fu[4+5*l]=lgrg(zz,ab[8+s0],d1,zh[i]);
            }
            //natural logarithm of normalized neutral NO and N densities
            uu[j][i].fu[20]=lgrg(zz,ab[49],d1,zh[i]);
            uu[j][i].fu[21]=lgrg(zz,ab[50],d1,zh[i]);
        }
        for (i=ip+1; i<xs+xm; i++) {
            for (l=0; l<3; l++) xx[j][i].fx[l]=xx[j][i-1].fx[l]*exp(-20.0*(rr[i]/rr[i-1]-1.0));
            xx[j][i].fx[3]=xx[j][i-1].fx[3];
            for (l=0; l<sl; l++) {
                s0=5*l;
                if (l==0) sc=20.0;
                else if (l==1) sc=5.0;
                else sc=65.0;
                if (l!=1) xx[j][i].fx[4+s0]=xx[j][i-1].fx[4+s0]*exp(-sc*(rr[i]/rr[i-1]-1.0));
                xx[j][i].fx[4+s0]=xx[j][i-1].fx[4+s0]*exp(-sc*(rr[i]/rr[i-1]-1.0));
                xx[j][i].fx[5+s0]=xx[j][i-1].fx[5+s0];
                xx[j][i].fx[6+s0]=xx[j][i-1].fx[6+s0];
                xx[j][i].fx[7+s0]=xx[j][i-1].fx[7+s0];
                xx[j][i].fx[8+s0]=xx[j][i-1].fx[8+s0];
            }
            for (l=0; l<sm; l++) {
                s0=5*l;
                if (l==0) sc=50.0;
                else if (l==1 || l==2) sc=70.0;
                else sc=2.0;
                uu[j][i].fu[s0]=exp(uu[j][i-1].fu[s0])*exp(-sc*(rr[i]/rr[i-1]-1.0));
                uu[j][i].fu[s0]=log(uu[j][i].fu[s0]);
                uu[j][i].fu[1+s0]=uu[j][i-1].fu[1+s0];
                uu[j][i].fu[2+s0]=uu[j][i-1].fu[2+s0];
                uu[j][i].fu[3+s0]=uu[j][i-1].fu[3+s0];
                uu[j][i].fu[4+s0]=uu[j][i-1].fu[4+s0];
            }
            uu[j][i].fu[20]=exp(uu[j][i-1].fu[20])*exp(-65.0*(rr[i]/rr[i-1]-1.0));
            uu[j][i].fu[20]=log(uu[j][i].fu[20]);
            uu[j][i].fu[21]=exp(uu[j][i-1].fu[21])*exp(-50.0*(rr[i]/rr[i-1]-1.0));
            uu[j][i].fu[21]=log(uu[j][i].fu[21]);
        }
    }

    pos=infstr.tellp();
    infstr.close();

    delete[] zz;
    for (l=0; l<a3; l++) delete[] ab[l];
    delete[] ab;

    /* 
     * tell the next processor that the file is free, and the position of the 
     * last character that has been read
     */
    if(prid != numpr-1) {
      MPI_Send(&file_free,1,MPI_INT,prid+1,tag1,comm);
      MPI_Send(&pos,1,MPI_LONG,prid+1,tag1,comm);
    }
  }
 
  return 0;
}
