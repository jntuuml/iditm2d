/*****************************************************************************
 *  output_final
 *   Output simulation results for the last time step. 
 *
 *   Jiannan Tu
 *   12/11/2013, 1/11/2014
 ****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <cmath>
#include <unistd.h>

using namespace std;

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"

#define master 0
#define tag1  1

PetscInt output_solution(DM da,Field **xx,Fieldu **uu,double t,PetscInt nt,char *workdir)
{
  int        i,j,l,s0,xs,xm,ys,ym,ierr,file_free,yj,xi;
  char       fname[150],num[60];
  fstream    outfstr,colfstr; //pdlfstr;
  PetscInt   prid,numpr,jm,im;
  double     ve1,ve2,ve3,ns,ne,J1,J2,J3;
  MPI_Comm   comm;
  MPI_Status status;
  const double bfc=-1.0;

  ierr=PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);

  MPI_Comm_rank(comm,&prid);
  MPI_Comm_size(comm,&numpr);

  /*
   Get local grid boundaries
  */
  ierr=DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  file_free=0;

  /* the first process receive results from other processes and does the output */
  if (prid==0) {
    file_free=1;

    /* open files to write */
    sprintf(num, "%d", nt);
    strncpy(fname, workdir, 150);
    strcat(fname, "/output/uvbenp");
    strcat(fname, num);
    strcat(fname, ".dat");
    outfstr.open(fname, fstream::out);
    if(!outfstr) {
      cout << "Can't open file " << fname << endl;
      return -1;
    }

    strncpy(fname, workdir, 150);
    strcat(fname, "/output/collision");
    strcat(fname, num);
    strcat(fname, ".dat");
    colfstr.open(fname, fstream::out);
    if(!colfstr) {
      cout << "Can't open file " << fname << endl;
      return -1;
    }

    /*strncpy(fname, workdir, 150);
    strcat(fname, "/output/prdlos");
    strcat(fname, num);
    strcat(fname, ".dat");
    pdlfstr.open(fname, fstream::out);
    if(!pdlfstr) {
      cout << "Can't open file " << fname << endl;
      return -1;
    }*/
  }
  else {
    MPI_Recv(&file_free,1,MPI_INT,prid-1,tag1,comm,&status);

    /* open the files to append contents*/
    sprintf(num, "%d", nt);
    strncpy(fname, workdir, 150);
    strcat(fname, "/output/uvbenp");
    strcat(fname, num);
    strcat(fname, ".dat");
    outfstr.open(fname);
    if(!outfstr) {
      cout << "Can't open file " << fname << endl;
      return -1;
    }
    outfstr.seekp(0,fstream::end);

    strncpy(fname, workdir, 150);
    strcat(fname, "/output/collision");
    strcat(fname, num);
    strcat(fname, ".dat");
    colfstr.open(fname);
    if(!colfstr) {
        cout << "Can't open file " << fname << endl;
        return -1;
    }
    colfstr.seekp(0,fstream::end);

    /*strncpy(fname, workdir, 150);
    strcat(fname, "/output/prdlos");
    strcat(fname, num);
    strcat(fname, ".dat");
    pdlfstr.open(fname);
    if(!pdlfstr) {
        cout << "Can't open file " << fname << endl;
        return -1;
    }
    pdlfstr.seekp(0,fstream::end);*/
  }

  if(file_free==1) {
    if (prid==0) outfstr<<"# t ="<<scientific<<setw(12)<<setprecision(3)<<t<<endl;

    for (j=ys; j<ys+ym; j++) {
      yj=j-ys; jm=j-1;

      outfstr << "# j = " << setw(3) << j << endl;

      for (i=xs; i<xs+xm; i++) {
        xi=i-xs; im=i-1;

        outfstr <<scientific<<setw(14)<<setprecision(5)<<zh[i]*1.0e-3;

        for (l=0; l<3; l++) outfstr<<scientific<<setw(14)<<setprecision(5)<<xx[j][i].fx[l]*B0;

        ve1=0.0; ve2=0.0; ve3=0.0; ne=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ns=xx[j][i].fx[4+s0];
            ne=ne+ns;
            ve1=ve1+ns*xx[j][i].fx[5+s0];
            ve2=ve2+ns*xx[j][i].fx[6+s0];
            ve3=ve3+ns*xx[j][i].fx[7+s0];
        }
        if (i==0) {
            J1=0.5*(1.0+bfc)*(xx[j][i].fx[2]-xx[jm][i].fx[2])/h2[i];
            J2=0.0; //0.5*(bfc-1.0)*(xx[j][i].fx[2]+xx[jm][i].fx[2])/h1[i];
            J3=0.0; //(1.0-bfc)*xx[j][i].fx[1]/h1[i]+0.5*(1.0+bfc)*xx[j][i].fx[1]/rr[i]
               //-(xx[j][i].fx[0]-xx[jm][i].fx[0])/h2[i];
        }
        else {
            J1=0.5*( xx[j][i].fx[2]+xx[j][im].fx[2]-xx[jm][i].fx[2]
                    -xx[jm][im].fx[2])/h2[i];
            J2=0.5*( xx[j][im].fx[2]+xx[jm][im].fx[2]-xx[j][i].fx[2]
                    -xx[jm][i].fx[2])/h1[i];
            J3= (xx[j][i].fx[1]-xx[j][im].fx[1])/h1[i]
               +0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])/rr[i]
               -(xx[j][i].fx[0]-xx[jm][i].fx[0])/h2[i];
        }
        ve1=ve1-J1/e;
        ve2=ve2-J2/e;
        ve3=ve3-J3/e;
        outfstr<<scientific<<setw(14)<<setprecision(5)<<ve1/ne*v0;
        outfstr<<scientific<<setw(14)<<setprecision(5)<<ve2/ne*v0;
        outfstr<<scientific<<setw(14)<<setprecision(5)<<ve3/ne*v0;
        outfstr<<scientific<<setw(14)<<setprecision(5)<<xx[j][i].fx[3]*T0;

/* O+ density, 3-comp velocity, temperature, followed by for O2+, N2+, H+, NO+ 
 * O, O2, N2, & H. Density in m^-3, velocity in m/s, temperature in K */
        for (l=0; l<sl; l++) {
          s0=5*l;

          outfstr <<scientific<<setw(14)<<setprecision(5)<<xx[j][i].fx[4+s0]*n0
                  <<scientific<<setw(14)<<setprecision(5)<<xx[j][i].fx[5+s0]*v0
                  <<scientific<<setw(14)<<setprecision(5)<<xx[j][i].fx[6+s0]*v0
                  <<scientific<<setw(14)<<setprecision(5)<<xx[j][i].fx[7+s0]*v0
                  <<scientific<<setw(14)<<setprecision(5)<<xx[j][i].fx[8+s0]*T0;
        }

        for (l=0; l<sm; l++) {
          s0=5*l;

          outfstr <<scientific<<setw(14)<<setprecision(5)<<exp(uu[j][i].fu[s0])*n0
                  <<scientific<<setw(14)<<setprecision(5)<<uu[j][i].fu[1+s0]*v0
                  <<scientific<<setw(14)<<setprecision(5)<<uu[j][i].fu[2+s0]*v0
                  <<scientific<<setw(14)<<setprecision(5)<<uu[j][i].fu[3+s0]*v0
                  <<scientific<<setw(14)<<setprecision(5)<<uu[j][i].fu[4+s0]*T0;
        }
/* NO & N densities */
        outfstr <<scientific<<setw(14)<<setprecision(5)<<exp(uu[j][i].fu[20])*n0
                <<scientific<<setw(14)<<setprecision(5)<<exp(uu[j][i].fu[21])*n0
                <<endl;

/* collision frequencies */
        colfstr << scientific << setw(9) << setprecision(3) << zh[i]*1.0e-3;
        for (l=0; l<90; l++) colfstr << scientific << setw(11) << setprecision(3)
                                       << uu[j][i].nu[l]/t0;
        colfstr << endl;

/* production rates and loss coefficients */
        /*pdlfstr << scientific << setw(9) << setprecision(3) << zh[i]*1.0e-3;
        pdlfstr << scientific << setw(11) << setprecision(3) << uu[j][i].lamdae*lamda0;
        pdlfstr << scientific << setw(11) << setprecision(3) << uu[j][i].betae*beta0;
        pdlfstr << scientific << setw(11) << setprecision(3) << Qe[yj][xi]*t0/p0;
        pdlfstr << scientific << setw(11) << setprecision(3) << Ce[yj][xi]*t0/p0;
        for (l = 0; l < slm; l++)
            pdlfstr << scientific << setw(11) << setprecision(3) << uu[j][i].lamdas[l]*lamda0;
        for (l=0; l<slm+2; l++) {
          pdlfstr<<scientific<<setw(11)<<setprecision(3)<<Ps[yj][xi][l]*n0/t0*1.0e-6;
          pdlfstr<<scientific<<setw(11)<<setprecision(3)<<Ls[yj][xi][l]/t0;
        }
        pdlfstr << endl;*/
      }
    }

    outfstr.close();
    colfstr.close();
    //pdlfstr.close();
  }

  if(prid != numpr-1) MPI_Send(&file_free,1,MPI_INT,prid+1,tag1,comm);
    
  return 0;
}
