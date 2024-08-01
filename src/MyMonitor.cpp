#include <ctime>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
using namespace std;

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "param.h"
#include "funcdef.h"

#undef __FUNCT__
#define __FUNCT__ "MyMonitor"
PetscErrorCode MyMonitor(SNES snes,PetscInt nt,PetscReal ftime,Vec X,void *ctx)
{
  PetscErrorCode ierr;
  AppCtx         *user  = (AppCtx*)ctx;
  MPI_Comm       comm;
  PetscMPIInt    rank;
  DM             da;
  Vec            localX;
  Field          **xx,**xn,**xn1;
  Fieldu         **uu;
  PetscInt       xs,xm,ys,ym,i,j,ir,ngnp,m,s0;
  PetscReal      t, tpass, fnorm;
  time_t         end_t;
  char           fname[150];
  fstream        btopf;

  SNESConvergedReason snes_reas=SNES_CONVERGED_FNORM_RELATIVE;

  PetscFunctionBeginUser;
  ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  ierr = SNESGetDM(snes,&da);;CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,localX,&xx);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->dau,user->U,&uu);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->Xn,&xn);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->Xn1,&xn1);CHKERRQ(ierr);

  ierr=DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  ngnp=0;

  if (user->ini==1) {
    //for (i=0; i<ntn; i++) advance_neutrals(user->dau,user->U,xx,sdt/t0);
  
//check if any negative temperatures (positive densities guaranteed)
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        if (xx[j][i].fx[3] <= 0.0) {
          cout<<"Negative or zero electron temperature Te = "<<xx[j][i].fx[3]*T0
              <<" at (i, j) = (" << i << ", " << j << ")"<<endl;
          ngnp=-1;
        }

        for (m=0; m<sl; m++) {
          s0=5*m;
          if (xx[j][i].fx[8+s0] <= 0.0) {
            cout<<"Negative or zero temperature of species "<< m <<" = "<<xx[j][i].fx[8+s0]*T0
            << " at (i, j) = (" << i << ", " << j << ")" << endl;
            ngnp=-1;
          }
          if (xx[j][i].fx[4+s0] <= 0.0) {
            cout << "Negative or zero density of species " << m << " = "
                 << xx[j][i].fx[4+s0] << " at (i, j) = (" << i << ", " << j << ")" << endl;
            ngnp=-2;
          }
        }
      }
    }
  }

  if (ngnp >= 0) {
    photoionization(uu,xs,xm,ys,ym);
    prod_loss_rates(xx,uu,xs,xm,ys,ym);
    collision_freq(xx,uu,xs,xm,ys,ym);
    ele_heating_rate(xx,uu,xs,xm,ys,ym);
    ele_cooling_rate(xx,uu,xs,xm,ys,ym);
    //neu_cooling_rate(uu,xs,xm,ys,ym);
    thermal_conductivity(xx,uu,xs,xm,ys,ym);
    jparallel(xx,uu,xs,xm,ys,ym);
  }

  t=ftime*t0;
  top_bc(t);

  if (user->ini == 1) {
    /* get converged/diverged information*/
    ierr = SNESGetConvergedReason(snes, &snes_reas);CHKERRQ(ierr);
    if (snes_reas < 0) cout << endl<< "SNES Failed with Reason "<< setw(3)<< snes_reas<< endl;
  }

  if (nvt!=-1 && ((t>=0.0 && t<=tsl+30.0) || (t>=TL && t<=TL+tsl+30.0)
      || (t>=2.0*(TL+tsl) && t<=2.0*(TL+tsl)+30.0)
      || (t>=3.0*(TL+tsl) && t<=3.0*(TL+tsl)+30.0))) user->nouta=nout;
  else user->nouta=20*nout;

  /* parallel output */
  if (nt % user->nouta == 0 || nt == ntot || nt % user->nosub == 0) {
    if (output_solution(da,xx,uu,t,nt,user->workdir)<0) return (-1);

/* output imposed v_theta or v_z at top boundary */
    if (nvt != -1) {
      strncpy(fname, user->workdir, 150);
      strcat(fname, "/output/top_bc.dat");
      if (nt==0) btopf.open(fname, fstream::out);
      else {
          btopf.open(fname);
          btopf.seekp(0,fstream::end);
      }
      if(!btopf) {
        cout << "Can't open file " << fname << endl;
        return -1;
      }

      if (smod == 0 && nt == 0) {
          if (nvt == 0) btopf<<"# Time (s)  vt(top) (m/s) on azimuthal grids"<<endl;
          if (nvt == 1) btopf<<"# Time (s)  vz(top) (m/s) on azimuthal grids"<<endl;
          btopf<<scientific<<setw(12)<<setprecision(3)<<0.0;

          for (j=0; j<=N2; j++) btopf<<scientific<<setw(12)<<setprecision(3)<<0.0;
          btopf<<endl;
      }
      else {
        //output imposed BC at top boundary
        btopf<<scientific<<setw(12)<<setprecision(3)<<t;
        for (j=0; j<=N2; j++) {
            if (btop[j]>-1.0e5) btopf<<scientific<<setw(12)<<setprecision(3)<<btop[j]*v0;
            else btopf<<scientific<<setw(12)<<setprecision(3)<<0.0;
        }
        btopf<<endl;
      }

      btopf.close();
    }
  }
  if (nt == ntot || (user->ini==1 && nt % user->nosub == 0)) {
    if (output_solution(da,xn,uu,t-dt*t0,nt-1,user->workdir) < 0) return -1;
    if (output_solution(da,xn1,uu,t-2.0*dt*t0,nt-2,user->workdir) < 0) return -1;
  }

  if (ngnp >=0 && snes_reas >=0 && user->ini == 1) {
    //smoothing
    filter(xx,xs,xm,ys,ym);

    //update arrays
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
          for (ir=0; ir<a3; ir++) {
              xn1[j][i].fx[ir]=xn[j][i].fx[ir];
              xn[j][i].fx[ir]=xx[j][i].fx[ir];
          }
      }
    }
  }

  ierr = DMDAVecRestoreArray(da,localX,&xx);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(user->dau,user->U,&uu);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->Xn,&xn);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->Xn1,&xn1);CHKERRQ(ierr);

  ierr = SNESGetFunctionNorm(snes, &fnorm);CHKERRQ(ierr);

  end_t=time(NULL);

  tpass=difftime(end_t, user->start_t);

  ierr = PetscPrintf(comm,"Time step %D Time elapsed %g (s) Function norm %g\n",
                     nt,tpass,fnorm);CHKERRQ(ierr);

  //stop execution if snes didn't converge or negative desnity/temperatures
  if (ngnp < 0 || snes_reas < 0) MPI_Abort(PETSC_COMM_WORLD,-1);

  PetscFunctionReturn(0);
}
