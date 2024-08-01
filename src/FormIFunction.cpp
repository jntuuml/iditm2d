#include <ctime>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "param.h"
#include "funcdef.h"

#undef __FUNCT__
#define __FUNCT__ "FormIFunction"
PetscScalar function_part(int,int,int,int,int,Field**,Field**,Field**,Fieldu**);

PetscErrorCode FormIFunction(SNES snes,Vec X,Vec F,void *ctx)
{
  AppCtx         *user = (AppCtx*)ctx;       /* user-defined application context */
  DM             da,daU=user->dau;
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,yj,xi,ir;
  Vec            localX,localU,localXn,localXn1;
  /* structures that contain variables of interest and left hand side of governing equations */
  Field          **xx,**ff,**xn,**xn1;
  Fieldu         **uu;

  PetscFunctionBeginUser;
  ierr = SNESGetDM(snes,&da);;CHKERRQ(ierr);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DAGlobalToLocalBegin(),DAGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&localXn);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,user->Xn,INSERT_VALUES,localXn);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,user->Xn,INSERT_VALUES,localXn);CHKERRQ(ierr);

  ierr = DMGetLocalVector(da,&localXn1);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,user->Xn1,INSERT_VALUES,localXn1);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,user->Xn1,INSERT_VALUES,localXn1);CHKERRQ(ierr);

  ierr = DMGetLocalVector(daU,&localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(daU,user->U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(daU,user->U,INSERT_VALUES,localU);CHKERRQ(ierr);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArray(da,F,&ff);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,localX,&xx);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(daU,localU,&uu);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,localXn,&xn);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,localXn1,&xn1);CHKERRQ(ierr);

  /* Get local grid boundaries */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);
 
    /* Compute function over the locally owned part of the grid */
    for (j=ys; j<ys+ym; j++) {
        yj=j-ys;
        for (i=xs; i<xs+xm; i++) {
            xi=i-xs;
            for (ir=0; ir<a3; ir++) {
                ff[j][i].fx[ir]=function_part(j,i,ir,yj,xi,xx,xn,xn1,uu);
            }
        }
    }

  /* Restore vectors */
  ierr = DMDAVecRestoreArray(da,F,&ff);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,localX,&xx);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(daU,localU,&uu);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,localXn,&xn);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,localXn1,&xn1);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(daU,&localU);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localXn);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localXn1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
