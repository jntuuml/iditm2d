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
#define __FUNCT__ "FormIJacobian"
PetscScalar function_part(int,int,int,int,int,Field**,Field**,Field**,Fieldu**);

PetscErrorCode FormIJacobian(SNES snes,Vec X,Mat Jac,Mat Jpre,void *ctx)
{
    AppCtx         *user = (AppCtx*)ctx;       /* user-defined application context */
    DM             da,daU;
    PetscErrorCode ierr;
    PetscInt       i,j,xs,ys,xm,ym,ir;
    Vec            localX,localXn,localXn1,localU;
    Field          **xx,**xn,**xn1;
    Fieldu         **uu;
    MatStencil     row,*col;
    int            l,m,jp,jm,ip,im,im2,nv,yj,xi,jj,ii,cc;
    double         *vals;
    int            s0,s4,s5,s6,s7,s8;
    double         ffm, ffp, absx, segm, x1;
    double         ne[6], cof, ve1[6], ve2[6], ve3[6], B01p[2], ns, vep[4];
    double         botn, botp;
    const double   dx=1.5e-8;
    const int      nzer=13+16*sl;
    const double   one6th=0.1666666666666667;
    const double   nc=1.0e8/n0;

    PetscFunctionBeginUser;
    ierr = SNESGetDM(snes,&da);CHKERRQ(ierr);

    daU=user->dau;
    vals=new double[nzer];
    col=new MatStencil[nzer];

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
    ierr = DMDAVecGetArray(da,localX,&xx);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(daU,localU,&uu);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,localXn,&xn);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,localXn1,&xn1);CHKERRQ(ierr);

    /* Get local grid boundaries */
    ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

    for (j=ys; j<ys+ym; j++) {
        jm=j-1; jp=j+1; yj=j-ys;

        row.j=j;

        for (i=xs; i<xs+xm; i++) {
            im=i-1; im2=i-2; ip=i+1; xi=i-xs;
            row.i=i;

            ne[0]=0.0;
            ve1[0]=0.0; ve2[0]=0.0; ve3[0]=0.0;
            for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
                ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
                ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
            }
            ve1[0]=ve1[0]/ne[0];
            ve2[0]=ve2[0]/ne[0];
            ve3[0]=ve3[0]/ne[0];

            for (ir=0; ir<a3; ir++) {
                row.c=ir; nv=0;

                if (ir==0) {
                    ne[2]=0.0; ve1[2]=0.0; ve2[2]=0.0;
                    for (l=0; l<sl; l++) {
                        s0=5*l;
                        ns=xx[jp][i].fx[4+s0];
                        ne[2]=ne[2]+ns;
                        ve1[2]=ve1[2]+ns*xx[jp][i].fx[5+s0];
                        ve2[2]=ve2[2]+ns*xx[jp][i].fx[6+s0];
                    }
                    ve1[2]=ve1[2]/ne[2];
                    ve2[2]=ve2[2]/ne[2];

                    if (i == 0) {
                        botn=-xx[j][i].fx[1];
                        botp=-xx[jp][i].fx[1];
                    }
                    else {
                        botn=xx[j][im].fx[1];
                        botp=xx[jp][im].fx[1];
                    }

                    cof=0.5*dt2/h2[i];

                    if (i>0) {
                        col[nv].j=j;  col[nv].i=im; col[nv].c=1;
                        vals[nv++]=cof*ve1[0];

                        col[nv].j=jp; col[nv].i=im; col[nv].c=1;
                        vals[nv++]=-cof*ve1[2];

                        col[nv].j=j;  col[nv].i=i;  col[nv].c=1;
                        vals[nv++]=cof*ve1[0];

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=1;
                        vals[nv++]=-cof*ve1[2];
                    }
                    col[nv].j=jm; col[nv].i=i;  col[nv].c=0;
                    vals[nv++]=-cof*ve2[0];

                    col[nv].j=j;  col[nv].i=i;  col[nv].c=0;
                    vals[nv++]=3.0+cof*(ve2[2]-ve2[0]);

                    col[nv].j=jp; col[nv].i=i;  col[nv].c=0;
                    vals[nv++]=cof*ve2[2];

                    cof=2.0*cof;

                    for (l=0; l<sl; l++) {
                        s0=5*l;
                        col[nv].j=j;  col[nv].i=i;  col[nv].c=4+s0;
                        vals[nv++]=cof/ne[0]
                                      *( (ve2[0]-xx[j][i].fx[6+s0])
                                         *( 0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])
                                           +uu[j][i].B0[0])
                                        -(ve1[0]-xx[j][i].fx[5+s0])
                                         *(0.5*(xx[j][i].fx[1]+botn)+uu[j][i].B0[1]));
                                      
                        col[nv].j=j;  col[nv].i=i;  col[nv].c=5+s0;
                        vals[nv++]=cof/ne[0]*xx[j][i].fx[4+s0]
                                      *(0.5*(xx[j][i].fx[1]+botn)+uu[j][i].B0[1]);
                                      
                        col[nv].j=j;  col[nv].i=i;  col[nv].c=6+s0;
                        vals[nv++]=-cof/ne[0]*xx[j][i].fx[4+s0]
                                       *(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]);

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=4+s0;
                        vals[nv++]=cof/ne[2]
                                      *( (ve1[2]-xx[jp][i].fx[5+s0])
                                         *(0.5*(xx[jp][i].fx[1]+botp)+uu[jp][i].B0[1])
                                        -(ve2[2]-xx[jp][i].fx[6+s0])
                                         *(0.5*(xx[jp][i].fx[0]+xx[j][i].fx[0])+uu[jp][i].B0[0]));

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=5+s0;
                        vals[nv++]=-cof/ne[2]*xx[jp][i].fx[4+s0]
                                       *(0.5*(xx[jp][i].fx[1]+botp)+uu[jp][i].B0[1]);

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=6+s0;
                        vals[nv++]=cof/ne[2]*xx[jp][i].fx[4+s0]
                                      *(0.5*(xx[jp][i].fx[0]+xx[j][i].fx[0])+uu[jp][i].B0[0]);
                    }
                }

                else if (ir==1) {
                    if (i<N1) {
                        ne[1]=0.0; ve1[1]=0.0;ve2[1]=0.0;
                        for (l=0; l<sl; l++) {
                            s0=5*l;
                            ns=xx[j][ip].fx[4+s0];
                            ne[1]=ne[1]+ns;
                            ve1[1]=ve1[1]+ns*xx[j][ip].fx[5+s0];
                            ve2[1]=ve2[1]+ns*xx[j][ip].fx[6+s0];

                        }
                        ve1[1]=ve1[1]/ne[1];
                        ve2[1]=ve2[1]/ne[1];

                        if (i == 0) botn=-xx[j][i].fx[1];
                        else botn=xx[j][im].fx[1];

                        cof=0.5*dt2/h1h[i];

                        if (i>0) {
                            col[nv].j=j;  col[nv].i=im; col[nv].c=1;
                            vals[nv++]=-cof*ve1[0];
                        }
                        col[nv].j=jm; col[nv].i=i;   col[nv].c=0;
                        vals[nv++]=cof*ve2[0];

                        col[nv].j=j;  col[nv].i=i;   col[nv].c=0;
                        vals[nv++]=cof*ve2[0];

                        col[nv].j=j;  col[nv].i=i;   col[nv].c=1;
                        if (i==0) vals[nv++]=3.0+cof*ve1[1];
                        else vals[nv++]=3.0+cof*(ve1[1]-ve1[0]);

                        col[nv].j=jm; col[nv].i=ip;  col[nv].c=0;
                        vals[nv++]=-cof*ve2[1];

                        col[nv].j=j;  col[nv].i=ip;  col[nv].c=0;
                        vals[nv++]=-cof*ve2[1];

                        col[nv].j=j;  col[nv].i=ip;  col[nv].c=1;
                        vals[nv++]=cof*ve1[1];

                        cof=2.0*cof;

                        for (l=0; l<sl; l++) {
                            s0=5*l;
                            col[nv].j=j;  col[nv].i=i;  col[nv].c=4+s0;
                            vals[nv++]=cof/ne[0]
                                          *( (ve1[0]-xx[j][i].fx[5+s0])
                                             *(0.5*(xx[j][i].fx[1]+botn)+uu[j][i].B0[1])
                                            -(ve2[0]-xx[j][i].fx[6+s0])
                                             *(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]));

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=5+s0;
                            vals[nv++]=-cof/ne[0]*xx[j][i].fx[4+s0]
                                           *(0.5*(xx[j][i].fx[1]+botn)+uu[j][i].B0[1]);

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=6+s0;
                            vals[nv++]=cof/ne[0]*xx[j][i].fx[4+s0]
                                          *(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]);

                            col[nv].j=j;  col[nv].i=ip; col[nv].c=4+s0;
                            vals[nv++]=cof/ne[1]
                                          *( (ve2[1]-xx[j][ip].fx[6+s0])
                                             *(0.5*(xx[j][ip].fx[0]+xx[jm][ip].fx[0])+uu[j][ip].B0[0])
                                            -(ve1[1]-xx[j][ip].fx[5+s0])
                                             *(0.5*(xx[j][ip].fx[1]+xx[j][i].fx[1])+uu[j][ip].B0[1]));

                            col[nv].j=j;  col[nv].i=ip; col[nv].c=5+s0;
                            vals[nv++]=cof/ne[1]*xx[j][ip].fx[4+s0]
                                          *(0.5*(xx[j][ip].fx[1]+xx[j][i].fx[1])+uu[j][ip].B0[1]);

                            col[nv].j=j;  col[nv].i=ip; col[nv].c=6+s0;
                            vals[nv++]=-cof/ne[1]*xx[j][ip].fx[4+s0]
                                           *(0.5*(xx[j][ip].fx[0]+xx[jm][ip].fx[0])+uu[j][ip].B0[0]);
                        }
                    }
                    else {
                        ne[3]=0.0; ve1[3]=0.0; ve2[3]=0.0;
                        for (l=0; l<sl; l++) {
                            s0=5*l;
                            ns=xx[j][im].fx[4+s0];
                            ne[3]=ne[3]+ns;
                            ve1[3]=ve1[3]+ns*xx[j][im].fx[5+s0];
                            ve2[3]=ve2[3]+ns*xx[j][im].fx[6+s0];
                        }
                        ve1[3]=ve1[3]/ne[3];
                        ve2[3]=ve2[3]/ne[3];

                        cof=0.5*dt2/h1h[i];

                        col[nv].j=j;  col[nv].i=im2; col[nv].c=1;
                        vals[nv++]=-cof*ve1[3];

                        col[nv].j=j;  col[nv].i=im;  col[nv].c=1;
                        vals[nv++]=cof*(ve1[0]-ve1[3]);

                        col[nv].j=j;  col[nv].i=i;   col[nv].c=1;
                        vals[nv++]=3.0+cof*ve1[0];

                        col[nv].j=jm; col[nv].i=im;  col[nv].c=0;
                        vals[nv++]=cof*ve2[3];

                        col[nv].j=j;  col[nv].i=im;  col[nv].c=0;
                        vals[nv++]=cof*ve2[3];

                        col[nv].j=jm; col[nv].i=i;   col[nv].c=0;
                        vals[nv++]=-cof*ve2[0];

                        col[nv].j=j;  col[nv].i=i;   col[nv].c=0;
                        vals[nv++]=-cof*ve2[0];

                        cof=2.0*cof;

                        for (l=0; l<sl; l++) {
                            s0=5*l;
                            col[nv].j=j;  col[nv].i=im; col[nv].c=4+s0;
                            vals[nv++]=cof/ne[3]
                                          *( (ve1[3]-xx[j][im].fx[5+s0])
                                             *(0.5*(xx[j][im].fx[1]+xx[j][im2].fx[1])+uu[j][im].B0[1])
                                            -(ve2[3]-xx[j][im].fx[6+s0])
                                             *(0.5*(xx[j][im].fx[0]+xx[jm][im].fx[0])+uu[j][im].B0[0]));

                            col[nv].j=j;  col[nv].i=im; col[nv].c=5+s0;
                            vals[nv++]=-cof/ne[3]*xx[j][im].fx[4+s0]
                                           *(0.5*(xx[j][im].fx[1]+xx[j][im2].fx[1])+uu[j][im].B0[1]);

                            col[nv].j=j;  col[nv].i=im; col[nv].c=6+s0;
                            vals[nv++]=cof/ne[3]*xx[j][im].fx[4+s0]
                                          *(0.5*(xx[j][im].fx[0]+xx[jm][im].fx[0])+uu[j][im].B0[0]);

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=4+s0;
                            vals[nv++]=cof/ne[0]
                                          *( (ve2[0]-xx[j][i].fx[6+s0])
                                             *(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0])
                                            -(ve1[0]-xx[j][i].fx[5+s0])
                                             *(0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])+uu[j][i].B0[1]));

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=5+s0;
                            vals[nv++]=cof/ne[0]*xx[j][i].fx[4+s0]
                                          *(0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])+uu[j][i].B0[1]);

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=6+s0;
                            vals[nv++]=-cof/ne[0]*xx[j][i].fx[4+s0]
                                           *(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]);
                        }
                    }
                }

                else if (ir==2) {
                    if (i<N1) {
                        ne[1]=0.0; ne[2]=0.0; ne[5]=0.0;
                        ve1[1]=0.0; ve1[2]=0.0; ve1[5]=0.0;
                        ve2[1]=0.0; ve2[2]=0.0; ve2[5]=0.0;
                        ve3[1]=0.0; ve3[2]=0.0; ve3[5]=0.0;
                        for (l=0; l<sl; l++) {
                            s0=5*l;
                            ns=xx[j][ip].fx[4+s0];
                            ne[1]=ne[1]+ns;
                            ve1[1]=ve1[1]+ns*xx[j][ip].fx[5+s0];
                            ve2[1]=ve2[1]+ns*xx[j][ip].fx[6+s0];
                            ve3[1]=ve3[1]+ns*xx[j][ip].fx[7+s0];

                            ns=xx[jp][i].fx[4+s0];
                            ne[2]=ne[2]+ns;
                            ve1[2]=ve1[2]+ns*xx[jp][i].fx[5+s0];
                            ve2[2]=ve2[2]+ns*xx[jp][i].fx[6+s0];
                            ve3[2]=ve3[2]+ns*xx[jp][i].fx[7+s0];

                            ns=xx[jp][ip].fx[4+s0];
                            ne[5]=ne[5]+ns;
                            ve1[5]=ve1[5]+ns*xx[jp][ip].fx[5+s0];
                            ve2[5]=ve2[5]+ns*xx[jp][ip].fx[6+s0];
                            ve3[5]=ve3[5]+ns*xx[jp][ip].fx[7+s0];
                        }
                        ve1[1]=ve1[1]/ne[1];
                        ve1[2]=ve1[2]/ne[2];
                        ve1[5]=ve1[5]/ne[5];

                        ve2[1]=ve2[1]/ne[1];
                        ve2[2]=ve2[2]/ne[2];
                        ve2[5]=ve2[5]/ne[5];

                        ve3[1]=ve3[1]/ne[1];
                        ve3[2]=ve3[2]/ne[2];
                        ve3[5]=ve3[5]/ne[5];

                        vep[0]=ve1[5]+ve1[1];          //v_{e1,j+1/2,i+1}
                        vep[1]=ve1[2]+ve1[0];          //v_{e1,j+1/2,i}
                        vep[2]=ve3[5]+ve3[1];          //v_{e3,j+1/2,i+1}
                        vep[3]=ve3[2]+ve3[0];          //v_{e3,j+1/2,i}
                        B01p[0]=0.5*(uu[jp][ip].B0[0]+uu[j][ip].B0[0]);   //B_{01,j+1/2,i+1}
                        B01p[1]=0.5*(uu[jp][i].B0[0]+uu[j][i].B0[0]);     //B_{01,j+1/2,i}

                        if (i == 0) botn=-xx[j][i].fx[2];
                        else botn=xx[j][im].fx[2];

                        if (i>0) {
                            col[nv].j=j;  col[nv].i=im; col[nv].c=2;
                            vals[nv++]=-0.25*dt2/h1h[i]*vep[1];
                        }
                        col[nv].j=jm; col[nv].i=i;  col[nv].c=2;
                        vals[nv++]=-0.25*dt2/h2h[i]*(ve2[1]+ve2[0]);

                        col[nv].j=j;  col[nv].i=i;  col[nv].c=0;
                        vals[nv++]=dt2*(0.5/h1h[i]*vep[3]-0.125/rh[i]*(vep[2]+vep[3]));

                        col[nv].j=j;  col[nv].i=i;  col[nv].c=1;
                        vals[nv++]=0.5*dt2/h2h[i]*(ve3[1]+ve3[0]);

                        col[nv].j=j;  col[nv].i=i;  col[nv].c=2;
                        if (i==0) 
                          vals[nv++]=3.0+0.25*dt2*( vep[0]/h1h[i]
                                                   +(vep[0]+vep[1])/rh[i]
                                                   +(ve2[5]+ve2[2]-ve2[1]-ve2[0])/h2h[i]);
                        else
                          vals[nv++]=3.0+0.25*dt2*( (vep[0]-vep[1])/h1h[i]
                                                   +(vep[0]+vep[1])/rh[i]
                                                   +(ve2[5]+ve2[2]-ve2[1]-ve2[0])/h2h[i]);

                        col[nv].j=j;  col[nv].i=i;  col[nv].c=3;
                        vals[nv++]=dt*(ne[2]-ne[1])/(e*(ne[5]+ne[1]+ne[2]+ne[0])*hh[i]);

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=1;
                        vals[nv++]=-0.5*dt2/h2h[i]*(ve3[5]+ve3[2]);

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=2;
                        vals[nv++]=0.25*dt2/h2h[i]*(ve2[5]+ve2[2]);

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=3;
                        vals[nv++]=dt*(ne[5]-ne[0])/(e*(ne[5]+ne[1]+ne[2]+ne[0])*hh[i]);

                        col[nv].j=j;  col[nv].i=ip; col[nv].c=0;
                        vals[nv++]=-dt2*(0.5/h1h[i]*vep[2]-0.125/rh[i]*(vep[2]+vep[3]));

                        col[nv].j=j;  col[nv].i=ip; col[nv].c=2;
                        vals[nv++]=0.25*dt2/h1h[i]*vep[0];

                        col[nv].j=j;  col[nv].i=ip; col[nv].c=3;
                        vals[nv++]=dt*(ne[0]-ne[5])/(e*(ne[5]+ne[1]+ne[2]+ne[0])*hh[i]);

                        col[nv].j=jp; col[nv].i=ip; col[nv].c=3;
                        vals[nv++]=dt*(ne[1]-ne[2])/(e*(ne[5]+ne[1]+ne[2]+ne[0])*hh[i]);

                        cof=(ne[5]+ne[1]+ne[2]+ne[0]);

                        for (l=0; l<sl; l++) {
                            s0=5*l;

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=4+s0;
                            vals[nv++]= 0.25*dt2/ne[0]
                                            *( (ve1[0]-xx[j][i].fx[5+s0])
                                               *((xx[j][i].fx[2]+botn)/h1h[i]-xx[j][i].fx[2]/rh[i])
                                              +2.0*(ve3[0]-xx[j][i].fx[7+s0])
                                                  *( 0.25*(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1]))/rh[i]
                                                    -(xx[j][i].fx[0]+B01p[1])/h1h[i]
                                                    -(xx[j][i].fx[1]+0.5*(uu[j][ip].B0[1]+uu[j][i].B0[1]))/h2h[i])
                                              +(ve2[0]-xx[j][i].fx[6+s0])
                                               *(xx[j][i].fx[2]+xx[jm][i].fx[2])/h2h[i])
                                       +dt*(xx[j][ip].fx[3]-xx[jp][i].fx[3])/(e*cof*hh[i])
                                       -dt2*( ( xx[jp][ip].fx[3]+xx[jp][i].fx[3]
                                               -xx[j][ip].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[1]-ne[2]-ne[0])
                                             -( xx[jp][ip].fx[3]+xx[j][ip].fx[3]
                                               -xx[jp][i].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[2]-ne[1]-ne[0])
                                            )/(e*cof*cof*hh[i]);

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=5+s0;
                            vals[nv++]= 0.25*dt2*xx[j][i].fx[4+s0]/ne[0]
                                            *(xx[j][i].fx[2]/rh[i]-(xx[j][i].fx[2]+botn)/h1h[i]);

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=6+s0;
                            vals[nv++]=-0.25*dt2/(h2h[i]*ne[0])*xx[j][i].fx[4+s0]
                                            *(xx[j][i].fx[2]+xx[jm][i].fx[2]);

                            col[nv].j=j;  col[nv].i=i;  col[nv].c=7+s0;
                            vals[nv++]=0.5*dt2*xx[j][i].fx[4+s0]/ne[0]
                                          *( (xx[j][i].fx[0]+B01p[1])/h1h[i]
                                            -0.25*(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1]))/rh[i]
                                            +(xx[j][i].fx[1]+0.5*(uu[j][ip].B0[1]+uu[j][i].B0[1]))/h2h[i]);

                            col[nv].j=j;  col[nv].i=ip; col[nv].c=4+s0;
                            vals[nv++]= 0.25*dt2/ne[1]
                                        *( 2.0*(ve3[1]-xx[j][ip].fx[7+s0])
                                           *( (xx[j][ip].fx[0]+B01p[0])/h1h[i]
                                             +0.25*dt2
                                              *(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1]))/rh[i]
                                             -(xx[j][i].fx[1]+0.5*(uu[j][ip].B0[1]+uu[j][i].B0[1]))/h2h[i])
                                         -(ve1[1]-xx[j][ip].fx[5+s0])
                                          *((xx[j][ip].fx[2]+xx[j][i].fx[2])/h1h[i]+xx[j][i].fx[2]/rh[i])
                                         +(ve2[1]-xx[j][ip].fx[6+s0])
                                          *(xx[j][i].fx[2]+xx[jm][i].fx[2])/h2h[i])
                                       +dt*(xx[jp][ip].fx[3]-xx[j][i].fx[3])/(e*cof*hh[i])
                                       -dt2*( ( xx[jp][ip].fx[3]+xx[jp][i].fx[3]
                                               -xx[j][ip].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[1]-ne[2]-ne[0])
                                             -( xx[jp][ip].fx[3]+xx[j][ip].fx[3]
                                               -xx[jp][i].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[2]-ne[1]-ne[0])
                                            )/(e*cof*cof*hh[i]);

                            col[nv].j=j;  col[nv].i=ip; col[nv].c=5+s0;
                            vals[nv++]=0.25*dt2*xx[j][ip].fx[4+s0]/ne[1]
                                           *( (xx[j][ip].fx[2]+xx[j][i].fx[2])/h1h[i]
                                             +xx[j][i].fx[2]/rh[i]);

                            col[nv].j=j;  col[nv].i=ip; col[nv].c=6+s0;
                            vals[nv++]=-0.25*dt2/(h2h[i]*ne[1])*xx[j][ip].fx[4+s0]
                                            *(xx[j][i].fx[2]+xx[jm][i].fx[2]);

                            col[nv].j=j;  col[nv].i=ip; col[nv].c=7+s0;
                            vals[nv++]=0.5*dt2*xx[j][ip].fx[4+s0]/ne[1]
                                          *( (xx[j][i].fx[1]+0.5*(uu[j][ip].B0[1]+uu[j][i].B0[1]))/h2h[i]
                                            -(xx[j][ip].fx[0]+B01p[0])/h1h[i]
                                            -0.25*(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1]))/rh[i]);

                            col[nv].j=jp; col[nv].i=i;  col[nv].c=4+s0;
                            vals[nv++]= 0.25*dt2/ne[2]
                                        *( (ve1[2]-xx[jp][i].fx[5+s0])
                                           *((xx[j][i].fx[2]+botn)/h1h[i]-xx[j][i].fx[2]/rh[i])
                                          +2.0*(ve3[2]-xx[jp][i].fx[7+s0])
                                              *( 0.25*(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1])
                                                      )/rh[i]
                                                -(xx[j][i].fx[0]+B01p[1])/h1h[i]
                                                +(xx[jp][i].fx[1]+0.5*(uu[jp][ip].B0[1]+uu[jp][i].B0[1]))
                                                 /h2h[i])
                                          -(ve2[2]-xx[jp][i].fx[6+s0])
                                           *(xx[jp][i].fx[2]+xx[j][i].fx[2])/h2h[i])
                                       -dt*(xx[jp][ip].fx[3]-xx[j][i].fx[3])/(e*cof*hh[i])
                                       -dt2*( ( xx[jp][ip].fx[3]+xx[jp][i].fx[3]
                                               -xx[j][ip].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[1]-ne[2]-ne[0])
                                             -( xx[jp][ip].fx[3]+xx[j][ip].fx[3]
                                               -xx[jp][i].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[2]-ne[1]-ne[0])
                                            )/(e*cof*cof*hh[i]);

                            col[nv].j=jp; col[nv].i=i;  col[nv].c=5+s0;
                            vals[nv++]=0.25*dt2*xx[jp][i].fx[4+s0]/ne[2]
                                           *(xx[j][i].fx[2]/rh[i]-(xx[j][i].fx[2]+botn)/h1h[i]);

                            col[nv].j=jp; col[nv].i=i;  col[nv].c=6+s0;
                            vals[nv++]=0.25*dt2/(h2h[i]*ne[2])*xx[jp][i].fx[4+s0]
                                           *(xx[jp][i].fx[2]+xx[j][i].fx[2]);

                            col[nv].j=jp; col[nv].i=i;  col[nv].c=7+s0;
                            vals[nv++]=0.5*dt2*xx[jp][i].fx[4+s0]/ne[2]
                                          *( (xx[j][i].fx[0]+B01p[1])/h1h[i]
                                            -0.25*(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1]))/rh[i]
                                            -(xx[jp][i].fx[1]+0.5*(uu[jp][ip].B0[1]+uu[jp][i].B0[1]))/h2h[i]);

                            col[nv].j=jp; col[nv].i=ip; col[nv].c=4+s0;
                            vals[nv++]= 0.25*dt2/ne[5]
                                        *( 2.0*(ve3[5]-xx[jp][ip].fx[7+s0])
                                           *( (xx[j][ip].fx[0]+B01p[0])/h1h[i]
                                             +0.25*(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1]))
                                                   /rh[i]
                                             +(xx[jp][i].fx[1]+0.5*(uu[jp][ip].B0[1]+uu[jp][i].B0[1]))
                                              /h2h[i])
                                          -(ve1[5]-xx[jp][ip].fx[5+s0])
                                           *( (xx[j][ip].fx[2]+xx[j][i].fx[2])/h1h[i]
                                              -xx[j][i].fx[2]/rh[i])
                                          -(ve2[5]-xx[jp][ip].fx[6+s0])
                                           *(xx[jp][i].fx[2]+xx[j][i].fx[2])/h2h[i])
                                       +dt*(xx[jp][i].fx[3]-xx[j][ip].fx[3])/(e*cof*hh[i])
                                       -dt2*( ( xx[jp][ip].fx[3]+xx[jp][i].fx[3]
                                               -xx[j][ip].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[1]-ne[2]-ne[0])
                                             -( xx[jp][ip].fx[3]+xx[j][ip].fx[3]
                                               -xx[jp][i].fx[3]-xx[j][i].fx[3])
                                              *(ne[5]+ne[2]-ne[1]-ne[0])
                                            )/(e*cof*cof*hh[i]);

                            col[nv].j=jp; col[nv].i=ip; col[nv].c=5+s0;
                            vals[nv++]=0.25*dt2*xx[jp][ip].fx[4+s0]/ne[5]
                                           *((xx[j][ip].fx[2]+xx[j][i].fx[2])/h1h[i]+xx[j][i].fx[2]/rh[i]);                                            

                            col[nv].j=jp; col[nv].i=ip; col[nv].c=6+s0;
                            vals[nv++]=0.25*dt2/(h2h[i]*ne[5])*xx[jp][ip].fx[4+s0]
                                           *(xx[jp][i].fx[2]+xx[j][i].fx[2]);

                            col[nv].j=jp; col[nv].i=ip; col[nv].c=7+s0;
                            vals[nv++]=-0.5*dt2*xx[jp][ip].fx[4+s0]/ne[5]
                                           *( (xx[j][ip].fx[0]+B01p[0])/h1h[i]
                                             +0.25*(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1]))/rh[i]
                                             +(xx[jp][i].fx[1]+0.5*(uu[jp][ip].B0[1]+uu[jp][i].B0[1]))/h2h[i]);
                        }
                    }
                    else {
                        col[nv].j=j;  col[nv].i=im; col[nv].c=2;
                        vals[nv++]=-h1[im]/h1[i];
                        col[nv].j=j;  col[nv].i=i;  col[nv].c=2;
                        vals[nv++]=1.0;
                    }
                }
                else if (ir==3) {
                    if (i==0) {
                        col[nv].j=j;  col[nv].i=i;  col[nv].c=3;
                        vals[nv++]=1.0;
                    }
                    else if (i>0 && i<N1) {
                        cof=one6th*dt2/(ne[0]+nc*nc/ne[0]);

                        col[nv].j=j;  col[nv].i=im; col[nv].c=3;
                        vals[nv++]=cof*( (uu[j][ip].lamdae-uu[j][im].lamdae)/h12[i]
                                        +2.0*uu[j][i].lamdae*dh1[i]
                                        -4.0*uu[j][i].lamdae/h12[i]);

                        col[nv].j=jm; col[nv].i=i;  col[nv].c=3;
                        vals[nv++]=cof*( uu[jp][i].lamdae-uu[jm][i].lamdae
                                        -4.0*uu[j][i].lamdae)/h22[i];

                        col[nv].j=j;  col[nv].i=i;  col[nv].c=3;
                        vals[nv++]=3.0+8.0*cof*uu[j][i].lamdae*(1.0/h12[i]+1.0/h22[i]);

                        col[nv].j=jp; col[nv].i=i;  col[nv].c=3;
                        vals[nv++]=-cof*( uu[jp][i].lamdae-uu[jm][i].lamdae
                                         +4.0*uu[j][i].lamdae)/h22[i];

                        col[nv].j=j;  col[nv].i=ip; col[nv].c=3;
                        vals[nv++]=-cof*( (uu[j][ip].lamdae-uu[j][im].lamdae)/h12[i]
                                         +2.0*uu[j][i].lamdae*dh1[i]
                                         +4.0*uu[j][i].lamdae/h12[i]);

                        for (l=0; l<sl; l++) {
                          col[nv].j=j;  col[nv].i=i;   col[nv].c=4+5*l;
                          vals[nv++]=cof/(ne[0]+nc*nc/ne[0])*(1.0-nc*nc/(ne[0]*ne[0]))
                            *( ( (uu[j][ip].lamdae-uu[j][im].lamdae)/h12[i]
                                +2.0*uu[j][i].lamdae*dh1[i])
                               *(xx[j][ip].fx[3]-xx[j][im].fx[3])
                              +(uu[jp][i].lamdae-uu[jm][i].lamdae)/h22[i]
                               *(xx[jp][i].fx[3]-xx[jm][i].fx[3])
                              +4.0*uu[j][i].lamdae
                               *( (xx[j][ip].fx[3]-2.0*xx[j][i].fx[3]+xx[j][im].fx[3])/h12[i]
                                 +(xx[jp][i].fx[3]-2.0*xx[j][i].fx[3]+xx[jm][i].fx[3])/h22[i]));
                        }
                    }
                    else {
                        col[nv].j=j;  col[nv].i=im;  col[nv].c=3;
                        vals[nv++]=-1.0;
                        col[nv].j=j;  col[nv].i=i;   col[nv].c=3;
                        vals[nv++]=1.0;

                        /*cof=one6th*dt2/(ne[0]+nc*nc/ne[0]);

                        col[nv].j=j;  col[nv].i=im2; col[nv].c=3;
                        vals[nv++]=-cof*( ( 3.0*uu[j][i].lamdae-4.0*uu[j][im].lamdae
                                           +uu[j][im2].lamdae)/h12[i]
                                         +2.0*uu[j][i].lamdae*dh1[i]);

                        col[nv].j=j;  col[nv].i=im;  col[nv].c=3;
                        vals[nv++]=4.0*cof*( ( 3.0*uu[j][i].lamdae-4.0*uu[j][im].lamdae
                                              +uu[j][im2].lamdae)/h12[i]
                                            +2.0*uu[j][i].lamdae*dh1[i]);

                        col[nv].j=jm; col[nv].i=i;   col[nv].c=3;
                        vals[nv++]=cof*( uu[jp][i].lamdae-uu[jm][i].lamdae
                                        -4.0*uu[j][i].lamdae)/h22[i];

                        col[nv].j=j;  col[nv].i=i;   col[nv].c=3;
                        vals[nv++]=3.0+cof*( 8.0*uu[j][i].lamdae/h22[i]
                                            -3.0*( ( 3.0*uu[j][i].lamdae-4.0*uu[j][im].lamdae
                                                     +uu[j][im2].lamdae)/h12[i]
                                                  +2.0*uu[j][i].lamdae*dh1[i]));

                        col[nv].j=jp; col[nv].i=i;   col[nv].c=3;
                        vals[nv++]=-cof*( uu[jp][i].lamdae-uu[jm][i].lamdae
                                         +4.0*uu[j][i].lamdae)/h22[i];

                        for (l=0; l<sl; l++) {
                          col[nv].j=j;  col[nv].i=i;   col[nv].c=4+5*l;
                          vals[nv++]=cof/(ne[0]+nc*nc/ne[0])*(1.0-nc*nc/(ne[0]*ne[0]))
                            *( ( ( 3.0*uu[j][i].lamdae-4.0*uu[j][im].lamdae
                                  +uu[j][im2].lamdae)/h12[i]
                                +2.0*uu[j][i].lamdae*dh1[i])
                               *(3.0*xx[j][i].fx[3]-4.0*xx[j][im].fx[3]+xx[j][im2].fx[3])
                              +(uu[jp][i].lamdae-uu[jm][i].lamdae)/h22[i]
                               *(xx[jp][i].fx[3]-xx[jm][i].fx[3])
                              +4.0*uu[j][i].lamdae
                               *(xx[jp][i].fx[3]-2.0*xx[j][i].fx[3]+xx[jm][i].fx[3])/h22[i]);
                        }*/
                    }
                }
                else if ((ir-4) % 5 == 0) {
                    m=(ir-4)/5;
                    s0=5*m;
                    s4=4+s0;

                    if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
                        col[nv].j=j; col[nv].i=im; col[nv].c=s4;
                        vals[nv++]=-h1[im]/h1[i];
                        col[nv].j=j; col[nv].i=i;  col[nv].c=s4;
                        vals[nv++]=1.0;
                    }
                    else {
                        col[nv].j=j;  col[nv].i=i;  col[nv].c=s4;
                        vals[nv++]=1.0;
                    }
                }
                else if ((ir-5) % 5 == 0) {
                    m=(ir-5)/5;
                    s0=5*m;
                    s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0;

                    if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
                        /*for (l=0; l<sl; l++) {
                            s0=5*l;
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=4+s0;
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=5+s0;
                        }*/
                        col[nv].j=j; col[nv].i=im; col[nv++].c=s5;
                        col[nv].j=j; col[nv].i=i;  col[nv++].c=s5;
                    }
                    else {
                        if (i==0) {
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=s5;
                            col[nv].j=j; col[nv].i=ip; col[nv++].c=s5;
                        }
                        else if (i>0 && i<N1) {
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=1;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=1;
                            col[nv].j=jm; col[nv].i=im; col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=2;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=2;

                            for (l=0; l<sl; l++) {
                                s0=5*l;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=4+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=5+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=6+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=7+s0;
                            }
                        }
                        else {
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=s5;
                            col[nv].j=j; col[nv].i=im; col[nv++].c=s5;
                        }
                    }
                }
                else if ((ir-6) % 5 == 0) {
                    m=(ir-6)/5;
                    s0=5*m;
                    s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0;

                    if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
                        /*for (l=0; l<sl; l++) {
                            s0=5*l;
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=4+s0;
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=6+s0;
                        }*/
                        col[nv].j=j; col[nv].i=im; col[nv++].c=s6;
                        col[nv].j=j; col[nv].i=i;  col[nv++].c=s6;
                    }
                    else {
                        if (i==0) {
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=s6;
                        }
                        else if (i==N1) {
                            if (nvt==0 && btop[j]>-1.0e5) {
                                col[nv].j=j; col[nv].i=i;  col[nv++].c=s6;
                            }
                            else {
                                col[nv].j=j; col[nv].i=im; col[nv++].c=s6;
                                col[nv].j=j; col[nv].i=i;  col[nv++].c=s6;
                            }
                        }
                        else {
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=1;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=1;
                            col[nv].j=jm; col[nv].i=im; col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=2;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=2;

                            for (l=0; l<sl; l++) {
                                s0=5*l;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=4+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=5+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=6+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=7+s0;
                            }
                        }
                    }
                }
                else if ((ir-7) % 5 == 0) {
                    m=(ir-7)/5;
                    s0=5*m;
                    s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0;

                    if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
                        /*for (l=0; l<sl; l++) {
                            s0=5*l;
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=4+s0;
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=7+s0;
                        }*/
                        col[nv].j=j; col[nv].i=im; col[nv++].c=s7;
                        col[nv].j=j; col[nv].i=i;  col[nv++].c=s7;
                    }
                    else {
                        if (i==0) {
                            col[nv].j=j; col[nv].i=i;  col[nv++].c=s7;
                        }
                        else if (i==N1) {
                            if (nvt==1 && btop[j]>-1.0e5) {
                                col[nv].j=j; col[nv].i=i;  col[nv++].c=s7;
                            }
                            else {
                                col[nv].j=j; col[nv].i=im; col[nv++].c=s7;
                                col[nv].j=j; col[nv].i=i;  col[nv++].c=s7;
                            }
                        }
                        else {
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=1;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=1;
                            col[nv].j=jm; col[nv].i=im; col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=2;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=2;

                            for (l=0; l<sl; l++) {
                                s0=5*l;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=4+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=5+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=6+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=7+s0;
                            }
                        }
                    }
                }
                else if ((ir-8) % 5 == 0) {
                    m=(ir-8)/5;
                    s0=5*m;
                    s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0; s8=8+s0;

                    if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
                        col[nv].j=j; col[nv].i=im; col[nv++].c=s8;
                        col[nv].j=j; col[nv].i=i;  col[nv++].c=s8;
                    }
                    else {
                        if (i == 0) {
                            col[nv].j=j; col[nv].i=i; col[nv++].c=s8;
                        }
                        else if (i>0 && i<N1) {
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=0;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=1;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=1;
                            col[nv].j=jm; col[nv].i=im; col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=2;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=2;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=3;

                            col[nv].j=j;  col[nv].i=im; col[nv++].c=s8;
                            col[nv].j=jm; col[nv].i=i;  col[nv++].c=s8;
                            col[nv].j=jp; col[nv].i=i;  col[nv++].c=s8;
                            //if (i < N1) {
                            col[nv].j=j;  col[nv].i=ip; col[nv++].c=s8;
                            //}
                            //else {
                            //    col[nv].j=j;  col[nv].i=im2; col[nv++].c=s8;
                            //}

                            for (l=0; l<sl; l++) {
                                s0=5*l;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=4+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=5+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=6+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=7+s0;
                                col[nv].j=j;  col[nv].i=i;  col[nv++].c=8+s0;
                            }
                        }
                        else {
                            col[nv].j=j;  col[nv].i=im; col[nv++].c=s8;
                            col[nv].j=j;  col[nv].i=i;  col[nv++].c=s8;
                        }
                    }
                }

                if (ir>3 && ((ir-4) % 5) != 0) {
                    for (m=0; m<nv; m++) {
                        jj=col[m].j; ii=col[m].i; cc=col[m].c;

                        x1=xx[jj][ii].fx[cc];
                        if (x1 != 0.0) {
                            absx=fabs(x1);
                            if (absx > 1.0) segm=x1*dx;
                            else segm=x1/absx*dx;
                        }
                        else segm=dx;

                        xx[jj][ii].fx[cc]=x1-segm;
                        ffm=function_part(j,i,ir,yj,xi,xx,xn,xn1,uu);

                        xx[jj][ii].fx[cc]=x1+segm;
                        ffp=function_part(j,i,ir,yj,xi,xx,xn,xn1,uu);
                        xx[jj][ii].fx[cc]=x1;

                        vals[m]=(ffp-ffm)/(2.0*segm);
                    }
                }

                ierr=MatSetValuesStencil(Jpre,1,&row,nv,col,vals,INSERT_VALUES);CHKERRQ(ierr);
            }
        }
    }

    ierr = MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    if (Jac != Jpre) {
        ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    delete[] vals;
    delete[] col;

    /* Restore vectors */
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
