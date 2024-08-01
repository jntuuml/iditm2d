/* initialize
 *   Call house keeping routines to input parameters, set up initial conditions,
 *   specifically, input state variables from files and so on.
 *
 * Jiannan Tu 9/13/2013, 12/6/2013
 */
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <petscdmda.h>

#include <cmath>
#include "param.h"
#include "funcdef.h"

using namespace std;

int initialize(DM da,Vec X,AppCtx user,char *workdir)
{
    int     xs,xm,ys,ym,ierr,j,i,l,s0;
    Field   **xx,**xn,**xn1;
    Fieldu  **uu;

    ierr=DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

    /* normalize constants  */
    const_normalize();

    /* allocate memory for arrays */
    array_alloc(xs,xm,ys,ym);

    /* calculate parameters of curvilinear coordinates. also magnetic
     * polar coordinates (r, theta) of the grids */
    curvicoords();

    /* magnetic longitude of the dayside magnetic meridian */
    phi[0]=mag_longitude();

    /* geographic coordinates of the grids in geomagnetic meridian */
    geocoords();

/* ---- solar EUV flux for given f107 & f107a */
    if (euvflux_seg(da,workdir) < 0) return -1;

    /* set distributions of ion density and temperature of given number
     * of species, and distribution of electron temperatures */
    ierr = DMDAVecGetArray(da,X,&xx);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user.Xn,&xn);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,user.Xn1,&xn1);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user.dau,user.U,&uu);CHKERRQ(ierr);

    if (input_prev(da,xx,uu,workdir,user.prefln) < 0) return -1;
    if (smod==0) {
        for (j=ys; j<ys+ym; j++) {
            for (i=xs; i<xs+xm; i++) {
                for (l=0; l<3; l++) {
                    xn[j][i].fx[l]=xx[j][i].fx[l];
                    xn1[j][i].fx[l]=xx[j][i].fx[l];
                }
                xn[j][i].fx[3]=xx[j][i].fx[3];
                xn1[j][i].fx[3]=xx[j][i].fx[3];
                for (l=0; l<sl; l++) {
                    s0=5*l;
                    xn[j][i].fx[4+s0]=xx[j][i].fx[4+s0];
                    xn[j][i].fx[5+s0]=xx[j][i].fx[5+s0];
                    xn[j][i].fx[6+s0]=xx[j][i].fx[6+s0];
                    xn[j][i].fx[7+s0]=xx[j][i].fx[7+s0];
                    xn[j][i].fx[8+s0]=xx[j][i].fx[8+s0];

                    xn1[j][i].fx[4+s0]=xx[j][i].fx[4+s0];
                    xn1[j][i].fx[5+s0]=xx[j][i].fx[5+s0];
                    xn1[j][i].fx[6+s0]=xx[j][i].fx[6+s0];
                    xn1[j][i].fx[7+s0]=xx[j][i].fx[7+s0];
                    xn1[j][i].fx[8+s0]=xx[j][i].fx[8+s0];
                }
            }
        }
    }
    else {
        if (input_prev(da,xn,uu,workdir,user.preflm) < 0) return -1;
        if (input_prev(da,xn1,uu,workdir,user.preflo) < 0) return -1;
    }

    ierr = DMDAVecRestoreArray(da,X,&xx);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user.Xn,&xn);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,user.Xn1,&xn1);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user.dau,user.U,&uu);CHKERRQ(ierr);
    
    return 0;
}