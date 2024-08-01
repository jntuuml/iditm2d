#include "param.h"
#include "funcdef.h"

#undef __FUNCT__
#define __FUNCT__ "jparalle"

void jparallel(Field **xx,Fieldu **uu,int xs, int xm, int ys, int ym)
{
    PetscInt       i,j,jm,im;
    double         BB, J1, J2;

    /* Compute function over the locally owned part of the grid */
    for (j=ys; j<ys+ym; j++) {
        jm=j-1;

        for (i=xs; i<xs+xm; i++) {
            //B^2
            BB=uu[j][i].B0[0]*uu[j][i].B0[0]+uu[j][i].B0[1]*uu[j][i].B0[1];

            if(i>0) {
                im=i-1;
                J1= 0.5*(xx[j][i].fx[2]+xx[j][im].fx[2]-xx[jm][i].fx[2]-xx[jm][im].fx[2])/h2[i];
                J2=-0.5*(xx[j][i].fx[2]+xx[jm][i].fx[2]-xx[j][im].fx[2]-xx[jm][im].fx[2])/h1[i];

                //current density strength
                J1=(J1*uu[j][i].B0[0]+J2*uu[j][i].B0[1])/BB;

                //parallel current components in x1 & x2 directions
                uu[j][i].Jpar[0]=uu[j][i].B0[0]*J1;
                uu[j][i].Jpar[1]=uu[j][i].B0[1]*J1;
            }
        }
    }
}
