/* array_dealloc.cpp
 *  dynamically allocate memory for arrays
 * 
 *  Jiannan Tu
 *  9/15/2013
 */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"

void array_dealloc(DM da)
{
    int i, j, xs, xm, ys, ym, yj, xi;
    MPI_Comm comm;

    PetscObjectGetComm((PetscObject)da,&comm);

    /*
      Get local grid boundaries
    */
    DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);

    delete[] btop;

    delete[] rr;
    delete[] rh;
    delete[] theta;
    delete[] phi;
    delete[] zh;

    delete[] h1;
    delete[] h2;
    delete[] h12;
    delete[] h22;
    delete[] h1h;
    delete[] h2h;
    delete[] hh;
    delete[] dh1;

    delete[] gr;

    delete[] Omega1;
    delete[] Omega2;
    delete[] Omega1h;
    delete[] Omega2h;

    delete[] thetag;
    delete[] phig;

    delete[] zenith;

    for(j=ys; j<ys+ym; j++) {
        yj=j-ys;

        for(i=xs; i<xs+xm; i++) {
            xi=i-xs;

            delete[] qi[yj][xi];
            delete[] qib[yj][xi];

            delete[] Ps[yj][xi];
            delete[] Ls[yj][xi];

            delete[] Qeuv[yj][xi];
            delete[] Qch[yj][xi];
        }
 
        delete[] qi[yj];
        delete[] qib[yj];
        delete[] qit[yj];

        delete[] Qeuv[yj];
        delete[] Qch[yj];
        delete[] Cn[yj];

        delete[] Ps[yj];
        delete[] Ls[yj];

        delete[] Qe[yj];
        delete[] Ce[yj];
    }

    delete[] qi;
    delete[] qib;
    delete[] qit;

    delete[] Ps;
    delete[] Ls;
    
    delete[] Qeuv;
    delete[] Qch;
    delete[] Cn;

    delete[] Qe;
    delete[] Ce;

    return;
}
