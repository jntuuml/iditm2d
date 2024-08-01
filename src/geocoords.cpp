/******************************************************************************
 *  geocoords
 *   Calculate geographic colatitude and longitude of grid points in magnetic
 *   meridian
 *
 *   Jiannan Tu
 *   12/5/2013, 1/10/2014
 ******************************************************************************/

//for test purpose
/*#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>

using namespace std;*/

#include "funcdef.h"
#include "param.h"

void geocoords()
{
    double  xmag, ymag, zmag, xgeo, ygeo, zgeo;
    double  rm, rg, colat, phin;
    int     j;

    rm=rr[0];
    for (j = 0; j <= N2; j++) {
        if (theta[j] < pi) {
            colat=theta[j];
            phin=phi[0];
        }
        else {
            colat=pi2-theta[j];
            phin=phi[0]+pi;
            if (phin < 0.0) phin=pi2+phin;
            else if (phin > pi2) phin=phin-pi2;
        }
        sphcar(rm, colat, phin, xmag, ymag, zmag, 1);
        geomag(xgeo, ygeo, zgeo, xmag, ymag, zmag, -1);
        sphcar(rg, colat, phin, xgeo, ygeo, zgeo, -1);

        thetag[j]=colat;
        phig[j]=phin;
    }

    return;
}
