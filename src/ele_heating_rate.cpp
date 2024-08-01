/*************************************************************************
 * ele_heating_rate.cpp  
 * calculate electron heating rate by photoelectrons, using parameterization
 * developed by Smithtro and Solomon, JGR, 113(8), A08307, 2008
 * 
   Input:  parms  - system parameters
              xx  - solution vector

    Output: non
 *  
 * Jiannan Tu
 * 1/22/2014
*************************************************************************/
#include <cmath>

using namespace std;

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"

void ele_heating_rate(Field **xx,Fieldu **uu,int xs,int xm,int ys,int ym)
{
    int    i, j, l, m, i0, yj, xi;
    double R, epsn, aveps1, aveps2, Q1, Q2;
    double y, BB, BB0, ne, ne0, nel;

    i0=0;

    for (i = 0; i <= N1; i++) if(zh[i] <= 300.0e3) i0=i;

    //average photon energy in 0-55 nm
    epsn=0.0;
    for (l = 0; l < 16; l++) epsn=epsn+epi[l];
    epsn=epsn/16.0;

    for (j=ys; j<ys+ym; j++) {
        yj=j-ys;

        BB0=sqrt(uu[j][i0].B0[0]*uu[j][i0].B0[0]+uu[j][i0].B0[1]*uu[j][i0].B0[1]);

        ne0=0.0;
        for (l = 0; l < sl; l++) ne0=ne0+xx[j][i0].fx[4+5*l];

        for (i=xs; i<xs+xm; i++) {
            xi=i-xs;

            if (i <= i0) {
                ne=0.0;
                for (l = 0; l < sl; l++) ne=ne+xx[j][i].fx[4+5*l];

                R=ne/(exp(uu[j][i].fu[0])+exp(uu[j][i].fu[5])+exp(uu[j][i].fu[10]));

                aveps1=0.0;
                aveps2=0.0;
                for (l = 0; l < 7; l++) {
                    aveps1=aveps1+cc1[l]*pow(log(R), double(l));
                    aveps2=aveps2+cc2[l]*pow(log(R), double(l));
                }
                if (aveps1 < -45.0) aveps1=-45.0;
                if (aveps2 < -45.0) aveps2=-45.0;
                aveps1=exp(aveps1);
                aveps2=exp(aveps2);

                //heating rate for wavelength 0-55 nm & 55-105 nm multiplied by t0
                Q1=0.0;
                for (l = 0; l < 16; l++) Q1=Q1+epi[l]*qib[yj][xi][l];
                if (Q1 < 1.0e-20) Q1=1.0e-20;

                Q1=aveps1*Q1;
                Q2=aveps2*qit[yj][xi];
                if (Q2 < 1.0e-20) Q2=1.0e-20;

                //normalized photoelectron heating rate in the local heating region
                Qe[yj][xi]=q*(Q1+Q2)*n0/p0;
            }
            else {
                y=0.0; nel=ne0;
                for (l = i0; l < i; l++) {
                    ne=0.0;
                    for (m = 0; m < sl; m++) ne=ne+xx[j][l+1].fx[4+5*m];
                    y=y+0.5*(nel+ne)*(rr[l+1]-rr[l]);
                    nel=ne;
                }
                y=y*r0*n0;

                //magnetic field strength
                BB=sqrt(uu[j][i].B0[0]*uu[j][i].B0[0]+uu[j][i].B0[1]*uu[j][i].B0[1]);

                ne=0.0;
                for (l = 0; l < sl; l++) ne=ne+xx[j][i].fx[4+5*l];

                //normalized photoelectron heating rate in transport-dominated region
                Qe[yj][xi]=ne/ne0*BB/BB0*Qe[yj][i0]*exp(-7.0e-18*y);
            }
        }
    }

    return;
}
