/*************************************************************************
 * neu_cooling_rate.cpp  
 *    calculate neutral cooling rates due to infrared radiative loss by
 * O at wavelength of 63 micro-meter and 147 micro-meter, NO at 5.3 micro-meter  
   Input:  parms  - system parameters
              xx  - solution vector

    Output: none
 *  
 * Jiannan Tu
 * 1/24/2014
*************************************************************************/
/*#include <iostream>
#include <fstream>
#include <iomanip>*/

#include <cmath>
#include "param.h"
#include "funcdef.h"

using namespace std;

void neu_cooling_rate(Fieldu **uu,int xs,int xm,int ys,int ym)
{
    int    i, j, xup;
    double n00, t0p0, TO, nO, nNO, Tx1, Tx2, de, chi, E21, E22, TOb[1], dr;

    t0p0=t0/p0;
    n00=n0*1.0e-6;

    xup=xs+xm;

    double *tao = new double[xup-xs];

    for (j = ys; j < ys+ym; j++) {

        for (i = xs; i < xup; i++) tao[i-xs]=0.0;

        TOb[0]=200.0;

        /* first loop i is to calculate optical depth on all grids */
        for (i = xup-2; i >= xs; i--) {
            /* density in cm^-3 & temperature in K (not using arrays for nO, nNO, TO 
             * because of better efficiency when using openMP */
            nO=exp(uu[j][i].fu[0])*n00;
            nNO=exp(uu[j][i].fu[20])*n00;
            TO=uu[j][i].fu[4]*T0;

            Tx1=exp(228.0/TO);

            //calculate optical depth
            if (i > 0) dr=(rr[i]-rr[i-1])*r0;
            else dr=(rr[i+1]-rr[i])*r0;
            tao[i-xs]=tao[i-xs+1]+1.0e-14/sqrt(TO)*nO*(Tx1-1.0)
                                         /(0.6+0.2*exp(-98.0/TO)+Tx1);
        }

        for (i = xs; i < xup; i++) {
            nO=exp(uu[j][i].fu[0])*n00;
            nNO=exp(uu[j][i].fu[20])*n00;
            TO=uu[j][i].fu[4]*T0;
            if (i == 0) TOb[0]=TO;

            Tx1=exp(228.0/TO);
            Tx2=exp(-326.0/TO);

            //calculate chi for reduction factor
            E21=expon_integral(tao[i-xs], 60);
            E22=expon_integral(tao[0]-tao[i-xs],60);
            chi=0.5*(Tx1-1.0)*((2.0-E21-E22)/(Tx1-1.0)+E22/(exp(228.0/TOb[0])-1.0));

            /* normalized cooling rate (in Joule m^-3 s^-1 before normalization) */
            de=(1.0+0.6/Tx1+0.2*Tx2);
            Cn[j-ys][i-xs]=nO*( (1.0-chi)/de*(1.69e-19/Tx1+4.59e-21*Tx2)
                         +3.24015e-23*nNO*exp(-2714.57/TO)/(6.5e-11*nO+13.3))*t0p0;
        }
    }

    delete[] tao;

    return;
}

