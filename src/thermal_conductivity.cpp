/*************************************************************************
 * thermal_conductivity.cpp  
 *    calculate thermal conductivities for electrons, ion species & neutral
 *  species
   Input:  parms  - system parameters
              xx  - solution vector

    Output: non
 *  
 * Jiannan Tu
 * 1/22/2014
*************************************************************************/
#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"
#include <cmath>
#include "param.h"

using namespace std;

void thermal_conductivity(Field **xx,Fieldu **uu,int xs,int xm,int ys,int ym)
{
    int    i, j, m, l, ll;
    double Te, ne, Td, qn[4], nqd, nuei, nuec, Dst, Tnn, Ti;
    double ni[7], nuss[7], nuen, n0t0, fq;
    const double z800=800.e3, nc=1.0e8/n0;

    n0t0=1.0e-6*n0*t0;

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            ne=0.0; nuei=0.0;
            for (l = 0; l< sl; l++) {
                ni[l]=xx[j][i].fx[4+5*l];
                ne=ne+ni[l];
                nuei=nuei+uu[j][i].nu[l];
            }

            nuen=0.0;
            for (l=sl; l<slm; l++) nuen=nuen+uu[j][i].nu[l];

            Te=xx[j][i].fx[3]*T0;
            Td=sqrt(Te);

            qn[0]=1.1e-16*(1.0+5.7e-4*Te);
            qn[1]=2.2e-16*(1.0+0.036*Td);
            qn[2]=2.82e-17*Td*(1.0-1.21e-4*Te);
            qn[3]=5.47e-15*(1.0-1.35e-4*Te);

            nqd= exp(uu[j][i].fu[0])*qn[0]+exp(uu[j][i].fu[5])*qn[1]
                +exp(uu[j][i].fu[10])*qn[2]+exp(uu[j][i].fu[15])*qn[3];

            //normalized electron thermal conductivity
            fq=1.0/(1.0+nc*nc/(ne*ne));
            uu[j][i].lamdae=1.233694e-11*fq*pow(Te, 2.5)/((1.0+3.32e4*Te*Te/ne*nqd)*lamda0);            

            //normalized thermal-electric conductivity (Shunk & Nagy, 2000, p. 134-135)
            nuec=3.85373195747e1*ne*n0t0/pow(Te, 1.5)+1.625*nuei+3.125*nuen;
            uu[j][i].betae=1.6157861066e-4*fq*(nuei/nuec)*Te/beta0;

            /* normalized ion thermal conductivity */
            for (m = 0; m < sl; m++) {
                nqd=0.0;
                for (l = 1; l < slm; l++) {
                    if (l<=m) ll=l-1; else ll=l;
                    if (l < sl) Dst=(3.0*ms[m]*ms[m]-0.2*ms[ll]*ms[ll]+0.1*ms[m]*ms[ll])
                                    /((ms[m]+ms[ll])*(ms[m]+ms[ll]));
                    else Dst=(3.0*ms[m]*ms[m]+ms[ll]*ms[ll]+1.6*ms[m]*ms[ll])
                             /((ms[m]+ms[ll])*(ms[m]+ms[ll]));

                    nqd=nqd+uu[j][i].nu[9+9*m+l]*(Dst+1.5*ms[ll]/(ms[m]+ms[ll]));
                }

                fq=1.0/(1.0+pow(nc/xx[j][i].fx[4+5*m],3));

                Ti=xx[j][i].fx[8+5*m]*T0;
                nuss[m]=1.27*ni[m]*n0t0/(sqrt(ms[m])*pow(Ti, 1.5));
                nqd=1.0+1.25*nqd/nuss[m];

                uu[j][i].lamdas[m]=4.96682e-13*fq*pow(Ti, 2.5)/sqrt(ms[m])/(nqd*lamda0);
            }

            /* normalized neutral thermal conductivity*/
            fq=1.0/(1.0+pow(zh[i]/z800,4));
            uu[j][i].lamdas[5] =7.59e-4*fq*pow(uu[j][i].fu[4]*T0, 0.69)/lamda0;            // O
            Tnn=uu[j][i].fu[9]*T0;
            uu[j][i].lamdas[6]=1.00e-4*fq*(3.93*pow(Tnn, 0.69)+0.255*Tnn-9.27)/lamda0;  // O2
            Tnn=uu[j][i].fu[14]*T0;
            uu[j][i].lamdas[7]=1.00e-4*fq*(3.82*pow(Tnn, 0.69)+0.190*Tnn+5.14)/lamda0;  // N2
            uu[j][i].lamdas[8]=3.79e-5*fq*pow(uu[j][i].fu[19]*T0, 0.69)/lamda0;            // H
        }
    }

    return;
}
