/*************************************************************************
 * ele_cooling_rate.cpp  
 * calculate electron cooling rate by excitation of neutral states through
 * inelastic electron - neutral collisions. 
 * See Schun and Nagy, Geophy. Review, 113(8), A08307, 2008
 * 
   Input:  parms  - system parameters
              xx  - solution vector

    Output: non
 *  
 * Jiannan Tu
 * 1/22/2014
*************************************************************************/
#include <cmath>
#include <ctime>

using namespace std;

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"

void ele_cooling_rate(Field **xx,Fieldu **uu,int xs,int xm,int ys,int ym)
{
    int    i, j, m;
    double Te, ne, Td, fc, gc, he, Z, Dxi[3], Exi[3], dd;
    double nO, nO2, nN2, nH, TO, TO2, TN2, TH, Le[10], n00, cci;

    const double ei[3]={0.02, 0.028, 0.008};
    const double Ei[3]={228.0, 326.0, 98.0};
    const double Ai[3]={8.58e-6, 7.201e-6, 2.463e-7};
    const double Bi[3]={1.008, 0.9617, 1.1448};
    const double Ci[3]={1.009, 0.9444, 1.466};

    n00=n0*1.0e-6;
    cci=q*t0/p0*1.0e6;

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            // density in cm^-3
            nO =exp(uu[j][i].fu[0])*n00;
            nO2=exp(uu[j][i].fu[5])*n00;
            nN2=exp(uu[j][i].fu[10])*n00;
            nH =exp(uu[j][i].fu[15])*n00;

            ne=0.0;
            for (m = 0; m < sl; m++) ne=ne+xx[j][i].fx[4+5*m];
            ne=ne*n00;

            Te=xx[j][i].fx[3]*T0;
            TO =uu[j][i].fu[4]*T0;
            TO2=uu[j][i].fu[9]*T0;
            TN2=uu[j][i].fu[14]*T0;
            TH =uu[j][i].fu[19]*T0;

            Td=sqrt(Te);

            if (Te > TN2) {
                /* rate (divided by ne) of cooling due to impact excitation of
                 * N2 rotation in eV s^-1 (not multiplied by ne yet) */
                Le[0]=2.9e-14*nN2*(Te-TN2)/Td;

                /* due to N2 vibration */
                fc=1.06e4+7.51e3*tanh(1.10e-3*(Te-1800.0));
                gc=3300.0+(1.233-2.056e-4*(Te-4000.0))*(Te-1000.0);

                //original expression regards cooling as negative heating
                Le[1]=2.99e-12*nN2*exp(fc*((Te-2000.0)/(2000.0*Te)))
                                  *(1.0-exp(gc*(TN2-Te)/(Te*TN2)));

                //cooling due to elastic collisions with N2
                Le[2]=1.77e-19*nN2*(1.0-1.21e-4*Te)*Te*(Te-TN2);
            }
            else {
                Le[0]=0.0; Le[1]=0.0; Le[2]=0.0;
            }

            if (Te > TO2) {
                /* due to O2 rotation */
                Le[3]=6.9e-14*nO2*(Te-TO2)/Td;

                /* due to O2 vibration */
                he=3300.0-839.0*sin(1.91e-4*(Te-2700.0));

                //original expression regards cooling as negative heating
                Le[4]=5.196e-13*nO2*exp(he*(Te-700.0)/(700.0*Te))
                                   *(1.0-exp(2770.0*(TO2-Te)/(Te*TO2)));

                //cooling due to elastic collisions with O2
                Le[5]=1.21e-18*nO2*(1.0+3.6e-2*Td)*Td*(Te-TO2);
            }
            else {
                Le[3]=0.0; Le[4]=0.0; Le[5]=0.0;
            }

            if (Te > TO) {
                /* due to O fine structure */
                Z=5.0+3.0*exp(-228.0/TO)+exp(-326.0/TO);
                Dxi[0]=exp(-228.0/TO);
                Dxi[1]=exp(-326.0/TO);
                Dxi[2]=exp(-326.0/TO);

                Exi[0]=exp(-228.0/Te);
                Exi[1]=exp(-326.0/Te);
                Exi[2]=exp(-(98.0/Te+228.0/TO));

                Le[6]=0.0;
                for (m = 0; m < 3; m++) {
                //original expression regards cooling as negative heating
                    Le[6]=Le[6]+Ai[m]*Ci[m]*pow(Te, (Bi[m]-0.5))
                                     *( ei[m]*(Exi[m]-Dxi[m])
                                       +5.91e-9*(Te-TO)*( (1.0+Bi[m])*Dxi[m]
                                                         +(Ei[m]/Te+1.0+Bi[m])*Exi[m]));
                }
                Le[6]=Le[6]*8.629e-6*nO/Z;

                /* due to O(1D) excitation */
                dd=2.4e4+(0.3-1.947e-5*(Te-4000.0))*(Te-1500.0);

                //original expression regards cooling as negative heating
                Le[7]=1.57e-12*nO*exp(dd*(Te-3000.0)/(3000.0*Te))
                                 *(1.0-exp(22713.0*(TO-Te)/(Te*TO)));

                //cooling due to elastic collisions with O
                Le[8]=7.9e-19*nO*(1.0+5.7e-4*Te)*Td*(Te-TO);
            }
            else {
                Le[6]=0.0; Le[7]=0.0; Le[8]=0.0;
            }

            //cooling due to elastic collisions with H
            if (Te > TH) Le[9]=9.63e-16*nH*(1.0-1.35e-4*Te)*Td*(Te-TH);
            else Le[9]=0.0;

            //electron Coulomb collision cooling rate has been included in the
            //electron temperature equation
            /*Ti=0.0; mst=0.0;
            for (m = 0; m < sl; m++) {
                Ti=Ti+ms[m]*xx[j][i].fx[11+5*m]; mst=mst+ms[m];
            }
            Ti=Ti/mst;

            if (Te > Ti) Le[10]=3.2e-8*(Te-Ti)/(Td*Td*Td)*logA
                                      *( xx[j][i].fx[7]+16.0*xx[j][i].fx[12]+0.5*xx[j][i].fx[17]
                                        +0.53*xx[j][i].fx[22])*n00;
            else Le[10]=0.0;*/

            //normalized cooling rate
            Ce[j-ys][i-xs]=(Le[0]+Le[1]+Le[2]+Le[3]+Le[4]+Le[5]+Le[6]+Le[7]+Le[8]+Le[9])*ne*cci;
            if (zh[i]<=110.0e3 && Ce[j-ys][i-xs]>1.0e-7) Ce[j-ys][i-xs]=1.0e-7; 
        }
    }

    return;
}
