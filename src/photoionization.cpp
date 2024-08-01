/*************************************************************************
 * photoionization.cpp  
 *    function to calculate neutral photo-ionization coefficients. The 
 * normalization is done by using normalized neutral densities and 
 * multiplying qi by t0. Also compute EUV heating rates for O, O2, & N2 
 *  
   Input:  parms  - system parameters
              xx  - solution vector

    Output: non
 *  
 * Jiannan Tu
 * 9/26/2013, 12/14/2013, 1/22/2014
*************************************************************************/
#include <cmath>
using namespace std;

#include "funcdef.h"
#include "param.h"

void photoionization(Fieldu **uu,int xs,int xm,int ys,int ym)
{
    int    i, j, m, l, i280;
    double pih, r0n0, GHy120, GH1120;
    double tsin, tcos, r, grr, mass, Tn0, x, y, dc, euvf, Hn[5], Ch[5], tao[37];
    double Nnn[5], A0[4], GHy, GH1, GH2, Fph[4], f107p, Ih[4], Adm;
    const double dem[5]={0.1,0.3,0.5,0.7,0.9};
    int   yj, xi;

    pih=double(0.5)*pi;

    i280=0;

    r0n0=r0*n0;

    //values of Ghy & Gh1 at 120 deg (GH2 at 120 deg =1)
    GHy120=3.62-1.7/pow(0.5, 0.37);
    GH1120=exp(-4.05);

    //altitude index at 280 km
    for (i=xs; i<xs+xm; i++) {
        if (zh[i] >= 280.0e3) {i280=i; break;}
    }

    solar_zenith();

    for (j=ys; j<ys+ym; j++) {
        yj=j-ys;

        tsin=sin(zenith[j]);
        tcos=abs(cos(zenith[j]));

        //calculate daytime photoionization rates up to height at z[ipc]
        for (i=xs; i<xs+xm; i++) {
          xi=i-xs;

          if (i <= ipc) {
            r=rr[i]*r0;
            grr=ge*pow(Re/r, 2.0);

            //optical depth due to absorption by O, O2, N2, & N
            for (m = 0; m < 4; m++) {
                /* first calculate scale height of neutral species */
                if (m < 3) {
                    mass=ms[sl+m];
                    Tn0=uu[j][i].fu[4+5*m]*T0;
                }
                else { //set N mass, & temperature that is not calculated 
                    mass=ms[sl+5];
                    Tn0=uu[j][i].fu[4]*T0;   //using TO for TN
                }
                Hn[m]=kbmp*Tn0/(mass*grr);

                x=r/Hn[m];
                y=sqrt(0.5*x)*tcos;
                if (y > 13.0) y=13.0;

                /* calculate Chapman integration along path from Sun to observation point */
                if (zenith[j] < pih) {
                    Ch[m]=sqrt(pih*x)*exp(y*y)*cerfc(y);
                }
                else {
                    dc=x*(1.0-tsin);
                    if(dc > 100.0) dc=100.0;
                    Ch[m]=sqrt(pi2*x)*(sqrt(tsin)*exp(dc)-0.5*exp(y*y)*cerfc(y));
                }
            }

            Nnn[0]=exp(uu[j][i].fu[0]);    //normalized O density
            Nnn[1]=exp(uu[j][i].fu[5]);    //normalized O2 density
            Nnn[2]=exp(uu[j][i].fu[10]);    //normalized N2 density
            Nnn[3]=exp(uu[j][i].fu[21]);    //normalized N density

            //finally, determine optical depth
            for (l = 0; l < 37; l++) {
                tao[l]=0.0;
                for (m = 0; m < 4; m++) tao[l]=tao[l]+segabs[l][m]*Ch[m]*Nnn[m]*Hn[m];
                tao[l]=tao[l]*n0;
            }

            // ionization rates normalized by n0/t0
            for (m = 0; m < 4; m++) {
                qi[yj][xi][m]=0.0;
                if(m < 3) Qeuv[yj][xi][m]=0.0;
            }
            qi[yj][xi][4]=0.0;
            qit[yj][xi]=0.0;

            for (l = 0; l < 37; l++) {
                if (tao[l] > 100.0) tao[l]=100.0;
                euvf=euvflux[l]*exp(-tao[l])*t0;

                /* normalized ionization rates over all wavelengths for
                   O  + hv --> O+ + e;  O2 + hv --> O+ + O + e;
                   O2 + hv --> O2+ + e; N2 + hv --> N2+ + e */
                for (m = 0; m < 4; m++) {
                    if (m < 2) qi[yj][xi][m]= qi[yj][xi][m]+segion[l][m]*Nnn[m]*euvf;
                    else qi[yj][xi][m]= qi[yj][xi][m]+segion[l][m]*Nnn[m-1]*euvf;
                }

                /* normalized EUV heating of O, O2, N2 with heating efficiency of 0.45 */
                for (m = 0; m < 3; m++) 
                    Qeuv[yj][xi][m]=Qeuv[yj][xi][m]+segabs[l][m]*Nnn[m]*n0*euvf*pene[l]/p0;

                /* normalized ionization rates over all species for individual
                 * wavelength bins. Used to calculate local photoelectron heating
                 * rate for the wavelength 0-55 nm */
                if (l < 16) {
                    qib[yj][xi][l]=( segion[l][0]*Nnn[0]+segion[l][1]*Nnn[1]
                                    +segion[l][2]*Nnn[1]+segion[l][3]*Nnn[2])*euvf;
                    if (qib[yj][xi][l] <= 1.0e-30) qib[yj][xi][l]=1.0e-30;
                }
                /* total normalized ionization rate in 55-105 nm used to evaluate
                 * photoelectron heating rate for the wavelength 55-105 nm */
                else qit[yj][xi]= qit[yj][xi]+( segion[l][0]*Nnn[0]+segion[l][1]*Nnn[1]
                                 +segion[l][2]*Nnn[1]+segion[l][3]*Nnn[2])*euvf;
            }

            /* total ionization rate. used to evaluate photoelectron heating rate for
             * the wavelength 55-105 nm */
            for (m = 0; m < 4; m++) qit[yj][xi]=qit[yj][xi]+qi[yj][xi][m];
            if (qit[yj][xi] < 1.0e-30) qit[yj][xi]=1.0e-30;

            /* night time photoionization rates */
            if (zenith[j] >= pih && i <= i280) {
                for (m = 0; m < 4; m++) {A0[m]=0.0; Ih[m]=0.0;}

                for (l = i; l < i280; l++) {
                    Nnn[0]=exp(uu[j][l+1].fu[0]);  //normalized O density
                    Nnn[1]=exp(uu[j][l+1].fu[5]);  //normalized O2 density
                    Nnn[2]=exp(uu[j][l+1].fu[10]);  //normalized N2 density
                    Nnn[3]=exp(uu[j][l+1].fu[20]);  //normalized NO density

                    /* altitude attenuation at 4 lines (He II 303.78 A, He I 584.3 A,
                      Ly-beta 1025.67 A, Ly-alpha 1215.7 A) */
                    for (m = 0; m < 4; m++) {
                        A0[m]=A0[m]+0.5*( (Nnn[0]+exp(uu[j][i].fu[0]))*segabsn[m][0]
                                         +(Nnn[1]+exp(uu[j][i].fu[5]))*segabsn[m][1]
                                         +(Nnn[2]+exp(uu[j][i].fu[10]))*segabsn[m][2]
                                         +(Nnn[3]+exp(uu[j][i].fu[20]))*segabsn[m][3])
                                       *(rr[l+1]-rr[l])*r0n0;

                    }
                }
                if (i == i280) for (m = 0; m < 4; m++) A0[m]=0.0;

                //solar zenith dependence of nighttime radiation
                if (tcos != 0.0) GHy=3.62+3.4*cos(zenith[j])/pow(tcos, 0.37)/GHy120;
                else GHy=3.96/GHy120;

                GH1=exp(0.135*(90.0-zenith[j]*deg))/GH1120;
                GH2=2.2-0.01*zenith[j]*deg;

                /* altitude variation factor at four lines (He II 303.78 A, He I 584.3 A,
                   Ly-beta 1025.67 A, Ly-alpha 1215.7 A) */
                for (m = 0; m < 4; m++) {
                    for (l = 0; l < 4; l++) {
                        Adm=A0[l]/dem[m];
                        if (Adm > 130.0) Adm=130.0;
                        Ih[l]=Ih[l]+exp(-Adm);
                    }
                }

                //solar activity dependence of photon flux for these four lines
                f107p=(f107+120.0)/195.0;

                //solar photo flux scattered into night sector at 4 lines (photons m^-2 s^-1))
                Fph[0]=GHy*f107p*euvfluxn[0]*Ih[0];
                Fph[1]=GHy*f107p*euvfluxn[1]*Ih[1];
                Fph[2]=GH1*f107p*euvfluxn[2]*Ih[2];
                Fph[3]=GH2*f107p*euvfluxn[3]*Ih[3];

                //recover densities at rr[i]
                Nnn[0]=exp(uu[j][i].fu[0]);   //normalized O density
                Nnn[1]=exp(uu[j][i].fu[5]);   //normalized O2 density
                Nnn[2]=exp(uu[j][i].fu[10]);   //normalized N2 density
                Nnn[3]=exp(uu[j][i].fu[20]);   //normalized NO density

                //add night time photoionization rates (normalized) for O+, O2+, N2+, & NO+
                qi[yj][xi][0]=qi[yj][xi][0]+Nnn[0]*( Fph[0]*segionn[0][0]+Fph[1]*segionn[1][0]
                                                    +Fph[2]*segionn[2][0]+Fph[3]*segionn[3][0])*t0;
                qi[yj][xi][2]=qi[yj][xi][2]+Nnn[1]*( Fph[0]*segionn[0][1]+Fph[1]*segionn[1][1]
                                                    +Fph[2]*segionn[2][1]+Fph[3]*segionn[3][1])*t0;
                qi[yj][xi][3]=qi[yj][xi][3]+Nnn[2]*( Fph[0]*segionn[0][2]+Fph[1]*segionn[1][2]
                                                    +Fph[2]*segionn[2][2]+Fph[3]*segionn[3][2])*t0;
                qi[yj][xi][4]=qi[yj][xi][4]+Nnn[3]*( Fph[0]*segionn[0][3]+Fph[1]*segionn[1][3]
                                                    +Fph[2]*segionn[2][3]+Fph[3]*segionn[3][3])*t0;

            }
          }
        }
    }

    return;
}
