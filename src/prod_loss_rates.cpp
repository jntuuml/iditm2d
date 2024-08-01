/*************************************************************************
  Routine prod_loss_rates
     Evaluate normalized production and loss rates for ion and neutral
 *   species considered. Also evaluate ion and neutral heating rates due
 *   to exothermic reactions

 * Jiannan Tu
 * 12/13/2013, 1/14/2014, 1/22/2014
*************************************************************************/
#include <cmath>

using namespace std;

#include "param.h"
#include "funcdef.h"

void prod_loss_rates(Field **xx,Fieldu **uu,int xs,int xm,int ys,int ym)
{
    int    i, j, m, sl1, sl2, sl3, sl4, sl5;
    double n0t0, n02t0, n0p0;
    double niO, niO2, niN2, niH, niNO, ne;
    double nO, nO2, nN2, nH, nNO, nN;
    double T1, T2, T300, Tr1, Tr2, Ti, Te;
    double kk[20], jj[5], k5, k6, k7, k11, k14, k15, ener;

/* energy released from chemical reactions (Joule) */
    const double cene[18]={2.4930232e-19, 1.74447536e-19, 6.985592e-19, 3.2044e-21,
                           1.490046e-19, 4.5069886e-19, 6.745262e-19, 3.140312e-19, 
                           1.12154e-19, 5.639744e-19, 1.0446344e-18, -3.2044e-21, 
                           8.203264e-19, 1.1199378e-18, 9.324804e-19, 4.40605e-19};

    int yj,xi;

    sl1=sl+1; sl2=sl+2; sl3=sl+3; sl4=sl+4; sl5=sl+5;

    /* normalized chemical and recombination reaction rates */
    n0t0=n0*t0;
    n02t0=n0*n0t0;
    n0p0=n0/p0;

    k5 =5.0e-22*n0t0;        // O2+ + N2 --> NO+ + NO
    k6 =4.5e-16*n0t0;        // O2+ + NO --> NO+ + O2
    k7 =1.5e-16*n0t0;        // O2+ + N  --> NO+ + O
    k11=3.3e-16*n0t0;        // N2+ + NO --> NO+ + N2
    k14=4.8e-18*n0t0;        // O+  + N2 --> N2+ + O + 3.02 eV
    k15=4.6e-16*n0t0;        // N2+ + O2 --> O+ + NO + 1.28 eV

    for (j=ys; j<ys+ym; j++) {
        yj=j-ys;

        for (i=xs; i<xs+xm; i++) {
          xi=i-xs;

          if (i <= ipc) {
            Ti=xx[j][i].fx[8]*T0;

            // O+ + O2--> O2+ + O (note kk[0] not used)
            T1=(0.667*Ti+0.333*uu[j][i].fu[9]*T0);
            T300=T1/300.0;
            kk[1]=( 2.820e-17-7.740e-18*T300+1.073e-18*T300*T300
                   -5.170e-20*pow(T300, 3.0)+9.650e-22*pow(T300, 4.0))*n0t0;

            // O+ + N2 --> NO+ + N
            T2=(0.6363*Ti+0.3637*uu[j][i].fu[14]*T0);
            T300=T2/300.0;
            if (T2 <= 1700.0) kk[2]=(1.533e-18-5.920e-19*T300+8.600e-20*T300*T300)*n0t0;
            else kk[2]=(2.730e-18-1.155e-18*T300+1.483e-19*T300*T300)*n0t0;

            // O+ + NO --> NO+ + O
            T300=Ti/300.0;
            if (Ti <= 1500.0) kk[3]=(8.36e-19-2.02e-19*T300+6.95e-20*T300*T300)*n0t0;
            else kk[3]=( 5.33e-19-1.64e-20*T300+4.72e-20*T300*T300
                        -7.05e-22*pow(T300,3.0))*n0t0;

            // O+ + H --> H+ + O
            kk[4]=2.5e-17*sqrt(uu[j][i].fu[19]*T0)*n0t0;

            Tr1=0.5*(xx[j][i].fx[23]+uu[j][i].fu[4])*T0;
            T300=Tr1/300.0;
            if (Tr1 <= 1500.0) {
                // N2+ + O --> O+ + N2
                kk[8]=1.0e-17*pow(T300, -0.23)*n0t0;

                // N2+ + O --> NO+ + N
                kk[9]=1.4e-16*pow(T300, -0.44)*n0t0;
            }
            else {
                kk[8]=3.6e-18*pow(T300, 0.41)*n0t0;
                kk[9]=5.2e-17*pow(T300, 0.2)*n0t0;
            }

            // N2+ + O2 --> O2+ + N2
            Tr2=0.5*(xx[j][i].fx[23]+uu[j][i].fu[9])*T0;
            kk[10]=5.0e-17*(300.0/Tr2)*n0t0;

            // H+ + O --> O+ + H
            kk[12]=2.2e-17*sqrt(xx[j][i].fx[13]*T0)*n0t0;

            // O + O + N2 --> O2 + N2
            kk[13]=9.59e-46*exp(480.0/(uu[j][i].fu[4]*T0))*n02t0;

            // O+ + e --> O
            Te=xx[j][i].fx[3]*T0;
            jj[0]=3.7e-18*pow(250.0/Te, 0.7)*n0t0;

            // O2+ + e --> O + O
            T300=300.0/Te;
            if (Te <= 1200.0) jj[1]=1.95e-13*pow(T300, 0.7)*n0t0;
            else jj[1]=7.38e-14*pow(1200.0/Te, 0.56)*n0t0;

            // N2+ + e --> N + N
            jj[2]=2.2e-13*pow(T300, 0.39)*n0t0;

            // NO+ + e --> N + O
            jj[3]=4.0e-13*sqrt(T300)*n0t0;

            // H+ + e --> H
            jj[4]=4.8e-18*pow(250.0/Te, 0.7)*n0t0;

            //normalized ion and neutral densities
            niO =xx[j][i].fx[4];
            niH =xx[j][i].fx[9];
            niO2=xx[j][i].fx[14];
            niN2=xx[j][i].fx[19];
            niNO=xx[j][i].fx[24];

            nO= exp(uu[j][i].fu[0]);
            nO2=exp(uu[j][i].fu[5]);
            nN2=exp(uu[j][i].fu[10]);
            nH =exp(uu[j][i].fu[15]);
            nNO=exp(uu[j][i].fu[20]);
            nN =exp(uu[j][i].fu[21]);

            //normalized electron density
            ne =niO+niH+niO2+niN2+niNO;

            /* O+ production and loss rates normalized by n0/t0 (reaction rates 
             * have been multiplied by n0t0) */
            Ps[yj][xi][0]=(kk[8]*niN2+kk[12]*niH)*nO;
            Ls[yj][xi][0]=(kk[1]*nO2+(kk[2]+k14)*nN2+kk[3]*nNO+kk[4]*nH+jj[0]*ne);
            if (Ls[yj][xi][0]*niO > Ps[yj][xi][0]*1.05) Ls[yj][xi][0]=1.05*Ps[yj][xi][0]/niO;

            //H+
            Ps[yj][xi][1]=kk[4]*niO*nH;
            Ls[yj][xi][1]=(kk[12]*nO+jj[4]*ne);
            if (Ls[yj][xi][1]*niH > Ps[yj][xi][1]*1.01) Ls[yj][xi][1]=1.01*Ps[yj][xi][1]/niH;

            // O2 +
            Ps[yj][xi][2]=(kk[1]*niO+kk[10]*niN2)*nO2;
            Ls[yj][xi][2]=(k5*nN2+k6*nNO+k7*nN+jj[1]*ne);
            if (Ps[yj][xi][2] > Ls[yj][xi][2]*1.05*niO2) Ps[yj][xi][2]=1.05*Ls[yj][xi][2]*niO2;

            //N2+
            Ps[yj][xi][3]=k14*niO*nN2;
            Ls[yj][xi][3]=((kk[8]+kk[9])*nO+k11*nNO+jj[2]*ne+kk[10]*nO2);
            if (Ls[yj][xi][3]*niN2 > Ps[yj][xi][3]*1.1) Ls[yj][xi][3]=1.1*Ps[yj][xi][3]/niN2;

            //NO+
            Ps[yj][xi][4]=( (kk[2]*nN2+kk[3]*nNO)*niO+(k5*nN2+k6*nNO+k7*nN)*niO2
                           +(kk[9]*nO+k11*nNO)*niN2)/100.0;
            Ls[yj][xi][4]=jj[3]*ne;

            //O
            Ps[yj][xi][5]= ( kk[1]*nO2+kk[3]*nNO+kk[4]*nH+jj[0]*ne)*niO
                            +k7*niO2*nN +(2.0*jj[1]*niO2+jj[3]*niNO)*ne+k14*niO*nN2;
            Ls[yj][xi][5]= ((kk[8]+kk[9])*niN2+kk[12]*niH)+2.0*kk[13]*nN2*nO;

            //O2
            Ps[yj][xi][6]=k6*niO2*nNO+kk[13]*nN2*nO*nO;
            Ls[yj][xi][6]=(kk[1]*niO+kk[10]*niN2)+k15*niN2;

            //N2
            Ps[yj][xi][7]=(kk[8]*nO+kk[10]*nO2+k11*nNO)*niN2;
            Ls[yj][xi][7]=(kk[2]*niO+k5*niO2)+k14*niO;

            //H
            Ps[yj][xi][8]=(kk[12]*nO+jj[4]*ne)*niH;
            Ls[yj][xi][8]=kk[4]*niO;

            //NO
            Ps[yj][xi][9]=k5*niO2*nN2;
            Ls[yj][xi][9]=(kk[3]*niO+k6*niO2+k11*niN2);

            //N
            Ps[yj][xi][10]=kk[2]*niO*nN2+kk[9]*niN2*nO+(2.0*jj[2]*niN2+jj[3]*niNO)*ne;
            Ls[yj][xi][10]=k7*niO2;

/* add photoionization rate & compute exothermic heating rates for neutrals */
            Ps[yj][xi][0]=Ps[yj][xi][0]+qi[yj][xi][0]+qi[yj][xi][1];
            Ps[yj][xi][2]=Ps[yj][xi][2]+qi[yj][xi][2];
            Ps[yj][xi][3]=Ps[yj][xi][3]+qi[yj][xi][3];
            Ps[yj][xi][4]=Ps[yj][xi][4]+qi[yj][xi][4];
            Ps[yj][xi][5]=Ps[yj][xi][5]+qi[yj][xi][1];
            Ls[yj][xi][5]=Ls[yj][xi][5]+qi[yj][xi][0]/nO;
            Ls[yj][xi][6]=Ls[yj][xi][6]+(qi[yj][xi][1]+qi[yj][xi][2])/nO2;
            Ls[yj][xi][7]=Ls[yj][xi][7]+qi[yj][xi][3]/nN2;
            Ls[yj][xi][9]=Ls[yj][xi][9]+qi[yj][xi][4]/nNO;

/* normalized exothermic heating rates */
            for (m = 0; m < sl+sm; m++) Qch[yj][xi][m]=0.0;

            //reaction k1 O+ + O2 --> O2+ + O
            ener=kk[1]*niO*nO2*cene[0]*ms[sl]/(ms[2]+ms[sl])*n0p0;
            Qch[yj][xi][1]=Qch[yj][xi][1]+ener;
            Qch[yj][xi][sl]=Qch[yj][xi][sl]+ener*ms[2]/ms[sl];

            //reaction k2  O+ + N2 --> NO+ + N
            ener=kk[2]*niO*nN2*cene[1]*ms[sl5]/(ms[4]+ms[sl5])*n0p0;
            Qch[yj][xi][4]=Qch[yj][xi][4]+ener;

            //reaction k3 O+ + NO --> NO+ + O
            ener=kk[3]*niO*nNO*cene[2]*ms[sl]/(ms[4]+ms[sl])*n0p0;
            Qch[yj][xi][4]=Qch[yj][xi][4]+ener;
            Qch[yj][xi][sl]=Qch[yj][xi][sl]+ener*ms[4]/ms[sl];

            //reaction k4 O+ + H --> H+ + O
            ener=kk[4]*niO*nH*cene[3]*ms[sl]/(ms[1]+ms[sl])*n0p0;
            Qch[yj][xi][3]=Qch[yj][xi][3]+ener;
            Qch[yj][xi][sl]=Qch[yj][xi][sl]+ener*ms[1]/ms[sl];

            //reaction k5 O2+ + N2 --> NO+ + NO
            ener=k5*niO2*nN2*cene[4]*ms[sl4]/(ms[4]+ms[sl4])*n0p0;
            Qch[yj][xi][4]=Qch[yj][xi][4]+ener;

            //reaction k6 O2+ + NO --> NO+ + O2
            ener=k6*niO2*nNO*cene[5]*ms[sl1]/(ms[4]+ms[sl1])*n0p0;
            Qch[yj][xi][4]=Qch[yj][xi][4]+ener;
            Qch[yj][xi][sl1]=Qch[yj][xi][sl1]+ener*ms[4]/ms[sl1];

            //reaction k7 O2+ + N --> NO+ + O
            ener=k7*niO2*nN*cene[6]*ms[sl]/(ms[4]+ms[sl])*n0p0;
            Qch[yj][xi][4]=Qch[yj][xi][4]+ener;
            Qch[yj][xi][sl]=Qch[yj][xi][sl]+ener*ms[4]/ms[sl];
 
            //reaction k8 N2+ + O --> O+ + N2
            ener=kk[8]*niN2*nO*cene[7]*ms[sl2]/(ms[0]+ms[sl2])*n0p0;
            Qch[yj][xi][0]=Qch[yj][xi][0]+ener;
            Qch[yj][xi][sl2]=Qch[yj][xi][sl2]+ener*ms[0]/ms[sl2];

            //reaction k9 N2+ + O --> NO+ + N
            ener=kk[9]*niN2*nO*cene[8]*ms[sl5]/(ms[4]+ms[sl5])*n0p0;
            Qch[yj][xi][4]=Qch[yj][xi][4]+ener;

            //reaction k10 N2+ + O2 --> O2+ + N2
            ener=kk[10]*niN2*nO2*cene[9]*ms[sl2]/(ms[2]+ms[sl2])*n0p0;
            Qch[yj][xi][1]=Qch[yj][xi][1]+ener;
            Qch[yj][xi][sl2]=Qch[yj][xi][sl2]+ener*ms[2]/ms[sl2];

            //reaction k11 N2+ + NO --> NO+ + N2
            ener=k11*niN2*nNO*cene[10]*ms[sl2]/(ms[4]+ms[sl2])*n0p0;
            Qch[yj][xi][4]=Qch[yj][xi][4]+ener;
            Qch[yj][xi][sl2]=Qch[yj][xi][sl2]+ener*ms[4]/ms[sl2];

            //reaction k12 H+ + O --> O+ + H
            ener=kk[12]*niH*nO*cene[11]*ms[sl3]/(ms[0]+ms[sl3])*n0p0;
            Qch[yj][xi][0]=Qch[yj][xi][0]-ener;
            Qch[yj][xi][sl3]=Qch[yj][xi][sl3]-ener*ms[0]/ms[sl3];

            //reaction k13 O + O + M --> O2 + M
            Qch[yj][xi][sl]=Qch[yj][xi][sl]+kk[13]*nO*nO*cene[12]*n0p0;

            //reaction j2 O2+ + e --> O + O
            Qch[yj][xi][sl]=Qch[yj][xi][sl]+jj[1]*niO2*ne*cene[13]*n0p0;

            //reaction j4 NO+ + e --> N + O
            Qch[yj][xi][sl]=Qch[yj][xi][sl]+jj[3]*niNO*ne*cene[15]*ms[sl5]
                                             /(ms[sl]+ms[sl5])*n0p0;
          }
        }
    }

    return;
}
