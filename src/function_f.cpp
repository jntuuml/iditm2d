#include <ctime>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "param.h"
#include "funcdef.h"

#undef __FUNCT__
#define __FUNCT__ "function_f"

PetscScalar function_f(int j,int i,int ir,int yj,int xi,Field **xx,Fieldu **uu)
{
    int    l, m, jp, ip, jm, im, im2;
    double ne[6], nve=0.0, nuv, ve1[6], ve2[6], ve3[6], ncs;
    double vd[3], dJ, dbeta, dvr, ns, vep[2], vip[2], vnp[2], nup[2];
    double ui[4], uj[4], Qi[4], Qj[4];
    int    s0, s1, s2, s4, s5, s6, s7, s8;
    PetscScalar ff=0.0;

    const double two3rd=0.6666666666666667;
    const double one6th=0.1666666666666667;
    const double nc=1.0e8/n0;

  /* Compute function fpart over the locally owned part of the grid */
    jm=j-1; jp=j+1;;

    im=i-1; ip=i+1; im2=i-2;;

    if (ir==0) {
        ne[0]=0.0; ne[2]=0.0;
        ve3[0]=0.0; ve3[2]=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ns=xx[j][i].fx[4+s0];
            ne[0]=ne[0]+ns;
            ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
            ns=xx[jp][i].fx[4+s0];
            ne[2]=ne[2]+ns;
            ve3[2]=ve3[2]+ns*xx[jp][i].fx[7+s0];
        }
        ve3[0]=ve3[0]/ne[0];
        ve3[2]=ve3[2]/ne[2];

        nve=0.0;
        for (m = 0; m < slm; m++) {
            if (m < sl) {
                s0=5*m;
                nve=nve+( uu[jp][i].nu[m]*(ve3[2]-xx[jp][i].fx[7+s0])
                         -uu[j][i].nu[m]*(ve3[0]-xx[j][i].fx[7+s0]));
            }
            else {
                s0=5*(m-sl);
                nve=nve+( uu[jp][i].nu[m]*(ve3[2]-uu[jp][i].fu[3+s0])
                         -uu[j][i].nu[m]*(ve3[0]-uu[j][i].fu[3+s0]));
            }
        }
        ff=-dt2*nve/(h2[i]*Omegae);
    }

    else if (ir==1) {
        ne[0]=0.0; ne[1]=0.0;
        ve3[0]=0.0; ve3[1]=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ns=xx[j][i].fx[4+s0];
            ne[0]=ne[0]+ns;
            ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
            ns=xx[j][ip].fx[4+s0];
            ne[1]=ne[1]+ns;
            ve3[1]=ve3[1]+ns*xx[j][ip].fx[7+s0];
        }
        ve3[0]=ve3[0]/ne[0];
        ve3[1]=ve3[1]/ne[1];

        nve=0.0;
        for (m = 0; m < slm; m++) {
            if (m < sl) {
                s0=5*m;
                nve=nve+( uu[j][ip].nu[m]*(ve3[1]-xx[j][ip].fx[7+s0])
                         -uu[j][i].nu[m]*(ve3[0]-xx[j][i].fx[7+s0]));
            }
            else {
                s0=5*(m-sl);
                nve=nve+( uu[j][ip].nu[m]*(ve3[1]-uu[j][ip].fu[3+s0])
                         -uu[j][i].nu[m]*(ve3[0]-uu[j][i].fu[3+s0]));
            }
        }
        ff=dt2*nve/(h1h[i]*Omegae);
    }

    else if (ir==2) {
        ve1[0]=0.0; ve1[1]=0.0; ve1[2]=0.0; ve1[5]=0.0;
        ve2[0]=0.0; ve2[1]=0.0; ve2[2]=0.0; ve2[5]=0.0;
        ne[0]=0.0; ne[1]=0.0; ne[2]=0.0; ne[5]=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ns=xx[j][i].fx[4+s0];
            ne[0]=ne[0]+ns;
            ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
            ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
            ns=xx[j][ip].fx[4+s0];
            ne[1]=ne[1]+ns;
            ve1[1]=ve1[1]+ns*xx[j][ip].fx[5+s0];
            ve2[1]=ve2[1]+ns*xx[j][ip].fx[6+s0];
            ns=xx[jp][i].fx[4+s0];
            ne[2]=ne[2]+ns;
            ve1[2]=ve1[2]+ns*xx[jp][i].fx[5+s0];
            ve2[2]=ve2[2]+ns*xx[jp][i].fx[6+s0];
            ns=xx[jp][ip].fx[4+s0];
            ne[5]=ne[5]+ns;
            ve1[5]=ve1[5]+ns*xx[jp][ip].fx[5+s0];
            ve2[5]=ve2[5]+ns*xx[jp][ip].fx[6+s0];
        }
        ve1[0]=ve1[0]/ne[0];
        ve1[1]=ve1[1]/ne[1];
        ve1[2]=ve1[2]/ne[2];
        ve1[5]=ve1[5]/ne[5];
        vep[0]=ve2[5]/ne[5]+ve2[1]/ne[1];       //v_{e2,j+1/2,i+1}
        vep[1]=ve2[2]/ne[2]+ve2[0]/ne[0];       //v_{e2,j+1/2,i}

        nve=0.0;
        for (l = 0; l < slm; l++) {
          nup[0]=uu[jp][ip].nu[l]+uu[j][ip].nu[l];   //nue_{j+1/2,i+1}
          nup[1]=uu[jp][i].nu[l]+uu[j][i].nu[l];     //nue_{j+1/2,i}

          if (l < sl) {
            s0=5*l;
            s5=5+s0; s6=6+s0;

            vip[0]=xx[jp][ip].fx[s6]+xx[j][ip].fx[s6];   //v_{s2,j+1/2,i+1}
            vip[1]=xx[jp][i].fx[s6]+xx[j][i].fx[s6];     //v_{s2,j+1/2,i}

            nve=nve+0.25*( (nup[0]*(vep[0]-vip[0])-nup[1]*(vep[1]-vip[1]))/h1h[i]
                          +0.25/rh[i]*(nup[0]+nup[1])*(vep[0]+vep[1]-(vip[0]+vip[1]))
                          -( (uu[jp][ip].nu[l]+uu[jp][i].nu[l])
                             *((ve1[5]+ve1[2])-(xx[jp][ip].fx[s5]+xx[jp][i].fx[s5]))
                            -(uu[j][ip].nu[l]+uu[j][i].nu[l])
                             *((ve1[1]+ve1[0])-(xx[j][ip].fx[s5]+xx[j][i].fx[s5])))/h2h[i]);
          }
          else {
            s0=5*(l-sl);
            s1=1+s0; s2=2+s0;
            vnp[0]=uu[jp][ip].fu[s2]+uu[j][ip].fu[s2];   //v_{n2,j+1/2,i+1}
            vnp[1]=uu[jp][i].fu[s2]+uu[j][i].fu[s2];     //v_{n2,j+1/2,i}

            nve=nve+0.25*( (nup[0]*(vep[0]-vnp[0])-nup[1]*(vep[1]-vnp[1]))/h1h[i]
                          +0.25/rh[i]*(nup[0]+nup[1])*(vep[0]+vep[1]-(vnp[0]+vnp[1]))
                          -( (uu[jp][ip].nu[l]+uu[jp][i].nu[l])
                             *((ve1[5]+ve1[2])-(uu[jp][ip].fu[s1]+uu[jp][i].fu[s1]))
                            -(uu[j][ip].nu[l]+uu[j][i].nu[l])
                             *((ve1[1]+ve1[0])-(uu[j][ip].fu[s1]+uu[j][i].fu[s1])))/h2h[i]);
          }
        }
        ff=-dt2*nve/Omegae;
    }

// ---- electron temperature equation for Te
    else if (ir==3) {
        ve1[0]=0.0; ve1[1]=0.0; ve1[3]=0.0;
        ve2[0]=0.0; ve2[2]=0.0; ve2[4]=0.0;
        ve3[0]=0.0;
        for (l=0; l<5; l++) ne[l]=0.0;
        for (l = 0; l < sl; l++) {
            s0=5*l;
            ns=xx[j][i].fx[4+s0];
            ne[0]=ne[0]+ns;
            ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
            ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
            ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
            if (i<N1) {
                ns=xx[j][ip].fx[4+s0];
                ne[1]=ne[1]+ns;
                ve1[1]=ve1[1]+ns*xx[j][ip].fx[5+s0];
            }
            ns=xx[jp][i].fx[4+s0];
            ne[2]=ne[2]+ns;
            ve2[2]=ve2[2]+ns*xx[jp][i].fx[6+s0];
            ns=xx[j][im].fx[4+s0];
            ne[3]=ne[3]+ns;
            ve1[3]=ve1[3]+ns*xx[j][im].fx[5+s0];
            ns=xx[jm][i].fx[4+s0];
            ne[4]=ne[4]+ns;
            ve2[4]=ve2[4]+ns*xx[jm][i].fx[6+s0];
        }
        ve1[0]=ve1[0]/ne[0];
        if (i<N1) ve1[1]=ve1[1]/ne[1];
        ve1[3]=ve1[3]/ne[3];
        ve2[0]=ve2[0]/ne[0];
        ve2[2]=ve2[2]/ne[2];
        ve2[4]=ve2[4]/ne[4];
        ve3[0]=ve3[0]/ne[0];

        /* Joule heating resulted from electron collisions with ions & neutrals */
        nuv=0.0; nve=0.0;
        for (l = 0; l < slm; l++) {
          if (l < sl) {
            s0=5*l;
            vd[0]=xx[j][i].fx[5+s0]-ve1[0];
            vd[1]=xx[j][i].fx[6+s0]-ve2[0];
            vd[2]=xx[j][i].fx[7+s0]-ve3[0];
            nve=nve+uu[j][i].nu[l]*(xx[j][i].fx[8+s0]-xx[j][i].fx[3])/ms[l];
          }
          else {
            s0=5*(l-sl);
            vd[0]=uu[j][i].fu[1+s0]-ve1[0];
            vd[1]=uu[j][i].fu[2+s0]-ve2[0];
            vd[2]=uu[j][i].fu[3+s0]-ve3[0];
            nve=nve+uu[j][i].nu[l]*(uu[j][i].fu[4+s0]-xx[j][i].fx[3])/ms[l];
          }
          nuv=nuv+uu[j][i].nu[l]*(vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2]);
        }

        uj[0]=(ve2[0]+ve2[4]);
        uj[1]=(ve2[2]+ve2[0]);
        if (uj[0] >= 0.0) Qj[0]=xx[jm][i].fx[3];
        else Qj[0]=xx[j][i].fx[3];
        if (uj[1] >= 0.0) Qj[1]=xx[j][i].fx[3];
        else Qj[1]=xx[jp][i].fx[3];

        if (i<N1) {
          //divergence of field-aligned current density (factor 0.5 not included)
          dJ= (uu[j][ip].Jpar[0]-uu[j][im].Jpar[0])/(2.0*h1[i])+uu[j][i].Jpar[0]/rr[i]
             +(uu[jp][i].Jpar[1]-uu[jm][i].Jpar[1])/(2.0*h2[i]);

          //gradient of thermal-electric conductivity dot product field-aligned current
          dbeta= (uu[j][ip].betae-uu[j][im].betae)/(2.0*h1[i])*uu[j][i].Jpar[0]
                +(uu[jp][i].betae-uu[jm][i].betae)/(2.0*h2[i])*uu[j][i].Jpar[1];

          ui[0]=(ve1[0]+ve1[3]);
          ui[1]=(ve1[1]+ve1[0]);
          if (ui[0] >= 0.0) Qi[0]=xx[j][im].fx[3];
          else Qi[0]=xx[j][i].fx[3];
          if (ui[1] >= 0.0) Qi[1]=xx[j][i].fx[3];
          else Qi[1]=xx[j][ip].fx[3];

          ff=dt2*( 0.5*(Qi[1]*ui[1]-Qi[0]*ui[0])/h1[i]+xx[j][i].fx[3]*ve1[0]/rr[i]
                  +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]
                  -one6th*xx[j][i].fx[3]
                         *((ve1[1]-ve1[3])/h1[i]+2.0*ve1[0]/rr[i]+(ve2[2]-ve2[4])/h2[i])
                  -two3rd/(ne[0]+nc*nc/ne[0])*(uu[j][i].betae*dJ+dbeta)
                  -2.0*me*nve-two3rd*(me*nuv+(Qe[yj][xi]-Ce[yj][xi])/ne[0]));
        }
        else {
          dJ= (3.0*uu[j][i].Jpar[0]-4.0*uu[j][im].Jpar[0]+uu[j][im2].Jpar[0])/(2.0*h1[i])
             +uu[j][i].Jpar[0]/rr[i]+(uu[jp][i].Jpar[1]-uu[jm][i].Jpar[1])/(2.0*h2[i]);

          dbeta= (3.0*uu[j][i].betae-4.0*uu[j][im].betae+uu[j][im2].betae)
                 /(2.0*h1[i])*uu[j][i].Jpar[0]
                +(uu[jp][i].betae-uu[jm][i].betae)/(2.0*h2[i])*uu[j][i].Jpar[1];

          ve1[5]=0.0; ncs=0.0;
          for (l=0; l<sl; l++) {
              s0=5*l;
              ns=xx[j][im2].fx[4+s0];
              ncs=ncs+ns;
              ve1[5]=ve1[5]+ns*xx[j][im2].fx[5+s0];
          }
          ve1[5]=ve1[5]/ncs;

          ff=dt2*( 0.5*( 3.0*xx[j][i].fx[3]*ve1[0]-4.0*xx[j][im].fx[3]*ve1[3]
                        +xx[j][im2].fx[3]*ve1[5])/h1[i]
                  +xx[j][i].fx[3]*ve1[0]/rr[i]
                  +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]
                  -one6th*xx[j][i].fx[3]
                         *( (3.0*ve1[0]-4.0*ve1[3]+ve1[5])/h1[i]
                           +2.0*ve1[0]/rr[i]+(ve2[2]-ve2[4])/h2[i])
                  -two3rd/(ne[0]+nc*nc/ne[0])*(uu[j][i].betae*dJ+dbeta)
                  -2.0*me*nve-two3rd*(me*nuv+(Qe[yj][xi]-Ce[yj][xi])/ne[0]));
        }
    }

/****************************************************************************************************
* elements from difference equations of densities, velocity components and temperatures */ 
// ---- difference equations of ns
    else if ((ir-4) % 5 == 0) {
        m=(ir-4)/5;
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0;

        uj[0]=(xx[j][i].fx[s6]+xx[jm][i].fx[s6]);
        uj[1]=(xx[jp][i].fx[s6]+xx[j][i].fx[s6]);
        if (uj[0] >= 0.0) Qj[0]=xx[jm][i].fx[s4];
        else Qj[0]=xx[j][i].fx[s4];
        if (uj[1] >= 0.0) Qj[1]=xx[j][i].fx[s4];
        else Qj[1]=xx[jp][i].fx[s4];

        if (i<N1) {
            ui[0]=(xx[j][i].fx[s5]+xx[j][im].fx[s5]);
            ui[1]=(xx[j][ip].fx[s5]+xx[j][i].fx[s5]);
            if (ui[0] >= 0.0) Qi[0]=xx[j][im].fx[s4];
            else Qi[0]=xx[j][i].fx[s4];
            if (ui[1] >= 0.0) Qi[1]=xx[j][i].fx[s4];
            else Qi[1]=xx[j][ip].fx[s4];

            ff=dt2*( 0.5*(Qi[1]*ui[1]-Qi[0]*ui[0])/h1[i]
                    +xx[j][i].fx[s4]*xx[j][i].fx[s5]/rr[i]
                    +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]);
        }
        else {
            ff=dt2*( 0.5*( 3.0*xx[j][i].fx[s4]*xx[j][i].fx[s5]
                          -4.0*xx[j][im].fx[s4]*xx[j][im].fx[s5]
                          +xx[j][im2].fx[s4]*xx[j][im2].fx[s5])/h1[i]
                    +xx[j][i].fx[s4]*xx[j][i].fx[s5]/rr[i]
                    +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]);
        }
    }

// ---- difference equations of vs1
    else if ((ir-5) % 5 == 0) {
        m=(ir-5)/5;
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s8=8+s0;

        ne[0]=0.0; ne[1]=0.0; ne[3]=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ne[0]=ne[0]+xx[j][i].fx[4+s0];
            ne[1]=ne[1]+xx[j][ip].fx[4+s0];
            ne[3]=ne[3]+xx[j][im].fx[4+s0];
        }

        ui[0]=(xx[j][i].fx[s5]+xx[j][im].fx[s5]);
        ui[1]=(xx[j][ip].fx[s5]+xx[j][i].fx[s5]);
        if (ui[0] >= 0.0) Qi[0]=xx[j][im].fx[s4]*xx[j][im].fx[s5];
        else Qi[0]=xx[j][i].fx[s4]*xx[j][i].fx[s5];
        if (ui[1] >= 0.0) Qi[1]=xx[j][i].fx[s4]*xx[j][i].fx[s5];
        else Qi[1]=xx[j][ip].fx[s4]*xx[j][ip].fx[s5];

        uj[0]=(xx[j][i].fx[s6]+xx[jm][i].fx[s6]);
        uj[1]=(xx[jp][i].fx[s6]+xx[j][i].fx[s6]);
        if (uj[0] >= 0.0) Qj[0]=xx[jm][i].fx[s4]*xx[jm][i].fx[s5];
        else Qj[0]=xx[j][i].fx[s4]*xx[j][i].fx[s5];
        if (uj[1] >= 0.0) Qj[1]=xx[j][i].fx[s4]*xx[j][i].fx[s5];
        else Qj[1]=xx[jp][i].fx[s4]*xx[jp][i].fx[s5];

        ff=dt2*( 0.5*(Qi[1]*ui[1]-Qi[0]*ui[0])/h1[i]
                +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]
                +xx[j][i].fx[s4]*( xx[j][i].fx[s5]*xx[j][i].fx[s5]
                                  -xx[j][i].fx[s6]*xx[j][i].fx[s6])/rr[i]
                +0.5/(ms[m]*h1[i])
                    *( xx[j][i].fx[s4]*(xx[j][ip].fx[s8]-xx[j][im].fx[s8])
                      +xx[j][i].fx[s8]*(xx[j][ip].fx[s4]-xx[j][im].fx[s4])
                      +xx[j][i].fx[s4]/ne[0]
                                      *( ne[0]*(xx[j][ip].fx[3]-xx[j][im].fx[3])
                                        +xx[j][i].fx[3]*(ne[1]-ne[3])))
                +xx[j][i].fx[s4]*gr[i]);
    }

// ---- difference equations of vs2
    else if ((ir-6) % 5 == 0) {
        m=(ir-6)/5;
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s8=8+s0;

        ne[0]=0.0; ne[2]=0.0; ne[4]=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ne[0]=ne[0]+xx[j][i].fx[4+s0];
            ne[2]=ne[2]+xx[jp][i].fx[4+s0];
            ne[4]=ne[4]+xx[jm][i].fx[4+s0];
        }

        if (i<N1) {
            ui[0]=(xx[j][i].fx[s5]+xx[j][im].fx[s5]);
            ui[1]=(xx[j][ip].fx[s5]+xx[j][i].fx[s5]);
            if (ui[0] >= 0.0) Qi[0]=xx[j][im].fx[s4]*xx[j][im].fx[s6];
            else Qi[0]=xx[j][i].fx[s4]*xx[j][i].fx[s6];
            if (ui[1] >= 0.0) Qi[1]=xx[j][i].fx[s4]*xx[j][i].fx[s6];
            else Qi[1]=xx[j][ip].fx[s4]*xx[j][ip].fx[s6];

            dvr=0.5*(Qi[1]*ui[1]-Qi[0]*ui[0])/h1[i];
        }
        else dvr=0.5*( 3.0*xx[j][i].fx[s4]*xx[j][i].fx[s6]*xx[j][i].fx[s5]
                      -4.0*xx[j][im].fx[s4]*xx[j][im].fx[s6]*xx[j][im].fx[s5]
                      +xx[j][im2].fx[s4]*xx[j][im2].fx[s6]*xx[j][im2].fx[s5])/h1[i];

        uj[0]=(xx[j][i].fx[s6]+xx[jm][i].fx[s6]);
        uj[1]=(xx[jp][i].fx[s6]+xx[j][i].fx[s6]);
        if (uj[0] >= 0.0) Qj[0]=xx[jm][i].fx[s4]*xx[jm][i].fx[s6];
        else Qj[0]=xx[j][i].fx[s4]*xx[j][i].fx[s6];
        if (uj[1] >= 0.0) Qj[1]=xx[j][i].fx[s4]*xx[j][i].fx[s6];
        else Qj[1]=xx[jp][i].fx[s4]*xx[jp][i].fx[s6];

        ff=dt2*( dvr+0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]
                +2.0*xx[j][i].fx[s4]*xx[j][i].fx[s5]*xx[j][i].fx[s6]/rr[i]
                +0.5/(ms[m]*h2[i])
                    *( xx[j][i].fx[s4]*(xx[jp][i].fx[s8]-xx[jm][i].fx[s8])
                      +xx[j][i].fx[s8]*(xx[jp][i].fx[s4]-xx[jm][i].fx[s4])
                      +xx[j][i].fx[s4]/ne[0]
                                      *( ne[0]*(xx[jp][i].fx[3]-xx[jm][i].fx[3])
                                        +xx[j][i].fx[3]*(ne[2]-ne[4]))));
    }

// ---- difference equations for vs3
    else if ((ir-7) % 5 == 0) {
        m=(ir-7)/5;
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0;

        if (i<N1) {
            ui[0]=(xx[j][i].fx[s5]+xx[j][im].fx[s5]);
            ui[1]=(xx[j][ip].fx[s5]+xx[j][i].fx[s5]);
            if (ui[0] >= 0.0) Qi[0]=xx[j][im].fx[s4]*xx[j][im].fx[s7];
            else Qi[0]=xx[j][i].fx[s4]*xx[j][i].fx[s7];
            if (ui[1] >= 0.0) Qi[1]=xx[j][i].fx[s4]*xx[j][i].fx[s7];
            else Qi[1]=xx[j][ip].fx[s4]*xx[j][ip].fx[s7];

            dvr=0.5*(Qi[1]*ui[1]-Qi[0]*ui[0])/h1[i];
        }
        else dvr=0.5*( 3.0*xx[j][i].fx[s4]*xx[j][i].fx[s7]*xx[j][i].fx[s5]
                      -4.0*xx[j][im].fx[s4]*xx[j][im].fx[s7]*xx[j][im].fx[s5]
                      +xx[j][im2].fx[s4]*xx[j][im2].fx[s7]*xx[j][im2].fx[s5])/h1[i];

        uj[0]=(xx[j][i].fx[s6]+xx[jm][i].fx[s6]);
        uj[1]=(xx[jp][i].fx[s6]+xx[j][i].fx[s6]);
        if (uj[0] >= 0.0) Qj[0]=xx[jm][i].fx[s4]*xx[jm][i].fx[s7];
        else Qj[0]=xx[j][i].fx[s4]*xx[j][i].fx[s7];
        if (uj[1] >= 0.0) Qj[1]=xx[j][i].fx[s4]*xx[j][i].fx[s7];
        else Qj[1]=xx[jp][i].fx[s4]*xx[jp][i].fx[s7];

        ff=dt2*( dvr+xx[j][i].fx[s4]*xx[j][i].fx[s7]*xx[j][i].fx[s5]/rr[i]
                +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]);
    }

// ---- difference equations for Ts. first calculate rates of heating (or cooling) by
//      electron, neutral & other ion species */
    else if ((ir-8) % 5 == 0) {
        m=(ir-8)/5;
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0; s8=8+s0;

        uj[0]=(xx[j][i].fx[s6]+xx[jm][i].fx[s6]);
        uj[1]=(xx[jp][i].fx[s6]+xx[j][i].fx[s6]);
        if (uj[0] >= 0.0) Qj[0]=xx[jm][i].fx[s8];
        else Qj[0]=xx[j][i].fx[s8];
        if (uj[1] >= 0.0) Qj[1]=xx[j][i].fx[s8];
        else Qj[1]=xx[jp][i].fx[s8];

        if (i < N1) {
            ui[0]=(xx[j][i].fx[s5]+xx[j][im].fx[s5]);
            ui[1]=(xx[j][ip].fx[s5]+xx[j][i].fx[s5]);
            if (ui[0] >= 0.0) Qi[0]=xx[j][im].fx[s8];
            else Qi[0]=xx[j][i].fx[s8];
            if (ui[1] >= 0.0) Qi[1]=xx[j][i].fx[s8];
            else Qi[1]=xx[j][ip].fx[s8];

            ff=dt2*( 0.5*(Qi[1]*ui[1]-Qi[0]*ui[0])/h1[i]
                    +xx[j][i].fx[s8]*xx[j][i].fx[s5]/rr[i]
                    +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]
                    -one6th*xx[j][i].fx[s8]
                           *( (xx[j][ip].fx[s5]-xx[j][im].fx[s5])/h1[i]
                             +2.0*xx[j][i].fx[s5]/rr[i]
                             +(xx[jp][i].fx[s6]-xx[jm][i].fx[s6])/h2[i]));
        }
        else {
            ff=dt2*( 0.5*( 3.0*xx[j][i].fx[s8]*xx[j][i].fx[s5]
                          -4.0*xx[j][im].fx[s8]*xx[j][im].fx[s5]
                          +xx[j][im2].fx[s8]*xx[j][im2].fx[s5])/h1[i]
                    +xx[j][i].fx[s8]*xx[j][i].fx[s5]/rr[i]
                    +0.5*(Qj[1]*uj[1]-Qj[0]*uj[0])/h2[i]
                    -one6th*xx[j][i].fx[s8]
                           *( ( 3.0*xx[j][i].fx[s5]-4.0*xx[j][im].fx[s5]
                               +xx[j][im2].fx[s5])/h1[i]
                             +2.0*xx[j][i].fx[s5]/rr[i]
                             +(xx[jp][i].fx[s6]-xx[jm][i].fx[s6])/h2[i]));

            //add heating and cooling rates
            /*if (i < ipc) {
              if (m == 0) {
                //heating due to exothermic reactions 
                ff=ff-two3rd*dt2/(cdt*ns[0])*Qch[yj][ci][m];
              }
            }*/
        }
    }

    return ff;
}
