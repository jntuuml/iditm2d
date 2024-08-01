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

PetscScalar function_f(int, int, int, int, int, Field**, Fieldu**);

#undef __FUNCT__
#define __FUNCT__ "function_part"

PetscScalar function_part(int j,int i,int ir,int yj,int xi,Field **xx,
                          Field **xn,Field **xn1,Fieldu **uu)
{
    int    l, m, ll, jp, ip, jm, im, im2;
    double ne[6], nve=0.0, nuv, ve1[6], ve2[6], ve3[6];
    double botn, botp;
    double vd[3], vep[4], B01p[2], Td;
    double J1=0.0, J2=0.0, J3=0.0, ns, ba1=0.0, ba2=0.0, ba3=0.0, vd2;
    int    s0, s4, s5, s6, s7, s8, c0;
    PetscScalar ff=0.0;

    const double two3rd=0.6666666666666667;
    const double one6th=0.1666666666666667;
    const double nc=1.0e8/n0;

  /* Compute function over the locally owned part of the grid */
    jm=j-1; jp=j+1;;
    im=i-1; ip=i+1; im2=i-2;;

/* ---- difference equation of B_r */
    if (ir==0) {
        ne[0]=0.0; ne[2]=0.0;
        ve1[0]=0.0; ve1[2]=0.0;
        ve2[0]=0.0; ve2[2]=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ns=xx[j][i].fx[4+s0];
            ne[0]=ne[0]+ns;
            ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
            ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
            ns=xx[jp][i].fx[4+s0];
            ne[2]=ne[2]+ns;
            ve1[2]=ve1[2]+ns*xx[jp][i].fx[5+s0];
            ve2[2]=ve2[2]+ns*xx[jp][i].fx[6+s0];
        }
        ve1[0]=ve1[0]/ne[0]; ve1[2]=ve1[2]/ne[2];
        ve2[0]=ve2[0]/ne[0]; ve2[2]=ve2[2]/ne[2];

        if (i == 0) {
            botn=-xx[j][i].fx[1];
            botp=-xx[jp][i].fx[1];
        }
        else {
            botn=xx[j][im].fx[1];
            botp=xx[jp][im].fx[1];
        }

      ff= 3.0*xx[j][i].fx[0]-4.0*xn[j][i].fx[0]+xn1[j][i].fx[0]
         +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
         +( ( ve2[2]*(0.5*(xx[jp][i].fx[0]+xx[j][i].fx[0])+uu[jp][i].B0[0])
             -ve2[0]*(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]))
           -( ve1[2]*(0.5*(xx[jp][i].fx[1]+botp)+uu[jp][i].B0[1])
             -ve1[0]*(0.5*(xx[j][i].fx[1]+botn)+uu[j][i].B0[1])))*dt2/h2[i];
    }

/* ---- difference equation for B_theta */
    else if (ir==1) {
      if(i<N1) {
        ne[0]=0.0; ne[1]=0.0;
        ve1[0]=0.0; ve1[1]=0.0;
        ve2[0]=0.0; ve2[1]=0.0;
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
        }
        ve1[0]=ve1[0]/ne[0]; ve1[1]=ve1[1]/ne[1];
        ve2[0]=ve2[0]/ne[0]; ve2[1]=ve2[1]/ne[1];

        if (i == 0) botn=-xx[j][i].fx[1];
        else botn=xx[j][im].fx[1];

        ff= 3.0*xx[j][i].fx[1]-4.0*xn[j][i].fx[1]+xn1[j][i].fx[1]
           +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
           +( ( ve1[1]*(0.5*(xx[j][ip].fx[1]+xx[j][i].fx[1])+uu[j][ip].B0[1])
               -ve1[0]*(0.5*(xx[j][i].fx[1]+botn)+uu[j][i].B0[1]))
             -( ve2[1]*(0.5*(xx[j][ip].fx[0]+xx[jm][ip].fx[0])+uu[j][ip].B0[0])
               -ve2[0]*(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]))
            )*dt2/h1h[i];
      }
      else {
        ne[0]=0.0; ne[3]=0.0;
        ve1[0]=0.0; ve1[3]=0.0;
        ve2[0]=0.0; ve2[3]=0.0;
        for (l=0; l<sl; l++) {
            s0=5*l;
            ns=xx[j][i].fx[4+s0];
            ne[0]=ne[0]+ns;
            ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
            ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
            ns=xx[j][im].fx[4+s0];
            ne[3]=ne[3]+ns;
            ve1[3]=ve1[3]+ns*xx[j][im].fx[5+s0];
            ve2[3]=ve2[3]+ns*xx[j][im].fx[6+s0];
        }
        ve1[0]=ve1[0]/ne[0]; ve1[3]=ve1[3]/ne[3];
        ve2[0]=ve2[0]/ne[0]; ve2[3]=ve2[3]/ne[3];

        ff= 3.0*xx[j][i].fx[1]-4.0*xn[j][i].fx[1]+xn1[j][i].fx[1]
           +( ( ve1[0]*(0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])+uu[j][i].B0[1])
               -ve1[3]*(0.5*(xx[j][im].fx[1]+xx[j][im2].fx[1])+uu[j][im].B0[1]))
             -( ve2[0]*(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0])
               -ve2[3]*(0.5*(xx[j][im].fx[0]+xx[jm][im].fx[0])+uu[j][im].B0[0]))
            )*dt2/h1h[i];
      }
    }

/* ---- difference equation for B_z */
    else if (ir==2) {
      if (i < N1) {
        ne[0]=0.0; ne[1]=0.0; ne[2]=0.0; ne[5]=0.0;
        ve1[0]=0.0; ve1[1]=0.0; ve1[2]=0.0; ve1[5]=0.0;
        ve2[0]=0.0; ve2[1]=0.0; ve2[2]=0.0; ve2[5]=0.0;
        ve3[0]=0.0; ve3[1]=0.0; ve3[2]=0.0; ve3[5]=0.0;
        for (l=0; l<sl; l++) {
          s0=5*l;
          ns=xx[j][i].fx[4+s0];
          ne[0]=ne[0]+ns;
          ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
          ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
          ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];

          ns=xx[j][ip].fx[4+s0];
          ne[1]=ne[1]+ns;
          ve1[1]=ve1[1]+ns*xx[j][ip].fx[5+s0];
          ve2[1]=ve2[1]+ns*xx[j][ip].fx[6+s0];
          ve3[1]=ve3[1]+ns*xx[j][ip].fx[7+s0];

          ns=xx[jp][i].fx[4+s0];
          ne[2]=ne[2]+ns;
          ve1[2]=ve1[2]+ns*xx[jp][i].fx[5+s0];
          ve2[2]=ve2[2]+ns*xx[jp][i].fx[6+s0];
          ve3[2]=ve3[2]+ns*xx[jp][i].fx[7+s0];

          ns=xx[jp][ip].fx[4+s0];
          ne[5]=ne[5]+ns;
          ve1[5]=ve1[5]+ns*xx[jp][ip].fx[5+s0];
          ve2[5]=ve2[5]+ns*xx[jp][ip].fx[6+s0];
          ve3[5]=ve3[5]+ns*xx[jp][ip].fx[7+s0];
        }
        ve1[0]=ve1[0]/ne[0]; ve1[1]=ve1[1]/ne[1];
        ve1[2]=ve1[2]/ne[2]; ve1[5]=ve1[5]/ne[5];
        ve2[0]=ve2[0]/ne[0]; ve2[1]=ve2[1]/ne[1];
        ve2[2]=ve2[2]/ne[2]; ve2[5]=ve2[5]/ne[5];
        ve3[0]=ve3[0]/ne[0]; ve3[1]=ve3[1]/ne[1];
        ve3[2]=ve3[2]/ne[2]; ve3[5]=ve3[5]/ne[5];

        vep[0]=ve1[5]+ve1[1];          //v_{e1,j+1/2,i+1}
        vep[1]=ve1[2]+ve1[0];          //v_{e1,j+1/2,i}
        vep[2]=ve3[5]+ve3[1];          //v_{e3,j+1/2,i+1}
        vep[3]=ve3[2]+ve3[0];          //v_{e3,j+1/2,i}
        B01p[0]=0.5*(uu[jp][ip].B0[0]+uu[j][ip].B0[0]);   //B_{01,j+1/2,i+1}
        B01p[1]=0.5*(uu[jp][i].B0[0]+uu[j][i].B0[0]);     //B_{01,j+1/2,i}

        if (i == 0) botn=-xx[j][i].fx[2];
        else botn=xx[j][im].fx[2];

        ff=3.0*xx[j][i].fx[2]-4.0*xn[j][i].fx[2]+xn1[j][i].fx[2]
         +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
         +( ( 0.25*(vep[0]*(xx[j][ip].fx[2]+xx[j][i].fx[2])-vep[1]*(xx[j][i].fx[2]+botn))
             -0.5*(vep[2]*(xx[j][ip].fx[0]+B01p[0])-vep[3]*(xx[j][i].fx[0]+B01p[1])))/h1h[i]
           +0.25/rh[i]*( (vep[0]+vep[1])*xx[j][i].fx[2]
                        -0.5*(vep[2]+vep[3])
                            *(xx[j][ip].fx[0]+xx[j][i].fx[0]+(B01p[0]+B01p[1])))
           +( 0.25*( (ve2[5]+ve2[2])*(xx[jp][i].fx[2]+xx[j][i].fx[2])
                    -(ve2[1]+ve2[0])*(xx[j][i].fx[2]+xx[jm][i].fx[2]))
             -0.5*( (ve3[5]+ve3[2])
                    *(xx[jp][i].fx[1]+0.5*(uu[jp][ip].B0[1]+uu[jp][i].B0[1]))
                   -(ve3[1]+ve3[0])
                    *(xx[j][i].fx[1]+0.5*(uu[j][ip].B0[1]+uu[j][i].B0[1]))))/h2h[i]
           +( (xx[jp][ip].fx[3]+xx[jp][i].fx[3]-xx[j][ip].fx[3]-xx[j][i].fx[3])
              *(ne[5]+ne[1]-ne[2]-ne[0])
             -(xx[jp][ip].fx[3]+xx[j][ip].fx[3]-xx[jp][i].fx[3]-xx[j][i].fx[3])
              *(ne[5]+ne[2]-ne[1]-ne[0])
            )/(e*(ne[5]+ne[1]+ne[2]+ne[0])*hh[i]))*dt2;
      }
      else ff=xx[j][i].fx[2]-h1[im]/h1[i]*xx[j][im].fx[2];
    }

// ---- electron temperature equation for Te
    else if (ir==3) {
      if (i==0) ff=xx[j][i].fx[3]-xn[j][i].fx[3];
      else {
        ne[0]=0.0;
        for (l=0; l<sl; l++) ne[0]=ne[0]+xx[j][i].fx[4+5*l];

        if (i<N1) {
          ff= 3.0*xx[j][i].fx[3]-4.0*xn[j][i].fx[3]+xn1[j][i].fx[3]
               +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
             -one6th*dt2/(ne[0]+nc*nc/ne[0])
                    *( ((uu[j][ip].lamdae-uu[j][im].lamdae)/h12[i]+2.0*uu[j][i].lamdae*dh1[i])
                       *(xx[j][ip].fx[3]-xx[j][im].fx[3])
                      +(uu[jp][i].lamdae-uu[jm][i].lamdae)/h22[i]
                       *(xx[jp][i].fx[3]-xx[jm][i].fx[3])
                      +4.0*uu[j][i].lamdae
                       *( (xx[j][ip].fx[3]-2.0*xx[j][i].fx[3]+xx[j][im].fx[3])/h12[i]
                         +(xx[jp][i].fx[3]-2.0*xx[j][i].fx[3]+xx[jm][i].fx[3])/h22[i]));
        }
        else {
          ff=xx[j][i].fx[3]-xx[j][im].fx[3];

          /*ff= 3.0*xx[j][i].fx[3]-4.0*xn[j][i].fx[3]+xn1[j][i].fx[3]
               +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
             -one6th*dt2/(ne[0]+nc*nc/ne[0])
                    *( ( (3.0*uu[j][i].lamdae-4.0*uu[j][im].lamdae+uu[j][im2].lamdae)/h12[i]
                        +2.0*uu[j][i].lamdae*dh1[i])
                       *(3.0*xx[j][i].fx[3]-4.0*xx[j][im].fx[3]+xx[j][im2].fx[3])
                      +(uu[jp][i].lamdae-uu[jm][i].lamdae)/h22[i]
                       *(xx[jp][i].fx[3]-xx[jm][i].fx[3])
                      +4.0*uu[j][i].lamdae
                       *(xx[jp][i].fx[3]-2.0*xx[j][i].fx[3]+xx[jm][i].fx[3])/h22[i]);*/
        }
      }
    }

/****************************************************************************************************
* elements from difference equations of densities, velocity components and temperatures */ 
// ---- difference equations of ns
    else if ((ir-4) % 5 == 0) {
        m=(ir-4)/5;
        s0=5*m;
        s4=4+s0;

        if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
            ff=xx[j][i].fx[s4]-h1[im]/h1[i]*xx[j][im].fx[s4];
        }
        else {
            if (i==0 || (m<2 && zh[i]<=120.0e3)) {
              if (xn[j][i].fx[s4]<=1.0e3/n0) ff=xx[j][i].fx[s4]-xn[j][i].fx[s4];
              else if (i==0) {
                ff= xx[j][i].fx[s4]
                   +( xn1[j][i].fx[s4]-4.0*xn[j][i].fx[s4]
                     -dt2*Ps[yj][xi][m])/(3.0+dt2*Ls[yj][xi][m]);
              }
              else {
                ff= xx[j][i].fx[s4]
                   +( xn1[j][i].fx[s4]-4.0*xn[j][i].fx[s4]
                     +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                     -dt2*Ps[yj][xi][m])/(3.0+dt2*Ls[yj][xi][m]);
              }
            }
            else {
              if (i<=ipc) {
                ff= xx[j][i].fx[s4]
                   +( xn1[j][i].fx[s4]-4.0*xn[j][i].fx[s4]
                     +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                     -dt2*Ps[yj][xi][m])/(3.0+dt2*Ls[yj][xi][m]);
              }
              else
                ff= xx[j][i].fx[s4]
                   +( xn1[j][i].fx[s4]-4.0*xn[j][i].fx[s4]
                     +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                    )/3.0;
            }
        }
    }

// ---- difference equations of vs1
    else if ((ir-5) % 5 == 0) {
        m=(ir-5)/5; c0=9*(m+1);
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0;

        if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
            /*ne[0]=0.0; ve1[0]=0.0;
            for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
            }
            ff=xx[j][i].fx[s5]-ve1[0]/ne[0];*/
            ff=xx[j][i].fx[s5]-h1[im]/h1[i]*xx[j][im].fx[s5];
        }
        else {
            if (i==0) ff=xx[j][i].fx[s5]-h1[i]/h1[ip]*xx[j][ip].fx[s5];
            else if (i==N1) ff=xx[j][i].fx[s5]-h1[im]/h1[i]*xx[j][im].fx[s5];
            else {
              ne[0]=0.0; ve1[0]=0.0; ve2[0]=0.0; ve3[0]=0.0;
              for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
                ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
                ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
              }
              J1=0.5*(xx[j][i].fx[2]+xx[j][im].fx[2]-xx[jm][i].fx[2]-xx[jm][im].fx[2])/h2[i];
              ve1[0]=ve1[0]/ne[0]-J1/(e*ne[0]);
              ve2[0]=ve2[0]/ne[0];
              ve3[0]=ve3[0]/ne[0];

              //collision terms in ion & neutral momentum equations
              nuv=uu[j][i].nu[c0]*(xx[j][i].fx[s5]-ve1[0]);
              nve=uu[j][i].nu[0]*(ve1[0]-xx[j][i].fx[5]);

              for (l=1; l<slm; l++) {
                if (l < sl) {
                    if (l <= m) ll=l-1; else ll=l;
                    nuv=nuv+uu[j][i].nu[l+c0]*(xx[j][i].fx[s5]-xx[j][i].fx[5+5*ll]);
                    nve=nve+uu[j][i].nu[l]*(ve1[0]-xx[j][i].fx[5+5*l]);
                }
                else {
                    s0=5*(l-sl);
                    nuv=nuv+uu[j][i].nu[l+c0]*(xx[j][i].fx[s5]-uu[j][i].fu[1+s0]);
                    nve=nve+uu[j][i].nu[l]*(ve1[0]-uu[j][i].fu[1+s0]);
                }
              }

              ba2=0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])+uu[j][i].B0[1];
              ba3=0.25*(xx[j][i].fx[2]+xx[jm][i].fx[2]+xx[j][im].fx[2]+xx[jm][im].fx[2]);

              J2=0.5*(xx[j][im].fx[2]+xx[jm][im].fx[2]-xx[j][i].fx[2]-xx[jm][i].fx[2])/h1[i];
              J3= (xx[j][i].fx[1]-xx[j][im].fx[1])/h1[i]
                 +0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])/rr[i]
                 -(xx[j][i].fx[0]-xx[jm][i].fx[0])/h2[i];

              ff= xx[j][i].fx[s4]*xx[j][i].fx[s5]
                 +( xn1[j][i].fx[s4]*xn1[j][i].fx[s5]-4.0*xn[j][i].fx[s4]*xn[j][i].fx[s5]
                   +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                   +dt2*xx[j][i].fx[s4]*( nuv+me/ms[m]*nve
                                         -qms[m]*( (xx[j][i].fx[s6]-ve2[0])*ba3
                                                  -(xx[j][i].fx[s7]-ve3[0])*ba2)
                                         -(J2*ba3-J3*ba2)/(ne[0]*ms[m])))
                  /(3.0+dt2*Ls[yj][xi][m]);
            }
        }
    }

// ---- difference equations of vs2
    else if ((ir-6) % 5 == 0) {
        m=(ir-6)/5; c0=9*(m+1);
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0;

        if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
            /*ne[0]=0.0; ve2[0]=0.0;
            for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
            }
            ff=xx[j][i].fx[s6]-ve2[0]/ne[0];*/
            ff=xx[j][i].fx[s6]-h1[im]/h1[i]*xx[j][im].fx[s6];
        }
        else {
            if (i==0) ff=xx[j][i].fx[s6];
            else if (i==N1) {
              if (nvt==0 && btop[j]>-1.0e5) ff=xx[j][i].fx[s6]-btop[j];
              else ff=xx[j][i].fx[s6]-h1[im]/h1[i]*xx[j][im].fx[s6];
            }
            else {
              ne[0]=0.0; ve1[0]=0.0; ve2[0]=0.0; ve3[0]=0.0;
              for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
                ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
                ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
              }
              J2=0.5*(xx[j][im].fx[2]+xx[jm][im].fx[2]-xx[j][i].fx[2]-xx[jm][i].fx[2])/h1[i];
              ve1[0]=ve1[0]/ne[0];
              ve2[0]=ve2[0]/ne[0]-J2/(e*ne[0]);
              ve3[0]=ve3[0]/ne[0];

              nuv=uu[j][i].nu[c0]*(xx[j][i].fx[s6]-ve2[0]);
              nve=uu[j][i].nu[0]*(ve2[0]-xx[j][i].fx[6]);

              for (l=1; l<slm; l++) {
                if (l < sl) {
                    if (l <= m) ll=l-1; else ll=l;
                    nuv=nuv+uu[j][i].nu[l+c0]*(xx[j][i].fx[s6]-xx[j][i].fx[6+5*ll]);
                    nve=nve+uu[j][i].nu[l]*(ve2[0]-xx[j][i].fx[6+5*l]);
                }
                else {
                    s0=5*(l-sl);
                    nuv=nuv+uu[j][i].nu[l+c0]*(xx[j][i].fx[s6]-uu[j][i].fu[2+s0]);
                    nve=nve+uu[j][i].nu[l]*(ve2[0]-uu[j][i].fu[2+s0]);
                }
              }

              J1=0.5*(xx[j][i].fx[2]+xx[j][im].fx[2]-xx[jm][i].fx[2]-xx[jm][im].fx[2])/h2[i];
              J3= (xx[j][i].fx[1]-xx[j][im].fx[1])/h1[i]
                 +0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])/rr[i]
                 -(xx[j][i].fx[0]-xx[jm][i].fx[0])/h2[i];

              ba1=(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]);                        
              ba3=0.25*(xx[j][i].fx[2]+xx[jm][i].fx[2]+xx[j][im].fx[2]+xx[jm][im].fx[2]);

              ff= xx[j][i].fx[s4]*xx[j][i].fx[s6]
                 +( xn1[j][i].fx[s4]*xn1[j][i].fx[s6]-4.0*xn[j][i].fx[s4]*xn[j][i].fx[s6]
                   +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                   +dt2*xx[j][i].fx[s4]*( nuv+me/ms[m]*nve
                                         -qms[m]*( (xx[j][i].fx[s7]-ve3[0])*ba1
                                                  -(xx[j][i].fx[s5]-ve1[0])*ba3)
                                         -(J3*ba1-J1*ba3)/(ne[0]*ms[m])))
                  /(3.0+dt2*Ls[yj][xi][m]);
            }
        }
    }

// ---- difference equations for vs3
    else if ((ir-7) % 5 == 0) {
        m=(ir-7)/5; c0=9*(m+1);
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0;

        if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
            /*ne[0]=0.0; ve3[0]=0.0;
            for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
            }
            ff=xx[j][i].fx[s7]-ve3[0]/ne[0];*/
            ff=xx[j][i].fx[s7]-h1[im]/h1[i]*xx[j][im].fx[s7];
        }
        else {
            if (i==0) ff=xx[j][i].fx[s7];
            else if (i==N1) {
              if (nvt==1 && btop[j]>-1.0e5) ff=xx[j][i].fx[s7]-btop[j];
              else ff=xx[j][i].fx[s7]-h1[im]/h1[i]*xx[j][im].fx[s7];
            }
            else {
              ne[0]=0.0; ve1[0]=0.0; ve2[0]=0.0; ve3[0]=0.0;
              for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
                ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
                ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
              }
              J3= (xx[j][i].fx[1]-xx[j][im].fx[1])/h1[i]
                 +0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])/rr[i]
                 -(xx[j][i].fx[0]-xx[jm][i].fx[0])/h2[i];
              ve1[0]=ve1[0]/ne[0];
              ve2[0]=ve2[0]/ne[0];
              ve3[0]=ve3[0]/ne[0]-J3/(e*ne[0]);

              nuv=uu[j][i].nu[c0]*(xx[j][i].fx[s7]-ve3[0]);
              nve=uu[j][i].nu[0]*(ve3[0]-xx[j][i].fx[7]);

              for (l=1; l<slm; l++) {
                if (l < sl) {
                    if (l <= m) ll=l-1; else ll=l;
                    nuv=nuv+uu[j][i].nu[l+c0]*(xx[j][i].fx[s7]-xx[j][i].fx[7+5*ll]);
                    nve=nve+uu[j][i].nu[l]*(ve3[0]-xx[j][i].fx[7+5*l]);
                }
                else {
                    s0=5*(l-sl);
                    nuv=nuv+uu[j][i].nu[l+c0]*(xx[j][i].fx[s7]-uu[j][i].fu[3+s0]);
                    nve=nve+uu[j][i].nu[l]*(ve3[0]-uu[j][i].fu[3+s0]);
                }
              }

              J1=0.5*(xx[j][i].fx[2]+xx[j][im].fx[2]-xx[jm][i].fx[2]-xx[jm][im].fx[2])/h2[i];
              J2=0.5*(xx[j][im].fx[2]+xx[jm][im].fx[2]-xx[j][i].fx[2]-xx[jm][i].fx[2])/h1[i];

              ba1=(0.5*(xx[j][i].fx[0]+xx[jm][i].fx[0])+uu[j][i].B0[0]);                        
              ba2=0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])+uu[j][i].B0[1];

              ff= xx[j][i].fx[s4]*xx[j][i].fx[s7]
                 +( xn1[j][i].fx[s4]*xn1[j][i].fx[s7]-4.0*xn[j][i].fx[s4]*xn[j][i].fx[s7]
                   +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                   +dt2*xx[j][i].fx[s4]*( nuv+me/ms[m]*nve
                                         -qms[m]*( (xx[j][i].fx[s5]-ve1[0])*ba2
                                                  -(xx[j][i].fx[s6]-ve2[0])*ba1)
                                         -(J1*ba2-J2*ba1)/(ne[0]*ms[m])))
                  /(3.0+dt2*Ls[yj][xi][m]);
            }
        }
    }

// ---- difference equations for Ts. first calculate rates of heating (or cooling) by
//      electron, neutral & other ion species */
    else if ((ir-8) % 5 == 0) {
        m=(ir-8)/5; c0=9*(m+1);
        s0=5*m;
        s4=4+s0; s5=5+s0; s6=6+s0; s7=7+s0; s8=8+s0;

        if (zh[i]>=1300.0e3 && xn[j][i].fx[s4]*n0<=1.0e-3) {
            ff=xx[j][i].fx[s8]-xx[j][im].fx[s8];
        }
        else {
            if (i == 0) ff=xx[j][i].fx[s8]-xn[j][i].fx[s8];
            else {
              ne[0]=0.0; ve1[0]=0.0; ve2[0]=0.0; ve3[0]=0.0;
              for (l=0; l<sl; l++) {
                s0=5*l;
                ns=xx[j][i].fx[4+s0];
                ne[0]=ne[0]+ns;
                ve1[0]=ve1[0]+ns*xx[j][i].fx[5+s0];
                ve2[0]=ve2[0]+ns*xx[j][i].fx[6+s0];
                ve3[0]=ve3[0]+ns*xx[j][i].fx[7+s0];
              }
              J1=0.5*(xx[j][i].fx[2]+xx[j][im].fx[2]-xx[jm][i].fx[2]-xx[jm][im].fx[2])/h2[i];
              J2=0.5*(xx[j][im].fx[2]+xx[jm][im].fx[2]-xx[j][i].fx[2]-xx[jm][i].fx[2])/h1[i];
              J3= (xx[j][i].fx[1]-xx[j][im].fx[1])/h1[i]
                 +0.5*(xx[j][i].fx[1]+xx[j][im].fx[1])/rr[i]
                 -(xx[j][i].fx[0]-xx[jm][i].fx[0])/h2[i];
              ve1[0]=ve1[0]/ne[0]-J1/(e*ne[0]);
              ve2[0]=ve2[0]/ne[0]-J2/(e*ne[0]);
              ve3[0]=ve3[0]/ne[0]-J3/(e*ne[0]);

              //calculate rates of heating (or cooling) by electron, neutral
              // & other ion species
              vd[0]=ve1[0]-xx[j][i].fx[s5];
              vd[1]=ve2[0]-xx[j][i].fx[s6];
              vd[2]=ve3[0]-xx[j][i].fx[s7];
              nuv=uu[j][i].nu[c0]*( 2.0*(xx[j][i].fx[3]-xx[j][i].fx[s8])
                                   +two3rd*me*(vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2]));

              for (l = 1; l < slm; l++) {
                if (l <= m) ll=l-1; else ll=l;

                if (l < sl) {
                    s0=5*ll;
                    vd[0]=xx[j][i].fx[5+s0]-xx[j][i].fx[s5];
                    vd[1]=xx[j][i].fx[6+s0]-xx[j][i].fx[s6];
                    vd[2]=xx[j][i].fx[7+s0]-xx[j][i].fx[s7];
                    Td=xx[j][i].fx[8+s0]-xx[j][i].fx[s8];
                }
                else {
                    s0=5*(l-sl);
                    vd[0]=uu[j][i].fu[1+s0]-xx[j][i].fx[s5];
                    vd[1]=uu[j][i].fu[2+s0]-xx[j][i].fx[s6];
                    vd[2]=uu[j][i].fu[3+s0]-xx[j][i].fx[s7];
                    Td=uu[j][i].fu[4+s0]-xx[j][i].fx[s8];
                }
                vd2=vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2];
 
                nuv=nuv+ms[m]*uu[j][i].nu[l+c0]/(ms[m]+ms[ll])
                                               *(2.0*Td+two3rd*ms[ll]*vd2);
              }

              if (i < N1) {
                ff= 3.0*xx[j][i].fx[s8]-4.0*xn[j][i].fx[s8]+xn1[j][i].fx[s8]
                   +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                   -dt2*( one6th/(xx[j][i].fx[s4]+nc*nc/xx[j][i].fx[s4])
                            *( ( (uu[j][ip].lamdas[m]-uu[j][im].lamdas[m])/h12[i]
                                +2.0*uu[j][i].lamdas[m]*dh1[i])
                               *(xx[j][ip].fx[s8]-xx[j][im].fx[s8])
                              +(uu[jp][i].lamdas[m]-uu[jm][i].lamdas[m])/h22[i]
                               *(xx[jp][i].fx[s8]-xx[jm][i].fx[s8])
                              +4.0*uu[j][i].lamdas[m]
                             *( (xx[j][ip].fx[s8]-2.0*xx[j][i].fx[s8]+xx[j][im].fx[s8])/h12[i]
                               +(xx[jp][i].fx[s8]-2.0*xx[j][i].fx[s8]+xx[jm][i].fx[s8])/h22[i]))
                         +nuv);
              }
              else {
                ff=xx[j][i].fx[s8]-xx[j][im].fx[s8];

                /*ff= 3.0*xx[j][i].fx[s8]-4.0*xn[j][i].fx[s8]+xn1[j][i].fx[s8]
                   +2.0*function_f(j,i,ir,yj,xi,xn,uu)-function_f(j,i,ir,yj,xi,xn1,uu)
                   -dt2*( one6th/(xx[j][i].fx[s4]+nc*nc/xx[j][i].fx[s4])
                            *( ( ( 3.0*uu[j][i].lamdas[m]-4.0*uu[j][im].lamdas[m]
                                  +uu[j][im2].lamdas[m])/h12[i]
                                +2.0*uu[j][i].lamdas[m]*dh1[i])
                               *(3.0*xx[j][i].fx[s8]-4.0*xx[j][im].fx[s8]+xx[j][im2].fx[s8])
                              +(uu[jp][i].lamdas[m]-uu[jm][i].lamdas[m])/h22[i]
                               *(xx[jp][i].fx[s8]-xx[jm][i].fx[s8])
                              +4.0*uu[j][i].lamdas[m]
                             *(xx[jp][i].fx[s8]-2.0*xx[j][i].fx[s8]+xx[jm][i].fx[s8])/h22[i])
                         +nuv);*/


                //add heating and cooling rates
                /*if (i < ipc) {
                  if (m == 0) {
                    //heating due to exothermic reactions 
                    ff=ff-two3rd*dt2/(cdt*ns[0])*Qch[yj][ci][m];
                  }
                }*/
              }
            }
        }
    }

    return ff;
}
