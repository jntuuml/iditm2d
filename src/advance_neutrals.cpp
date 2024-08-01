#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "param.h"
#include "funcdef.h"

PetscErrorCode advance_neutrals(DM dau,Vec U, Field **xx,double sdt)
{
    Vec    localU;
    Fieldu **uu;
    int    i, j, l, m, ll, msl, lsl, ierr;
    int    s0, s, s1, s2, s3, s4;
    double ns[5], rhon, nuv=0.0;
    int    jp, ip, ip2, jm, im, im2, c0, xs, xm, ys, ym, yj, xi;
    double dNr, dvr, dFr, dlamr, dTr, ddTr, vd[3], Td, vd2;

    const double one3rd=0.3333333333333333;
    const double two3rd=0.6666666666666667;
    const double one6th=0.1666666666666667;

    ierr = DMGetLocalVector(dau,&localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dau,U,INSERT_VALUES,localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dau,U,INSERT_VALUES,localU);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dau,localU,&uu);CHKERRQ(ierr);

    ierr=DMDAGetCorners(dau,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

    for (j=ys; j<ys+ym; j++) {
      jp=j+1; jm=j-1; yj=j-ys;

      for (i=xs; i<xs+xm; i++) {
        ip=i+1; ip2=i+2; im=i-1; im2=i-2;

        xi=i-xs;

        for (m = 0; m < sm; m++) {
          s=5*m; s1=s+1; s2=s+2; s3=s+3; s4=s+4;
	  msl=m+sl; c0=9*(msl+1);

          if (zh[i] >= 1500.0e3) {
            uu[j][i].fu[s] =h2[im]/h2[i]*uu[j][im].fu[s];
            uu[j][i].fu[s1]=h2[im]/h2[i]*uu[j][im].fu[s1];
            uu[j][i].fu[s2]=h2[im]/h2[i]*uu[j][im].fu[s2];
            uu[j][i].fu[s3]=h2[im]/h2[i]*uu[j][im].fu[s3];
            uu[j][i].fu[s4]=h2[i]/h2[im]*uu[j][im].fu[s4];
          }
          else {
// ---- difference equations of Nn */
            if (i==0) {
              dNr=0.5*(3.0*uu[j][ip2].fu[s]-4.0*uu[j][ip].fu[s]+uu[j][i].fu[s]);
              dvr=0.5*(3.0*uu[j][ip2].fu[s1]-4.0*uu[j][ip].fu[s1]+uu[j][i].fu[s1]);
            }
            else if(i>0 && i<N1) {
              dNr=0.5*(uu[j][ip].fu[s]-uu[j][im].fu[s]);
              dvr=0.5*(uu[j][ip].fu[s1]-uu[j][im].fu[s1]);
            }
            else {
              dNr=0.5*(3.0*uu[j][i].fu[s]-4.0*uu[j][im].fu[s]+uu[j][im2].fu[s]);
              dvr=0.5*(3.0*uu[j][i].fu[s1]-4.0*uu[j][im].fu[s1]+uu[j][im2].fu[s1]);
            }

            uu[j][i].fu[s]= uu[j][i].fu[s]
                 -sdt*( (dvr+uu[j][i].fu[s1]*dNr)/h1[i]+uu[j][i].fu[s1]/rr[i]
                       +0.5*( (uu[jp][i].fu[s2]-uu[jm][i].fu[s2])+
                             +uu[j][i].fu[s2]*(uu[jp][i].fu[s]-uu[jm][i].fu[s]))/h2[i]);
            if (i<=ipc) {
                cout<<j<<", "<<i<<", "<<m<<", "<<Ps[yj][xi][msl]<<", "<<Ls[yj][xi][msl]<<", "<<uu[j][i].nu[s]<<endl;
                uu[j][i].fu[s]=uu[j][i].fu[s]
                           +sdt*(Ps[yj][xi][msl]/exp(uu[j][i].fu[s])-Ls[yj][xi][msl]);
            }

// ---- difference equations of vs1
            if (i==0) uu[j][i].fu[s1]=h2[i]/h2[ip]*uu[j][ip].fu[s1];
            else if (i>0 && i<N1) {
              //collision terms in ion & neutral momentum equations
              nuv= uu[j][i].nu[c0]*(uu[j][i].fu[s1]-xx[j][i].fx[3]);

              for (l = 1; l < slm; l++) {
                if (l <= sl) {
                  nuv=nuv+uu[j][i].nu[c0+l]*(uu[j][i].fu[s1]-xx[j][i].fx[5+5*(l-1)]);
                }
                else {
                  lsl=l-sl;
                  if (lsl <= m) ll=lsl-1; else ll=lsl;
                  nuv=nuv+uu[j][i].nu[c0+l]*(uu[j][i].fu[s1]-uu[j][i].fu[1+5*ll]);
                }
              }

              uu[j][i].fu[s1]=uu[j][i].fu[s1]
                   -sdt*( 0.5*uu[j][i].fu[s1]*(uu[j][ip].fu[s1]-uu[j][im].fu[s1])/h1[i]
                         -uu[j][i].fu[s2]*uu[j][i].fu[s2]/rr[i]
                         +0.5*uu[j][i].fu[s2]*(uu[jp][i].fu[s1]-uu[jm][i].fu[s1])/h2[i]
                         +0.5/(ms[msl]*h1[i])
                             *( uu[j][i].fu[s4]*(uu[j][ip].fu[s]-uu[j][im].fu[s])
                               +(uu[j][ip].fu[s4]-uu[j][im].fu[s4]))
                         +gr[i]+Omega2[j]*(2.0*uu[j][i].fu[s3]-Omega2[j]*rr[i])+nuv);
            }
            else uu[j][i].fu[s1]=h2[im]/h2[i]*uu[j][im].fu[s1];

// ---- difference equations of vs2
            if (i==0) uu[j][i].fu[s2]=h2[i]/h2[ip]*uu[j][ip].fu[s2];
            else {
              nuv=uu[j][i].nu[c0]*(uu[j][i].fu[s2]-xx[j][i].fx[4]);

              for (l = 1; l < slm; l++) {
                if (l <= sl) {
                  nuv=nuv+uu[j][i].nu[c0+l]*(uu[j][i].fu[s2]-xx[j][i].fx[6+5*(l-1)]);
                }
                else {
                  lsl=l-sl;
                  if (lsl <= m) ll=lsl-1; else ll=lsl;
                  nuv=nuv+uu[j][i].nu[c0+l]*(uu[j][i].fu[s2]-uu[j][i].fu[2+5*ll]);
                }
              }

              if (i<N1) dvr=0.5*(uu[j][ip].fu[s2]-uu[j][im].fu[s2]);
              else dvr=0.5*(3.0*uu[j][i].fu[s2]-4.0*uu[j][im].fu[s2]+uu[j][im2].fu[s2]);

              uu[j][i].fu[s2]=uu[j][i].fu[s2]
                   -sdt*( uu[j][i].fu[s1]*dvr/h1[i]+uu[j][i].fu[s1]*uu[j][i].fu[s2]/rr[i]
                         +0.5*uu[j][i].fu[s2]*(uu[jp][i].fu[s2]-uu[jm][i].fu[s2])/h2[i]
                         +0.5/(ms[msl]*h2[i])
                             *( uu[j][i].fu[s4]*(uu[jp][i].fu[s]-uu[jm][i].fu[s])
                               +(uu[jp][i].fu[s4]-uu[jm][i].fu[s4]))
                         +Omega1[j]*(Omega2[j]*rr[i]-2.0*uu[j][i].fu[s3])+nuv);

            }

// ---- difference equations for vs3. vs3,i,0= 0
            if (i==0) uu[j][i].fu[s3]=h2[i]/h2[ip]*uu[j][ip].fu[s3];
            else {
              nuv=uu[j][i].nu[c0]*(uu[j][i].fu[s3]-xx[j][i].fx[5]);

              ll=sl+1;
              for (l = 1; l < slm; l++) {
                if (l <= sl) {
                  nuv=nuv+uu[j][i].nu[c0+l]*(uu[j][i].fu[s3]-xx[j][i].fx[7+5*(l-1)]);
                }
                else {
                  lsl=l-sl;
                  if (lsl <= m) ll=lsl-1; else ll=lsl;
                  nuv=nuv+uu[j][i].nu[c0+l]*(uu[j][i].fu[s3]-uu[j][i].fu[3+5*ll]);
                }
              }

              if (i<N1) dvr=0.5*(uu[j][ip].fu[s3]-uu[j][im].fu[s3]);
              else dvr=0.5*(3.0*uu[j][i].fu[s3]-4.0*uu[j][im].fu[s3]+uu[j][im2].fu[s3]);

              uu[j][i].fu[s3]=uu[j][i].fu[s3]
                   -sdt*( uu[j][i].fu[s1]*dvr/h1[i]
                         +0.5*uu[j][i].fu[s2]*(uu[jp][i].fu[s3]-uu[jm][i].fu[s3])/h2[i]
                         +2.0*(Omega1[j]*uu[j][i].fu[s2]-Omega2[j]*uu[j][i].fu[s1])
			 +nuv);
            }

// ---- difference equations for Ts. first calculate rates of heating (or cooling) by
//      electron, neutral & other ion species */
            if (i > 0) {
              //first calculate rates of heating (or cooling) by electron, neutral
              // & other ion species
              vd[0]=xx[j][i].fx[3]-uu[j][i].fu[s1];
              vd[1]=xx[j][i].fx[4]-uu[j][i].fu[s2];
              vd[2]=xx[j][i].fx[5]-uu[j][i].fu[s3];
              nuv= uu[j][i].nu[c0]*( 2.0*(xx[j][i].fx[3]-uu[j][i].fu[s4])
                                    +two3rd*me*(vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2]));

              for (l = 1; l < slm; l++) {
                if (l <= sl) {
                  s0=5*(l-1);
                  vd[0]=xx[j][i].fx[5+s0]-uu[j][i].fu[s1];
                  vd[1]=xx[j][i].fx[6+s0]-uu[j][i].fu[s2];
                  vd[2]=xx[j][i].fx[7+s0]-uu[j][i].fu[s3];
                  vd2  =vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2];

                  Td=uu[j][i].nu[c0+l]*ms[msl]/(ms[msl]+ms[l-1]);
                  nuv=nuv+Td*(2.0*(xx[j][i].fx[8+s0]-uu[j][i].fu[s4])+two3rd*ms[l-1]*vd2);
                }
                else {
                  lsl=l-sl;
                  if (lsl <= m) ll=lsl-1; else ll=lsl;
                  vd[0]=uu[j][i].fu[1+5*ll]-uu[j][i].fu[s1];
                  vd[1]=uu[j][i].fu[2+5*ll]-uu[j][i].fu[s2];
                  vd[2]=uu[j][i].fu[3+5*ll]-uu[j][i].fu[s3];
                  vd2  =vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2];

                  Td=uu[j][i].nu[c0+l]*ms[msl]/(ms[msl]+ms[ll+sl]);
 
                  nuv=nuv+Td*(2.0*(uu[j][i].fu[4+5*ll]-uu[j][i].fu[s4])+two3rd*ms[ll]*vd2);
                }
              }

              ns[0]=exp(uu[j][i].fu[s]);

              if (i < N1) {
                dFr=0.5*( uu[j][i].fu[s1]*(uu[j][ip].fu[s4]-uu[j][im].fu[s4])
                         +uu[j][i].fu[s4]*(uu[j][ip].fu[s1]-uu[j][im].fu[s1]));
                dvr=0.5*(uu[j][ip].fu[s1]-uu[j][im].fu[s1]);
                dlamr=(uu[j][ip].lamdas[msl]-uu[j][im].lamdas[msl]);
                dTr=(uu[j][ip].fu[s4]-uu[j][im].fu[s4]);
                ddTr=(uu[j][ip].fu[s4]-2.0*uu[j][i].fu[s4]+uu[j][im].fu[s4]);
              }
              else {
                dFr=0.5*( 3.0*uu[j][i].fu[s4]*uu[j][i].fu[s1]
                         -4.0*uu[j][im].fu[s4]*uu[j][im].fu[s1]
                         +uu[j][im2].fu[s4]*uu[j][im2].fu[s1]);
                dvr=0.5*(3.0*uu[j][i].fu[s1]-4.0*uu[j][im].fu[s1]+uu[j][im2].fu[s1]);
                dlamr=( 3.0*uu[j][i].lamdas[msl]-4.0*uu[j][im].lamdas[msl]
                       +uu[j][im2].lamdas[msl]);
                dTr=(3.0*uu[j][i].fu[s4]-4.0*uu[j][im].fu[s4]+uu[j][im2].fu[s4]);
                ddTr=0.0;
              }

              uu[j][i].fu[s4]=uu[j][i].fu[s4]
                   +sdt*( one3rd*uu[j][i].fu[s4]
                                *( dvr/h1[i]+uu[j][i].fu[s1]/rr[i]
                                  +0.5*(uu[jp][i].fu[s2]-uu[jm][i].fu[s2])/h2[i])
                         -dFr/h1[i]-uu[j][i].fu[s4]*uu[j][i].fu[s1]/rr[i]
                         -0.5*( uu[j][i].fu[s2]*(uu[jp][i].fu[s4]-uu[jm][i].fu[s4])
                               +uu[j][i].fu[s4]*(uu[jp][i].fu[s2]-uu[jm][i].fu[s2]))/h2[i]
                         +one6th/ns[0]
                                *( (dlamr/h12[i]+2.0*uu[j][i].lamdas[msl]*dh1[i])*dTr
                                  +(uu[jp][i].lamdas[msl]-uu[jm][i].lamdas[msl])/h22[i]
                                   *(uu[jp][i].fu[s4]-uu[jm][i].fu[s4])
                                  +4.0*uu[j][i].lamdas[msl]
                                      *(ddTr/h12[i]+( uu[jp][i].fu[s4]-2.0*uu[j][i].fu[s4]
                                                     +uu[jm][i].fu[s4])/h22[i]))
                         +nuv);

              if(i<=ipc && m<sm-1) {
                rhon=0.0;
                for (l=0; l<sm; l++) rhon=rhon+exp(uu[j][i].fu[5*l])*ms[l+sl];

                uu[j][i].fu[s4]=uu[j][i].fu[s4]
		              +two3rd*(Qeuv[yj][xi][m]/ns[0]-ms[msl]/rhon*Cn[yj][xi]);
              }
            }
          }
        }

        if (i<=ipc) {
          uu[j][i].fu[20]= uu[j][i].fu[20]
                          +sdt*(Ps[yj][xi][9]/exp(uu[j][i].fu[20])-Ls[yj][xi][9]);
          uu[j][i].fu[21]= uu[j][i].fu[21]
                          +sdt*(Ps[yj][xi][10]/exp(uu[j][i].fu[21])-Ls[yj][xi][10]);
        }
        else {
          uu[j][i].fu[20]=uu[j][im].fu[20];
          uu[j][i].fu[21]=uu[j][im].fu[21];
        }
      }
    }

    ierr = DMDAVecRestoreArray(dau,localU,&uu);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dau,&localU);CHKERRQ(ierr);

    return 0;
}
