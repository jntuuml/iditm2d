/************************************************************************
! Routine: collision_freq
! Purpose: calculate electron, ion, and neutral frequencies
!
! input : non
! output: non
!
! By   Jiannan Tu
! Date 12/10/2013, 12/31/2013, 1/17/2014
!*************************************************************************/
#include <cmath>
#include "param.h"
#include "funcdef.h"

using namespace std;

void collision_freq(Field **xx,Fieldu **uu,int xs,int xm,int ys,int ym)
{
    double Te, Te32, Te12, ne, Tr, Trn;
    double niO, niO2, niN2, niH, niNO;
    double nnO, nnO2, nnN2, nnH;
    double TiO, TiO2, TiN2, TiH, TiNO;
    double TnO, TnO2, TnN2, TnH;
    int    i, j;

    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            Te=xx[j][i].fx[3]*T0;
            Te32=pow(Te, 1.5);
            Te12=sqrt(Te);

            //Densities (cm^-3) & temperatures of neutral species
            niO =xx[j][i].fx[4];
            niH =xx[j][i].fx[9];
            niO2=xx[j][i].fx[14];
            niN2=xx[j][i].fx[19];
            niNO=xx[j][i].fx[24];
            ne=niO+niH+niO2+niN2+niNO;

            nnO =exp(uu[j][i].fu[0]);
            nnO2=exp(uu[j][i].fu[5]);
            nnN2=exp(uu[j][i].fu[10]);
            nnH =exp(uu[j][i].fu[15]);

            TiO =xx[j][i].fx[8]*T0;
            TiH =xx[j][i].fx[13]*T0;
            TiO2=xx[j][i].fx[18]*T0;
            TiN2=xx[j][i].fx[23]*T0;
            TiNO=xx[j][i].fx[28]*T0;

            TnO =uu[j][i].fu[4]*T0;
            TnO2=uu[j][i].fu[9]*T0;
            TnN2=uu[j][i].fu[14]*T0;
            TnH =uu[j][i].fu[19]*T0;

/* ---- electron collision frequencies: with O+, O2+, N2+, H+ ------------------------------*/
            // electron Coulomb collisions
            uu[j][i].nu[0]=coe[0]*niO/Te32;     // e - O+
            uu[j][i].nu[1]=coe[0]*niH/Te32;     // e - H+
            uu[j][i].nu[2]=coe[0]*niO2/Te32;    // e - O2+
            uu[j][i].nu[3]=coe[0]*niN2/Te32;    // e - N2+
            uu[j][i].nu[4]=coe[0]*niNO/Te32;    // e - NO+

            /* electron - neutral collision frequencies */
            uu[j][i].nu[5]=coe[1]*nnO *(1.0+5.7e-4*Te)*Te12;        // e - O
            uu[j][i].nu[6]=coe[2]*nnO2*(1.0+3.6e-2*Te12)*Te12;      // e - O2
            uu[j][i].nu[7]=coe[3]*nnN2*(1.0-1.21e-4*Te)*Te;         // e - N2
            uu[j][i].nu[8]=coe[4]*nnH *(1.0-1.35e-4*Te)*Te12;       // e - H

/* ---- H - neutral collision frequency calculate first since will be used below -----------*/
            Trn=0.5*(TnH+TnO);
            uu[j][i].nu[87]=con[3]*(nnH+16.0*nnO)*pow(Trn, 0.292);    // H - O

            Trn=0.5*(TnH+TnO2);
            uu[j][i].nu[88]=con[4]*(nnH+32.0*nnO2)*pow(Trn, 0.289);   // H - O2

            Trn=0.5*(TnH+TnN2);
            uu[j][i].nu[89]=con[5]*(nnH+28.0*nnN2)*pow(Trn, 0.302);   // H - N2

/* ---- O+ collision frequencies -----------------------------------------------------------*/
            /* O+ Coulomb collision frequencies */
            uu[j][i].nu[9]=ne*me/(niO*ms[0])*uu[j][i].nu[0];           // O+ - e
            uu[j][i].nu[10]=coiO[0]*niH/pow(TiH, 1.5);                 // O+ - H+
            uu[j][i].nu[11]=coiO[1]*niO2/pow(TiO2, 1.5);               // O+ - O2+
            uu[j][i].nu[12]=coiO[2]*niN2/pow(TiN2, 1.5);               // O+ - N2+
            uu[j][i].nu[13]=coiO[3]*niNO/pow(TiNO, 1.5);               // O+ - NO+

            /* O+ - neutral collision frequencies */
            Tr=0.5*(TiO+TnO);
            uu[j][i].nu[14]=coiO[4]*nnO*sqrt(Tr)*pow(1.0-0.064*log10(Tr), 2.0);   // O+ - O
            uu[j][i].nu[15]=coiO[5]*nnO2;                                         // O+ - O2
            uu[j][i].nu[16]=coiO[6]*nnN2;                                         // O+ - N2
            uu[j][i].nu[17]=coiO[7]*nnH*sqrt(TiO)*pow(1.0-0.047*log10(TiO), 2.0); // O+ - H

/* ---- H+ collision frequencies -----------------------------------------------------------*/
            /* H+ Coulomb collision frequencies */
            uu[j][i].nu[18]=ne*me/(niH*ms[1])*uu[j][i].nu[1];          // H+ - e
            uu[j][i].nu[19]=coiH[0]*niO/pow(TiO, 1.5);                 // H+ - O+
            uu[j][i].nu[20]=coiH[1]*niO2/pow(TiO2, 1.5);               // H+ - O2+
            uu[j][i].nu[21]=coiH[2]*niN2/pow(TiN2, 1.5);               // H+ - N2+
            uu[j][i].nu[22]=coiH[3]*niNO/pow(TiNO, 1.5);               // H+ - NO+

           /* H+ - neutral collision frequencies */
            uu[j][i].nu[23]=coiH[4]*nnO*sqrt(TiH)*pow(1.0-0.047*log10(TiH), 2.0); // H+ - O
            uu[j][i].nu[24]=coiH[5]*nnO2;                                         // H+ - O2
            uu[j][i].nu[25]=coiH[6]*nnN2;                                         // H+ - N2

            Tr=0.5*(TiH+TnH);
            uu[j][i].nu[26]=coiH[7]*nnH*sqrt(Tr)*pow(1.0-0.083*log10(Tr), 2.0);   // H+ - H

/* ---- O2+ collision frequencies ----------------------------------------------------------*/
            /* O2+ Coulomb collision frequencies */
            uu[j][i].nu[27]=ne*me/(niO2*ms[2])*uu[j][i].nu[2];          // O2+ - e
            uu[j][i].nu[28]=coiO2[0]*niO/pow(TiO, 1.5);                 // O2+ - O+
            uu[j][i].nu[29]=coiO2[1]*niH/pow(TiH, 1.5);                 // O2+ - H+
            uu[j][i].nu[30]=coiO2[2]*niN2/pow(TiN2, 1.5);               // O2+ - N2+
            uu[j][i].nu[31]=coiO2[3]*niNO/pow(TiNO, 1.5);               // O2+ - NO+

            /* O2+ - neutral collision frequencies */
            uu[j][i].nu[32]=coiO2[4]*nnO;                                         // O2+ - O

            Tr=0.5*(TiO2+TnO2);
            uu[j][i].nu[33]=coiO2[5]*nnO2*sqrt(Tr)*pow(1.0-0.073*log10(Tr), 2.0); // O2+ - O2
            uu[j][i].nu[34]=coiO2[6]*nnN2;                                        // O2+ - N2
            uu[j][i].nu[35]=coiO2[7]*nnH;                                         // O2+ - H

/* ---- N2+ collision frequencies ----------------------------------------------------------*/
            /* N2+ Coulomb collision frequencies */
            uu[j][i].nu[36]=ne*me/(niN2*ms[3])*uu[j][i].nu[3];           // N2+ - e
            uu[j][i].nu[37]=coiN2[0]*niO/pow(TiO, 1.5);                  // N2+ - O+
            uu[j][i].nu[38]=coiN2[1]*niH/pow(TiH, 1.5);                  // N2+ - H+
            uu[j][i].nu[39]=coiN2[2]*niO2/pow(TiO2, 1.5);                // N2+ - N2+
            uu[j][i].nu[40]=coiN2[3]*niNO/pow(TiNO, 1.5);                // N2+ - NO+

            /* N2+ - neutral collision frequencies */
            uu[j][i].nu[41]=coiN2[4]*nnO;                                // N2+ - O
            uu[j][i].nu[42]=coiN2[5]*nnO2;                               // N2+ - O2

            Tr=0.5*(TiN2+TnN2);
            uu[j][i].nu[43]=coiN2[6]*nnN2*sqrt(Tr)*pow(1.0-0.069*log10(Tr), 2.0); // N2+ - N2
            uu[j][i].nu[44]=coiN2[7]*nnH;                                         // N2+ - H

/* ---- NO+ collision frequencies ----------------------------------------------------------*/
           /* NO+ Coulomb collision frequencies */
            uu[j][i].nu[45]=ne*me/(niNO*ms[4])*uu[j][i].nu[4];      // NO+ - e
            uu[j][i].nu[46]=coiNO[0]*niO/pow(TiO, 1.5);             // NO+ - O+
            uu[j][i].nu[47]=coiNO[1]*niH/pow(TiH, 1.5);             // NO+ - H+
            uu[j][i].nu[48]=coiNO[2]*niO2/pow(TiO2, 1.5);           // NO+ - O2+
            uu[j][i].nu[49]=coiNO[3]*niN2/pow(TiN2, 1.5);           // NO+ - N2+

           /* NO+ - neutral collision frequencies */
            uu[j][i].nu[50]=coiNO[4]*nnO;                            // NO+ - O
            uu[j][i].nu[51]=coiNO[5]*nnO2;                           // NO+ - O2
            uu[j][i].nu[52]=coiNO[6]*nnN2;                           // NO+ - N2
            uu[j][i].nu[53]=coiNO[7]*nnH;                            // NO+ - H

/* ---- O collision frequencies ------------------------------------------------------------*/
            /* O - e, O - ion collision frequencies */
            uu[j][i].nu[54]=ne*me/(nnO*ms[sl])*uu[j][i].nu[5];         // O - e
            uu[j][i].nu[55]=niO/nnO*uu[j][i].nu[14];                   // O - O+
            uu[j][i].nu[56]=niH*ms[1]/(nnO*ms[sl])*uu[j][i].nu[23];    // O - H+
            uu[j][i].nu[57]=niO2*ms[2]/(nnO*ms[sl])*uu[j][i].nu[32];   // O - O2+
            uu[j][i].nu[58]=niN2*ms[3]/(nnO*ms[sl])*uu[j][i].nu[41];   // O - N2+
            uu[j][i].nu[59]=niNO*ms[4]/(nnO*ms[sl])*uu[j][i].nu[50];   // O - NO+

            /* O - neutral collision frequencies */
            Trn=0.5*(TnO+TnO2);
            uu[j][i].nu[60]=con[0]*(nnO+2.0*nnO2)*pow(Trn,0.226);      // O - O2

            Trn=0.5*(TnO+TnN2);
            uu[j][i].nu[61]=con[1]*(4.0*nnO+7.0*nnN2)*pow(Trn,0.226);     // O - N2
            uu[j][i].nu[62]=nnH*ms[8]/(nnO*ms[sl])*uu[j][i].nu[87];    // O - H
            if (i>0 && uu[j][i].nu[62]>uu[j][i-1].nu[62]) uu[j][i].nu[62]=uu[j][i-1].nu[62];

/* ---- O2 collision frequencies -----------------------------------------------------------*/
           /* O2 - e, O2 - ion collision frequencies  */
            uu[j][i].nu[63]=ne*me/(nnO2*ms[6])*uu[j][i].nu[6];         // O2 - e
            uu[j][i].nu[64]=niO*ms[0]/(nnO2*ms[6])*uu[j][i].nu[15];    // O2 - O+
            uu[j][i].nu[65]=niH*ms[1]/(nnO2*ms[6])*uu[j][i].nu[24];    // O2 - H+
            uu[j][i].nu[66]=niO2/nnO*uu[j][i].nu[33];                  // O2 - O2+
            uu[j][i].nu[67]=niN2*ms[3]/(nnO2*ms[6])*uu[j][i].nu[42];   // O2 - N2+
            uu[j][i].nu[68]=niNO*ms[4]/(nnO2*ms[6])*uu[j][i].nu[51];   // O2 - NO+

            /* O2 - neutrals collision frequencies */
            uu[j][i].nu[69]=nnO*ms[sl]/(nnO2*ms[6])*uu[j][i].nu[60];   // O2 - O

            Trn=0.5*(TnO2+TnN2);
            uu[j][i].nu[70]=con[2]*(8.0*nnO2+7.0*nnN2)*pow(Trn, 0.25);    // O2 - N2
            uu[j][i].nu[71]=nnH*ms[8]/(nnO2*ms[6])*uu[j][i].nu[88];    // O2 - H
            if (i>0 && uu[j][i].nu[71]>uu[j][i-1].nu[71]) uu[j][i].nu[71]=uu[j][i-1].nu[71];

/* ---- N2 collision frequencies -----------------------------------------------------------*/
            /* N2 - e, N2 -ion collision frequencies */
            uu[j][i].nu[72]=ne*me/(nnN2*ms[7])*uu[j][i].nu[7];         // N2 - e
            uu[j][i].nu[73]=niO*ms[0]/(nnN2*ms[7])*uu[j][i].nu[16];    // N2 - O+
            uu[j][i].nu[74]=niH*ms[1]/(nnN2*ms[7])*uu[j][i].nu[25];    // N2 - H+
            uu[j][i].nu[75]=niO2*ms[2]/(nnN2*ms[7])*uu[j][i].nu[34];   // N2 - O2+
            uu[j][i].nu[76]=niN2/nnN2*uu[j][i].nu[43];                 // N2 - N2+
            uu[j][i].nu[77]=niNO*ms[4]/(nnN2*ms[7])*uu[j][i].nu[52];   // N2 - NO+

            /* N2 - neutral collision frequencies */
            uu[j][i].nu[78]=nnO*ms[5]/(nnN2*ms[7])*uu[j][i].nu[61];    // N2 - O
            uu[j][i].nu[79]=nnO2*ms[6]/(nnN2*ms[7])*uu[j][i].nu[70];   // N2 - O2
            uu[j][i].nu[80]=nnH*ms[8]/(nnN2*ms[7])*uu[j][i].nu[89]; // N2 - H
            if (i>0 && uu[j][i].nu[80]>uu[j][i-1].nu[80]) uu[j][i].nu[80]=uu[j][i-1].nu[80];

/*----- H collision frequencies ------------------------------------------------------------*/
            // H - e, H - ion collision frequencies
            uu[j][i].nu[81]=ne*me/(nnH*ms[8])*uu[j][i].nu[8];          // H - e
            uu[j][i].nu[82]=niO*ms[0]/(nnH*ms[8])*uu[j][i].nu[17];     // H - O+
            uu[j][i].nu[83]=niH*ms[1]/nnH*uu[j][i].nu[26];             // H - H+
            uu[j][i].nu[84]=niO2*ms[2]/(nnH*ms[8])*uu[j][i].nu[35];    // H - O2+
            uu[j][i].nu[85]=niN2*ms[3]/(nnH*ms[8])*uu[j][i].nu[44];    // H - N2+
            uu[j][i].nu[86]=niNO*ms[4]/(nnH*ms[8])*uu[j][i].nu[53];    // H - NO+
        }
    }

    return;
}
