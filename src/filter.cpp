#include "param.h"
#include "funcdef.h"

#define maxim(a, b) ((a > b)?a:b)

void filter(Field **xx, int xs, int xm, int ys, int ym)
{
    int i, j, s, xi, yj;

    double *vOr, *vOt, *vOz, *vHr, *vHt, *vHz, *vO2r, *vO2t, *vO2z;
    double *vN2r, *vN2t, *vN2z, *vNOr, *vNOt, *vNOz, *dBr, *dBt, *dBz;

    s=max(xm,ym)+2;
    dBr=new double[s];
    dBt=new double[s];
    dBz=new double[s];
    vOr=new double[s];
    vOt=new double[s];
    vOz=new double[s];
    vHr=new double[s];
    vHt=new double[s];
    vHz=new double[s];
    vO2r=new double[s];
    vO2t=new double[s];
    vO2z=new double[s];
    vN2r=new double[s];
    vN2t=new double[s];
    vN2z=new double[s];
    vNOr=new double[s];
    vNOt=new double[s];
    vNOz=new double[s];

    for (j=ys; j<ys+ym; j++) {
        yj=j-ys;

        for (i=xs; i<xs+xm; i++) {
            xi=i-xs;

            dBr[xi]=xx[j][i].fx[0]*B0;
            dBt[xi]=xx[j][i].fx[1]*B0;
            dBz[xi]=xx[j][i].fx[2]*B0;
            vOr[xi]=xx[j][i].fx[5]*v0;
            vOt[xi]=xx[j][i].fx[6]*v0;
            vOz[xi]=xx[j][i].fx[7]*v0;
            vHr[xi]=xx[j][i].fx[10]*v0;
            vHt[xi]=xx[j][i].fx[11]*v0;
            vHz[xi]=xx[j][i].fx[12]*v0;
            vO2r[xi]=xx[j][i].fx[15]*v0;
            vO2t[xi]=xx[j][i].fx[16]*v0;
            vO2z[xi]=xx[j][i].fx[17]*v0;
            vN2r[xi]=xx[j][i].fx[20]*v0;
            vN2t[xi]=xx[j][i].fx[21]*v0;
            vN2z[xi]=xx[j][i].fx[22]*v0;
            vNOr[xi]=xx[j][i].fx[25]*v0;
            vNOt[xi]=xx[j][i].fx[26]*v0;
            vNOz[xi]=xx[j][i].fx[27]*v0;
        }

        shapiro(xm, dBr, 0);
        shapiro(xm, dBt, 0);
        shapiro(xm, dBz, 0);
        shapiro(xm, vOr, 0);
        shapiro(xm, vOt, 0);
        shapiro(xm, vOz, 0);
        shapiro(xm, vHr, 0);
        shapiro(xm, vHt, 0);
        shapiro(xm, vHz, 0);
        shapiro(xm, vO2r, 0);
        shapiro(xm, vO2t, 0);
        shapiro(xm, vO2z, 0);
        shapiro(xm, vN2r, 0);
        shapiro(xm, vN2t, 0);
        shapiro(xm, vN2z, 0);
        shapiro(xm, vNOr, 0);
        shapiro(xm, vNOt, 0);
        shapiro(xm, vNOz, 0);

        for (i=xs; i<xs+xm; i++) {
            xi=i-xs;

            xx[j][i].fx[0]=dBr[xi]/B0;
            xx[j][i].fx[1]=dBt[xi]/B0;
            xx[j][i].fx[2]=dBz[xi]/B0;
            xx[j][i].fx[5]=vOr[xi]/v0;
            xx[j][i].fx[6]=vOt[xi]/v0;
            xx[j][i].fx[7]=vOz[xi]/v0;
            xx[j][i].fx[10]=vHr[xi]/v0;
            xx[j][i].fx[11]=vHt[xi]/v0;
            xx[j][i].fx[12]=vHz[xi]/v0;
            xx[j][i].fx[15]=vO2r[xi]/v0;
            xx[j][i].fx[16]=vO2t[xi]/v0;
            xx[j][i].fx[17]=vO2z[xi]/v0;
            xx[j][i].fx[20]=vN2r[xi]/v0;
            xx[j][i].fx[21]=vN2t[xi]/v0;
            xx[j][i].fx[22]=vN2z[xi]/v0;
            xx[j][i].fx[25]=vNOr[xi]/v0;
            xx[j][i].fx[26]=vNOt[xi]/v0;
            xx[j][i].fx[27]=vNOz[xi]/v0;
        }
    }

    s=ym+2;
    for (i=xs; i<xs+xm; i++) {
        xi=i-xs;

        for (j=ys-1; j<=ys+ym; j++) {
            yj=j-(ys-1);

            dBr[yj]=xx[j][i].fx[0]*B0;
            dBt[yj]=xx[j][i].fx[1]*B0;
            dBz[yj]=xx[j][i].fx[2]*B0;
            vOr[yj]=xx[j][i].fx[5]*v0;
            vOt[yj]=xx[j][i].fx[6]*v0;
            vOz[yj]=xx[j][i].fx[7]*v0;
            vHr[yj]=xx[j][i].fx[10]*v0;
            vHt[yj]=xx[j][i].fx[11]*v0;
            vHz[yj]=xx[j][i].fx[12]*v0;
            vO2r[yj]=xx[j][i].fx[15]*v0;
            vO2t[yj]=xx[j][i].fx[16]*v0;
            vO2z[yj]=xx[j][i].fx[17]*v0;
            vN2r[yj]=xx[j][i].fx[20]*v0;
            vN2t[yj]=xx[j][i].fx[21]*v0;
            vN2z[yj]=xx[j][i].fx[22]*v0;
            vNOr[yj]=xx[j][i].fx[25]*v0;
            vNOt[yj]=xx[j][i].fx[26]*v0;
            vNOz[yj]=xx[j][i].fx[27]*v0;
        }

        shapiro(s, dBr, 0);
        shapiro(s, dBt, 0);
        shapiro(s, dBz, 0);
        shapiro(s, vOr,0);
        shapiro(s, vOt,0);
        shapiro(s, vOz,0);
        shapiro(s, vHr,0);
        shapiro(s, vHt,0);
        shapiro(s, vHz,0);
        shapiro(s, vO2r,0);
        shapiro(s, vO2t,0);
        shapiro(s, vO2z,0);
        shapiro(s, vN2r,0);
        shapiro(s, vN2t,0);
        shapiro(s, vN2z,0);
        shapiro(s, vNOr,0);
        shapiro(s, vNOt,0);
        shapiro(s, vNOz,0);

        for (j=ys-1; j<ys+ym; j++) {
            yj=j-(ys-1);
            if (yj==0) continue;

            xx[j][i].fx[0]=dBr[yj]/B0;
            xx[j][i].fx[1]=dBt[yj]/B0;
            xx[j][i].fx[2]=dBz[yj]/B0;
            xx[j][i].fx[5]=vOr[yj]/v0;
            xx[j][i].fx[6]=vOt[yj]/v0;
            xx[j][i].fx[7]=vOz[yj]/v0;
            xx[j][i].fx[10]=vHr[yj]/v0;
            xx[j][i].fx[11]=vHt[yj]/v0;
            xx[j][i].fx[12]=vHz[yj]/v0;
            xx[j][i].fx[15]=vO2r[yj]/v0;
            xx[j][i].fx[16]=vO2t[yj]/v0;
            xx[j][i].fx[17]=vO2z[yj]/v0;
            xx[j][i].fx[20]=vN2r[yj]/v0;
            xx[j][i].fx[21]=vN2t[yj]/v0;
            xx[j][i].fx[22]=vN2z[yj]/v0;
            xx[j][i].fx[25]=vNOr[yj]/v0;
            xx[j][i].fx[26]=vNOt[yj]/v0;
            xx[j][i].fx[27]=vNOz[yj]/v0;
        }
    }

    delete[] dBr;
    delete[] dBt;
    delete[] dBz;
    delete[] vOr;
    delete[] vOt;
    delete[] vOz;
    delete[] vHr;
    delete[] vHt;
    delete[] vHz;
    delete[] vO2r;
    delete[] vO2t;
    delete[] vO2z;
    delete[] vN2r;
    delete[] vN2t;
    delete[] vN2z;
    delete[] vNOr;
    delete[] vNOt;
    delete[] vNOz;

    return;
}
