#include "param.h"

void array_init(int xs, int xm, int ys, int ym)
{
    int i, j, m;

    for (j = 0; j <= N2; j++) {
        btop[j]=-1.0e6; theta[j]=0.0; Omega1[j]=0.0; Omega2[j]=0.0;
        Omega1h[j]=0.0; Omega2h[j]=0.0; thetag[j]=0.0; phig[j]=0.0;
        zenith[j]=0.0;
    }
    for (i = 0; i <= N1; i++) {
        rr[i]=0.0; h1[i]=0.0; h2[i]=0.0; h12[i]=0.0; h22[i]=0.0; h1h[i]=0.0;
        if (i < N1) {
            gr[i]=0.0; h2h[i]=0.0; hh[i]=0.0;
        }
    }

    for (j=0; j<ym; j++) {
        for (i=0; i<xm; i++) {
            qit[j][i]=0.0;
            Qe[j][i]=0.0;
            Ce[j][i]=0.0;
            Cn[j][i]=0.0;

            for (m=0; m<5; m++) qi[j][i][m]=0.0;
            for (m=0; m< 16; m++) qib[j][i][m]=0.0;
            for (m=0; m<3; m++) Qeuv[j][i][m]=0.0;
            for (m=0; m<slm; m++) Qch[j][i][m]=0.0;
            for (m=0; m<slm+2; m++) {
                Ps[j][i][m]=0.0; Ls[j][i][m]=0.0;
            }
        }
    }
}
