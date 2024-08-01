/********************************************************************** 
 * function curvicoords.cpp
 *   Calculate scale factors of the curvilinear coordinates
 *
 * Jiannan Tu 
 * 9/13/2013, 1/10/2014, 1/28/2014
 **********************************************************************/
#include <cmath>
#include "param.h"

/*#include <iostream>
#include <fstream>
#include <iomanip>*/

using namespace std;

void curvicoords()
{
    int    i, j;
    double rN2, x1, pi2N2, rNc, thetah;

    rN2=(double)a2;
    pi2N2=2.0*pi/rN2;
    rNc=sinh(N1/d);

    //normalized radial distance of the lower & upper boundaries
    rb=rb/r0;
    ru=ru/r0;

    for(j = 0; j <= N2; j++) {
        theta[j]=pi2N2*(0.5+(double)j);

        Omega1[j]=w0n*(a13*sin(theta[j])+a33*cos(theta[j]));
        Omega2[j]=w0n*(a13*cos(theta[j])-a33*sin(theta[j]));
    }

    for (j=0; j<N2; j++) {
        thetah=0.5*(theta[j+1]+theta[j]);
        Omega1h[j]=w0n*(a13*sin(thetah)+a33*cos(thetah));
        Omega2h[j]=w0n*(a13*cos(thetah)-a33*sin(thetah));
    }

    for (i = 0; i <= N1; i++) {
        x1=double(i);
        rr[i]=rb+(ru-rb)*sinh(x1/d)/rNc;
        rh[i]=rb+(ru-rb)*sinh((x1+0.5)/d)/rNc;

        /*--- scale factor of curvilinear coordinates */
        h1[i]=(ru-rb)/(d*rNc)*cosh(x1/d);
        h2[i]=pi2N2*rr[i];
        h12[i]=h1[i]*h1[i];
        h22[i]=h2[i]*h2[i];

        //scale factors at half grids (j+1/2)
        h1h[i]=(ru-rb)/(d*rNc)*cosh((x1+0.5)/d);
        h2h[i]=pi2N2*rh[i];
        hh[i] =h1h[i]*h2h[i];
        dh1[i]=1.0/(rr[i]-h1[i])-(rr[i]-rb)/(d*d*h1[i]*h1[i]*h1[i]);

        gr[i]=gen*pow(Re/(rr[i]*r0), 2.0);

        zh[i]=(rr[i]*r0-Re);
        if (i == 0 && zh[i] < rb*r0-Re) zh[i]=rb*r0-Re;
        if (i == N1 && zh[N1] > ru*r0-Re) zh[N1]=ru*r0-Re;
    }

    for (i=0; i<=N1; i++) if (zh[i]<=2000.0e3) ipc=i;

    return;
}
