 /************************************************************************
  Function: Vxtop
  Purpose: calculate time varying imposed Vx (presummedly anti-sunward
           convection) at the top boundary
 
  input : t    - time in sec
          Vimp - peak value of imposed Vx
          TL   - time to keep Vx at the value of Vimp
  output: non

  By Jiannan Tu 5/20/2011
************************************************************************/
//#include "petscsnes.h"
#include <math.h>

#include "param.h"

using namespace std;

void top_bc(double t)
{
    int    j;
    double colat, bt=0.0;

    for (j = 0; j <= N2; j++) {
        colat=theta[j]*deg;

        //Bt or Bz (normalized) increases from 0 to bimp in ts s; keeps at bimp for TL sec;
        if (t < tsl) bt=0.5*bimp*(1.0-cos(t*pi/tsl));
        else if (t >= tsl && t < TL+tsl) bt=bimp;
        //change from bimp to bimpf in ts sec
        else if (t >= TL+tsl && t <= TL+2.0*tsl) bt=bimp+0.5*(bimpf-bimp)*(1.0-cos((t-TL-tsl)*pi/tsl));
        else bt=bimp;

        if (colat<=lat[0] || colat>=lat[3]) btop[j]=bt;
        else if (colat>=lat[1] && colat<=lat[2]) btop[j]=-bt;
        else btop[j]=-1.0e5;
        btop[j]=-1.0e5;
    }

    return;
}
