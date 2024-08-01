/************************************************************************** 
 * const_normalize
 *   evaluate certain constants and perform normalization
 *
 * Jiannan Tu
 * 12/6/2013, 1/17/2014
 **************************************************************************/
#include <cmath>
#include "param.h"

using namespace std;

void const_normalize()
{
    int   i;

    pi=acos(double(-1.0));
    pi2=2.0*pi;
    deg=double(180.0)/pi;
    rad=pi/double(180.0);

/* evaluate normalization parameters from 4 basic parameters: r0, n0, B0, & mp */
    v0=B0/sqrt(mu0*n0*mp);
    t0=r0/v0;
    E0=v0*B0;
    g0=v0/t0;
    p0=B0*B0/mu0;
    j00=B0/(mu0*r0);
    e0=j00/(n0*v0);
    T0=p0/(kb*n0);
    q0=p0*v0;
    lamda0=q0*r0/T0;
    beta0=q0/j00;

    bimp=bimp/v0;
    bimpf=bimpf/v0;

    //normalize elementary charge
    e=q/e0;

    //normalize time steps
    dt=dt/t0;
    dt2=2.0*dt;

    //normalize gravitational acceleration on the surface of the earth
    gen=ge/g0;

    //normalize earth's rotational frequency
    w0n=w0*t0;

    //normalized charge/ion mass ratio for various ion species
    for (i = 0; i < sl; i++) qms[i]=e/ms[i];

    //normalized (multiplied by t0) electron gyro-frequency using normalization B0
    Omegae=q*B0*t0/(me*mp);
    
    cshap=0.5;
    
/* coefficients of collision frequencies, multiplied by n0*t0*1.0e-6 */
    coe[0]=54.5*n0*t0*1.0e-6;         //electron - ion
    coe[1]=8.90e-11*n0*t0*1.0e-6;        //electron - O
    coe[2]=1.82e-10*n0*t0*1.0e-6;        //electron - O2
    coe[3]=2.33e-11*n0*t0*1.0e-6;        //electron - N2
    coe[4]=4.50e-9*n0*t0*1.0e-6;        //electron - H

    coiO[0]=7.70e-2*n0*t0*1.0e-6;        //O+ - H+
    coiO[1]=2.60e-1*n0*t0*1.0e-6;        //O+ - O2+
    coiO[2]=2.50e-1*n0*t0*1.0e-6;        //O+ - N2+
    coiO[3]=2.60e-1*n0*t0*1.0e-6;        //O+ - NO+
    coiO[4]=1.5*3.67e-11*n0*t0*1.0e-6;   //O+ - O
    coiO[5]=6.64e-10*n0*t0*1.0e-6;       //O+ - O2
    coiO[6]=6.82e-10*n0*t0*1.0e-6;       //O+ - N2
    coiO[7]=6.61e-11*n0*t0*1.0e-6;       //O+ - H

    coiH[0]=1.23*n0*t0*1.0e-6;        //H+ - O+
    coiH[1]=1.25*n0*t0*1.0e-6;        //H+ - O2+
    coiH[2]=1.25*n0*t0*1.0e-6;        //H+ - N2+
    coiH[3]=1.25*n0*t0*1.0e-6;        //H+ - NO+
    coiH[4]=6.61e-11*n0*t0*1.0e-6;       //H+ - O
    coiH[5]=3.20e-9*n0*t0*1.0e-6;       //H+ - O2
    coiH[6]=3.36e-9*n0*t0*1.0e-6;       //H+ - N2
    coiH[7]=2.65e-10*n0*t0*1.0e-6;       //H+ - H

    coiO2[0]=1.30e-1*n0*t0*1.0e-6;       //O2+ - O+
    coiO2[1]=3.90e-2*n0*t0*1.0e-6;       //O2+ - H+
    coiO2[2]=1.50e-1*n0*t0*1.0e-6;       //O2+ - N2+
    coiO2[3]=1.60e-1*n0*t0*1.0e-6;       //O2+ - NO+
    coiO2[4]=2.31e-10*n0*t0*1.0e-6;      //O2+ - O
    coiO2[5]=2.59e-11*n0*t0*1.0e-6;      //O2+ - O2
    coiO2[6]=4.13e-10*n0*t0*1.0e-6;      //O2+ - N2
    coiO2[7]=6.50e-11*n0*t0*1.0e-6;      //O2+ - H
 
    coiN2[0]=1.50e-1*n0*t0*1.0e-6;       //N2+ - O+
    coiN2[1]=4.50e-2*n0*t0*1.0e-6;       //N2+ - H+
    coiN2[2]=1.80e-1*n0*t0*1.0e-6;       //N2+ - O2+
    coiN2[3]=1.70e-1*n0*t0*1.0e-6;       //N2+ - NO+
    coiN2[4]=2.58e-10*n0*t0*1.0e-6;      //N2+ - O
    coiN2[5]=4.49e-10*n0*t0*1.0e-6;      //N2+ - O2
    coiN2[6]=5.14e-11*n0*t0*1.0e-6;      //N2+ - N2
    coiN2[7]=7.40e-11*n0*t0*1.0e-6;      //N2+ - H
 
    coiNO[0]=1.40e-1*n0*t0*1.0e-6;       //NO+ - O+
    coiNO[1]=4.20e-2*n0*t0*1.0e-6;       //NO+ - H+
    coiNO[2]=1.70e-1*n0*t0*1.0e-6;       //NO+ - O2+
    coiNO[3]=1.60e-1*n0*t0*1.0e-6;       //NO+ - N2+
    coiNO[4]=2.44e-10*n0*t0*1.0e-6;      //NO+ - O
    coiNO[5]=4.27e-10*n0*t0*1.0e-6;      //NO+ - O2
    coiNO[6]=4.34e-10*n0*t0*1.0e-6;      //NO+ - N2
    coiNO[7]=6.90e-11*n0*t0*1.0e-6;      //NO+ - H
 
    con[0]=2.26216e-11*n0*t0*1.0e-6;    //O - O2
    con[1]=7.60616e-12*n0*t0*1.0e-6;    //O - N2
    con[2]=5.15410e-12*n0*t0*1.0e-6;    //O2 - N2
    con[3]=9.05133e-12*n0*t0*1.0e-6;    //H - O
    con[4]=5.43880e-12*n0*t0*1.0e-6;    //H - O2
    con[5]=6.05369e-12*n0*t0*1.0e-6;    //H - N2

    return;
}

