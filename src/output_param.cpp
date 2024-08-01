/* output_param
 *   Output various system constants & parameters, including background magnetic
 *   field. All vectors are in curvilinear coordinates. 
 *
 *   Jiannan Tu
 *   12/9/2013, 1/29/2014
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>

using namespace std;

#include "param.h"

int output_param(ofstream &outfstr,char *workdir,int nprev)
{
    int      i, j;
    char     fname[150];
    ofstream facfstr;

/* basic parameters*/
    outfstr << a1 << " " << a2 << "    a1 & a2: # of grids along x1 & x2" << endl;
    outfstr << iyr << " " << mon << " " << idate << " "
            << sec << " " << LT << "    iyr, mon, idate, sec & LT" << endl;
    outfstr << Ap << " " << f107 << " " << f107a
            << "    Ap, f107 & f107a" << endl;
    outfstr << r0 << " " << n0 << " " << B0 << " " << v0 << " "
            << E0 << "   r0, n0, B0, v0, E0" <<endl;
    outfstr << t0 << " " << g0 << " " << p0 << " " << j00 << " "
            << e0 << "   t0, g0, p0, j0, e0" <<endl;
    outfstr << T0 << " " << q0 << " " << lamda0 << " " << beta0 << " "
            << "   T0, q0, lamda0, beta0" <<endl;
    outfstr << dt*t0 << " " << (rb*r0-Re)*1.0e-3 << " " << (ru*r0-Re)*1.0e-3
            << "   dt (s), zmin & zmax (km)" << endl;
    outfstr << ntot << " " << nout
            << "   number of time step to run & number of time step between outputs" << endl;
    outfstr << sl << " " << sm << "   sl & sm # of ion and neutral species" << endl;
    outfstr << dth << "   co-latitude width " << endl;
    outfstr << d << "  radial grids stretching factor" << endl;
    outfstr << lat[0] << " " << lat[1]
            << " co-latitudes (deg) of the perturbation magnetic field reversal (north)" << endl;
    outfstr << lat[2] << " " << lat[3]
            << " co-latitudes (deg) of the perturbation magnetic field reversal (south)" << endl;
    outfstr << bimp*v0 << " imposed velocity (m/s) at the top boundary" << endl;
    outfstr << bimpf*v0 << " new value of bimp (m/s) after time period of TL s" << endl;
    outfstr << tsl << " time (s) to increase imposed velocity at the top bndry from 0 to bimp"
                        << endl;
    outfstr << TL << " time period (s) of keeping the imposed velocity" << endl;
    outfstr << nvt << " if = 0: imposing v_theta; else = 1: imposing v_z" << endl;
    if (smod == 1) outfstr << nprev << " time step at which continuous run starts" << endl;

/* magnetic longitude in deg */
    outfstr<<fixed<<setw(7)<<setprecision(2)<<phi[0]*deg<<" Magnetic longitude (deg)"<<endl;

    strncpy(fname, workdir, 150);
    strcat(fname, "/output/alt_corolis.dat");
    facfstr.open(fname, fstream::out);
    if(!facfstr) {
        cout << "Can't open file " << fname << endl;
        return -1;
    }

/* normalized radial distance, altitude (km) & delta_r in km */
    facfstr << setw(15) << "Normalized    r" << setw(28) << "Altitude (km)"
            << setw(19) << "delta_r (km)" << endl; 
    for (i = 0; i <= N1; i++) {
        facfstr << scientific << setw(23) << setprecision(16) << rr[i];
        facfstr << scientific << setw(24) << setprecision(16) << zh[i]/1.0e3;
        if (i > 0) facfstr << scientific << setw(16) << setprecision(8)
                           << (zh[i]-zh[i-1])/1.0e03; 
        facfstr << endl;
    }

/*  normalized angular velocity of the earth: Omega1 & Omega2 at (i,j)
 *  in curvilinear coordinates */
    facfstr << setw(15) << "Omega1" << setw(23) << "Omega2" << endl;
    for (j = 0; j <= N2; j++) {
        facfstr << scientific << setw(23) << setprecision(16) << Omega1[j]
                << scientific << setw(24) << setprecision(16) << Omega2[j];
        facfstr << endl;
    }
    facfstr.close();

    strncpy(fname, workdir, 150);
    strcat(fname, "/output/lat_lon.dat");
    facfstr.open(fname, fstream::out);
    if(!facfstr) {
        cout << "Can't open file " << fname << endl;
        return -1;
    }

/* geographic colatitude and longitude (rad) of grids in the magnetic meridian */
    facfstr<<"UT in hour "<<fixed<<setw(5)<<setprecision(2)<<sec/3600.0<< endl;
    facfstr<<setw(17)<<"Geo Colat"<<setw(23)<<"Geo long"<<setw(25)<<"Mag theta" 
           <<setw(14)<<"GEO LT"<<endl;
    for (j = 0; j <= N2; j++) {
        LT=sec/3600.0+phig[j]*deg/15.0;
        if (LT >= 24.0) LT=LT-24.0;
        facfstr << scientific << setw(24) << setprecision(16) << thetag[j]
                << scientific << setw(24) << setprecision(16) << phig[j]
                << scientific << setw(24) << setprecision(16) << theta[j]
                << fixed << setw(7) << setprecision(2) << LT << endl;
    }
    facfstr.close();

    strncpy(fname, workdir, 150);
    strcat(fname, "/output/factors_h.dat");
    facfstr.open(fname, fstream::out);
    if(!facfstr) {
        cout << "Can't open file " << fname << endl;
        return -1;
    }

/* scale factors h1, h2, & h, as well as dh1/dx2, 2/h*d(h1/h2)/dx2 */
    facfstr << setw(14) << "h1" << setw(24) << "h2" << setw(23)
            << "h2h" << setw(24) << "h1h" << setw(24) << "hh" << endl;
    for (i = 0; i <= N1; i++) {
        facfstr << scientific << setw(23) << setprecision(16) << h1[i]
                << scientific << setw(24) << setprecision(16) << h2[i]
                << scientific << setw(24) << setprecision(16) << h2h[i]
                << scientific << setw(24) << setprecision(16) << h1h[i]
                << scientific << setw(24) << setprecision(16) << hh[i];
        facfstr << endl;
    }

    facfstr.close();

    return 0;
}
