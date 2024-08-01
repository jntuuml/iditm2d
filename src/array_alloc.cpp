/*******************************************************************************
 *  array_alloc
 *  allocate memory for multi-dimension arrays
 * 
 *  Jiannan Tu
 *  9/15/2013, 1/29/2014
 ******************************************************************************/
#include <cmath>
using namespace std;

#include "param.h"
#include "funcdef.h"

void array_alloc(int xs,int xm,int ys,int ym)
{
    int      i,j,xi,yj;

/* top BC values of V */
    btop = new double[a2];

//spherical coordinates
    rr = new double[a1];
    rh = new double[a1];
    theta=new double[a2];
    phi =new double[1];
    zh =new double[a1];

// ---- scale factors of curvilinear coordinates (at both j and j+1/2)
    h1 = new double[a1];
    h2 = new double[a1];
    h12= new double[a1];
    h22= new double[a1];
    h1h= new double[a1];
    h2h= new double[a1];
    hh = new double[a1];
    dh1 = new double[a1];

//sum of acceleration associated with gravity, earth's rotation
    gr  = new double[a1];
    
// ---- components of earth's rotation frequency
    Omega1 = new double[a2];
    Omega2 = new double[a2];
    Omega1h= new double[a2];
    Omega2h= new double[a2];

/* geographic colatitude and longitude */
    thetag = new double[a2];
    phig =   new double[a2];

// ---- solar zenith angles
    zenith = new double[a2];

// ---- neutral ionization rates
    qi = new double**[ym];
    qib= new double**[ym];
    qit= new double*[ym];

    //production and loss rates
    Ps = new double**[ym];
    Ls = new double**[ym];

//neutral EUV heating rates
    Qeuv = new double**[ym];

//exothermic heating rates
    Qch = new double**[ym];

//neutral cooling rate
    Cn = new double*[ym];

// ---- Electron photoelectron heating & electron cooling rates
    Qe = new double*[ym];
    Ce = new double*[ym];

    for (j=ys; j<ys+ym; j++) {
        yj=j-ys;

        qi[yj] = new double*[xm];
        qib[yj]= new double*[xm];
        qit[yj]= new double[xm];

        Qeuv[yj] = new double*[xm];
        Qch[yj] = new double*[xm];
        Cn[yj] = new double[xm];

        Ps[yj] = new double*[xm];
        Ls[yj] = new double*[xm];

        Qe[yj] = new double[xm];
        Ce[yj] = new double[xm];

        for (i=xs; i<xs+xm; i++) {
            xi=i-xs;

            qi[yj][xi] = new double[5];
            qib[yj][xi]= new double[16];

            //EUV heating of O, O2, N2
            Qeuv[yj][xi] = new double[3];

            //chemical heating of O+, O2+, N2+, H+, NO+, O, O2, N2, H
            Qch[yj][xi] = new double[slm];

            //sl ion species, sm major neutral & 2 minor neutral species NO & N
            Ps[yj][xi] = new double[slm+2];
            Ls[yj][xi] = new double[slm+2];
        }
    }

    return;
}
