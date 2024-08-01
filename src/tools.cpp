#include <cmath>
#include "param.h"

using namespace std;

/* conversion between geographic and geomagnetic coordinates 
 if dir > 0 from (xgeo, ygeo, zgeo) to (xmag, ymag, zmag)
 else vice verse
*/

void geomag(double &xgeo, double &ygeo, double &zgeo, double &xmag, double &ymag,
            double &zmag, int dir)
{
    if(dir > 0) {
        xmag=a11*xgeo+a12*ygeo+a13*zgeo;
        ymag=a21*xgeo+a22*ygeo;
        zmag=a31*xgeo+a32*ygeo+a33*zgeo;
    }
    else {
       xgeo=a11*xmag+a21*ymag+a31*zmag;
       ygeo=a12*xmag+a22*ymag+a32*zmag;
       zgeo=a13*xmag+a33*zmag;
    }

    return;
}

/*   CALCULATES CARTESIAN FIELD COMPONENTS FROM LOCAL SPHERICAL ONES
-----INPUT:   thet, phi - SPHERICAL ANGLES OF THE POINT IN RADIANS
              Br, Bth, Bph -  LOCAL SPHERICAL COMPONENTS OF THE FIELD
-----OUTPUT:  Bx, By, Bz - CARTESIAN COMPONENTS OF THE FIELD
*/

void vec_sphcar(double thet, double phis, double Br, double Bth, double Bph, 
              double &Bx, double &By, double &Bz)
{
    double S, C, SF, CF, BE;
      
    S=sin(thet);
    C=cos(thet);
    SF=sin(phis);
    CF=cos(phis);
    BE=Br*S+Bth*C;
    Bx=BE*CF-Bph*SF;
    By=BE*SF+Bph*CF;
    Bz=Br*C-Bth*S;

    return;
}

void sphcar(double &ra, double &thet, double &phis, double &x, double &y,
            double &z, int dir)
{
/*  CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICE VERSA
    (THETA AND PHI IN RADIANS).

                 dir > 0           dir < 0
-----INPUT:   dir,R,THETA,PHI     dir,X,Y,Z
-----OUTPUT:      X,Y,Z           R,THETA,PHI

  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 WHEN CONVERTING
        FROM CARTESIAN TO SPHERICAL COORDS (I.E., FOR J<0)
*/

    double sq;

    if (dir < 0) {
        sq=x*x+y*y;
        ra=sqrt(sq+z*z);

        if (sq == 0.0) {
            phis=0.0;
            if (z >= 0.0) thet=0.0;
            else thet=pi;
        }
        else {
            sq=sqrt(sq);
            phis=atan2(y, x);
            thet=atan2(sq, z);
            if (phis < 0.0) phis=phis+2.0*pi;
        }
    }
    else {
        sq=ra*sin(thet);
        x=sq*cos(phis);
        y=sq*sin(phis);
        z=ra*cos(thet);
    }

    return;
}
/***************************************************************************
 * Function sun
 *    Calculate sidereal, time and position of the sun
 *    Good for years 1901 through 2099. Accuracy 0.006 degree

 * Input: iyr (yyyy), IDAY (integers), and sec universal time in sec
 * output: gst - Greenwich mean sidereal time (radian). It's different from the 
 *               Greenwich mean time.
 *         slong - longitude along ecliptic (radian)
 *         srasn - right ascension (radian)
 *         sdec - declination of the sun (radian)

 * ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
 * RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
 * From book: In Introduction to Space Physics edited by M. G. Kivelson
             and C. T. Russell, pp.568, Cambridge University Press, 1995
***************************************************************************/
#include <iostream>
#include <cmath>

using namespace std;

int sun(int yr, int day, double secs, double &gst, double &slong, double &srasn,
        double &sdec)
{
    double dj, fday, t, vl, g, obliq, slp, sind, cosd, sc, sob;

    if (yr < 1901 || yr > 2099) {
        cout << "Year must be within the range of 1901 - 2099";
        return -1;
    }

    fday=secs/86400.0;
    dj=365.0*(yr-1900.0)+(yr-1901.0)/4.0+day+fday-0.5;
    t=dj/36525.0;
    vl=fmod(279.696678+0.9856473354*dj, 360.0);
    gst=fmod(279.690983+0.9856473354*dj+360.0*fday+180.0, 360.0)/deg;
    g=fmod(358.475845+0.985600267*dj, 360.0)/deg;
    slong=(vl+(1.91946-0.004789*t)*sin(g)+0.020094*sin(2.0*g))/deg;
    if (slong > pi2) slong=slong-pi2;
    if (slong < 0.0) slong=slong+pi2;
    obliq=(23.45229-0.0130125*t)/deg;
    sob=sin(obliq);
    slp=slong-9.923942126839758e-5;
    sind=sob*sin(slp);
    cosd=sqrt(1.0-sind*sind);
    sc=sind/cosd;
    sdec=atan(sc);
    srasn=pi-atan2(cos(obliq)/sob*sc, -cos(slp)/cosd);

    return 0;
}
/***********************************************************************
 function iditm_dayno
   get number of days in a year from Jan 1

   Author Jiannan Tu
   Date March 15, 2007, coverted to C, October 22, 2012
***********************************************************************/
int dayno(int yr, int mn, int dd)
{
    int iday=1, i;

    const int id[11]={31,59,90,120,151,181,212,243,273,304,334};

    if(mn ==1) {
        iday=dd;
        return iday;
    }
    else {
        for (i = 1; i < 11; i++) if (mn == i+2) iday=id[i]+dd;
    }

    if (mn > 2 && (yr % 4) == 0) iday=iday+1;

    return iday;
}
/***********************************************************************
 function lgrg(X,Y,N,T,Z)
    Eight-point Lagrange interpolation and extrapolation
    Given arrays X and Y, each of length N, and given a value of T,
    this routine returns a value Z. The interpolation or extrapolation
    is performed using lagrange's formula.
***********************************************************************/
double lgrg(double x[], double y[], int n, double t)
{
    int i, j, k, m;
    double s, z;

    for (i = 0; i < n; i++) {
        if (t == x[i]) {
            z=y[i];
            return z;
        }
    }

    z=0.0;

    i=0;
    while (x[i] < t && i < n) {
      i=i+1;
    }

    k=i-3;
    if (k < 0) k=0;
    m=i+3;
    if (m > n) m=n;

    for (i = k; i < m; i++) {
        s=1.0;
        for (j = k; j < m; j++) {
            if (j != i) s=s*(t-x[j])/(x[i]-x[j]);
        }
        z=z+s*y[i];
    }

    return z;
}
/*************************************************************************
  Function cerfc: calculate value of complementary error function

   Input: y

  Numerical recipes, Press (1992)
*************************************************************************/
#include <cmath>

double cerfc(double y)
{
    const double a=1.265512230, b=1.000023680,  c=0.374091960, dd=0.096784180;
    const double ee=-0.186288060, f=0.278868070, g=-1.135203980, hc=1.488515870;
    const double p=-0.822152230, qc=0.170872770;

    double t, x;

    t=1.0/(1.0+0.5*fabs(y));

    x=t*exp(-y*y-a+t*(b+t*(c+t*(dd+t*(ee+t*(f+t*(g+t*(hc+t*(p+t*qc)))))))));

    if (y < 0.0) x=2.0-x;

    return x;
}

/********************************************************************************
  Function expon_integral
    compute value of the second exponential integral E2(z) with five-point Bode's
    rule.

    see, e.g., Abramowitz and Stegun, Handbook of Mathematical Functions with
    Formulas, Graphs, and Mathematical Tables, 9th printing, National Bureau
    of Standards, Applied Mathematics Series, 55, Washington, D.C., 1970

   Input: z -  argument of function E2(z)
          n  - number of grid cells, must be multiple of 5

   return value: value of E2(z)

   Jiannan Tu
   1/25/2014
*************************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

double expon_integral(double z, int n)
{
    int    i;
    double h, x, E2;

    double *y = new double[n+1];
 
    if ((n != 4) && n % 20 != 0) {
        cout << "Number of data points must be 4 or multiple of 20!" << endl;
        delete[] y; return -1.0;
    }

    //E2(0)=1.0
    if(z == 0.0) {
        delete[] y; return 1.0;
    }

    h=1.0/double(n);

    /* values of integrand for the second exponential integral E2(z)=integrating exp(-z/x)
       from x=0 to x=1 */
    y[0]=0.0;
    for (i = 1; i <= n; i++) {
        x=h*double(i);
        y[i]=exp(-z/x);
    }

    //five-point Bode's rule for integration
    E2=0.0;
    for (i = 0; i < n-3; i=i+4) {
        E2=E2+7.0*(y[i]+y[i+4])+32.0*(y[i+1]+y[i+3])+12.0*y[i+2];
    }
    E2=E2*2.0*h/45.0;

    delete[] y;

    return E2;
}