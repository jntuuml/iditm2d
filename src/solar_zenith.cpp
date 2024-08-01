/*************************************************************************
  Function solar_zenith
      calculate solar zenith angle in rad

  Input: parms - structure for system parameters
 *       t     - time in sec relative to the start of the simulation
 * Output: none. Result is solar zenith angles on all the grids in the global
 *               pointer *zenith to an one dimensional array

  By Jiannan Tu
  March 1, 2011, 1/13/2014
*************************************************************************/
#include <cmath>

#include "funcdef.h"
#include "param.h"

using namespace std;

void solar_zenith()
{
    int    j;
    double eot, bt, phin, lt, lst, sha;
    double gst, slong, srasn, sdec, hr;
    int    iday;

/* day number of the year*/
    iday=dayno(iyr, mon, idate);

/* solar declination angle */
    sun(iyr, iday, sec, gst, slong, srasn, sdec);

    hr=sec/3600.0;

    //equation of time in minutes (bt in radians & constant rai = 360/365/deg)
    bt=rad*(iday-81.0);
    eot=9.87*sin(bt)-7.53*cos(bt)-1.5*sin(bt);

//#pragma omp parallel for private (i,phin,lt,lst,sha)
    for (j = 0; j <= N2; j++) {
        //geographic longitude in degrees
        phin=phig[j]*deg;

        //local time in hours = UT+longitude/15
        lt=hr+phin/15.0;
        if (lt >= 24.0) lt=lt-24.0;

        /* local solar time in hours by adjusting local time with time correction
         * factor */
        lst=lt+eot/60.0;

        //solar hour angle in radians
        sha=0.2617993878*(lst-12.0);

        //finally calculate solar zenith angle in radians
        zenith[j]=acos(sin(sdec)*cos(thetag[j])+cos(sdec)*sin(thetag[j])*cos(sha));
    }
    
    return;
}

