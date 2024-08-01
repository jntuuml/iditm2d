/***********************************************************************
 function mag_longitude
   calculate magnetic longitude of the dayside magnetic meridian

   Author Jiannan Tu
   Date 9/24/2013, rev. 12/4/2013
***********************************************************************/
#include "param.h"
#include "funcdef.h"

double mag_longitude()
{
    double rg, thet, phin;
    double xgeo, ygeo, zgeo, xmag, ymag, zmag;

/* find the magnetic meridian passing the point with given local time, on the 
 * earth's surface, and at the equator. */

    //calculate geographic longitude at the point [Re, theta=pi/2, long=15*(LT-UT)]
    phin=15.0*(LT-sec)/3600.0;
    if (phin < 0.0) phin=phin+360.0;
    phin=phin*rad;

    //convert that location to Cartesian coordinates
    rg=Re;
    thet=pi/2.0;
    sphcar(rg, thet, phin, xgeo, ygeo, zgeo, 1);

    //GEO Cartesian to GEOMAG Cartesian
    geomag(xgeo, ygeo, zgeo, xmag, ymag, zmag, 1);

    //GEOMAG Cartesian to GEOMAG spherical
    sphcar(rg, thet, phin, xmag, ymag, zmag, -1);

    return phin;
}

