/************************************************************************
 header file: param.h
 Purpose: declare variables and arrays for global use (as externals)

 By Jiannan Tu 9/12/2013
************************************************************************/
#ifndef INC_PARAM_H
#define INC_PARAM_H

#include "petscsnes.h"

#define kbmp 8.25481286619634103e3

/* variables of interest to be solved at each grid point */
typedef struct {
   PetscScalar fx[29];
} Field;

/* parameters to be used in solving variables of interest */
typedef struct {
   double  B0[2];
   double  nu[90];
   double  betae;
   double  Jpar[2];
   double  lamdae;
   double  lamdas[9];
   double  fu[22];
} Fieldu;

/* User defined variables. Used in solving for variables of interest */
typedef struct {
  DM          dau;
  Vec         U, Xn, Xn1;
  PetscInt    nt;
  PetscReal   ti;
  PetscInt    nouta,nosub;
  char        workdir[150];
  PetscInt    ini;
  char        prefln[150];
  char        preflm[150];
  char        preflo[150];
  time_t      start_t;
} AppCtx;

//physical constants
const double q=1.6022e-19, ge=9.80665, Re=6371.2e3;
const double kb=1.3807e-23, eps0=8.85e-12, mu0=1.25663706e-6;
const double gm=3.984997527e14;
const double w0=7.27220517e-5;

//coefficients for transformation between the geographic and geomagnetic coordinates
const double a11=0.306892902, a12=-0.935005784, a13=-0.177709982;
const double a21=0.950129092, a22=0.311856776;
const double a31=5.54200634e-2, a32=-0.168847427, a33=0.984082878;

/* proton mass in kg */
const double mp=1.6726e-27;

//electron mass normalized by H mass, i.e., mass in AMU
const double me=5.455077356e-4;

/* mass of O+, O2+, N2+, H+, NO+, O, O2, N2, H, NO, N normalized 
 * by proton mass, i.e. mass in atomic mass unit (AMU) */
const double ms[11]={16.0, 1.0, 32.0, 28.0, 30.0, 16.0, 32.0, 28.0, 1.0, 30.0, 14.0};

//average photoelectron energies (eps, in eV) in the parameterization of 
//photoelectron heating
const double epi[16]={169.4, 99.07, 69.43, 54.50, 48.37, 43.63, 45.93, 40.88,
                      40.81, 37.38, 33.68, 34.22, 29.62, 29.65, 25.38, 23.69};

//coefficients used in the parameterization of photoelectron heating
const double cc1[7]={1.468, 9.229e-1, 4.956e-2, -1.897e-2, -3.934e-3, -2.634e-4, 
                     -5.980e-6};
const double cc2[7]={1.020, 1.540e-2, -6.858e-3, -8.528e-3, -2.052e-3, -1.634e-4,
                     -4.314e-6};

//year, month, date, and UT in second
extern int    iyr, mon, idate;
    
extern double sec;

extern double f107, f107a, Ap;

extern int    ntot, nout, npre;

//lower and upper radial distance
extern double rb, ru;

extern double LT;

//number of ion and major neutral species
extern int sl, sm, slm;

//# of mesh cells along x1 & x2 curvilinear coordinates
extern int  N1, N2;

//normalization parameters
extern double r0, n0, B0, v0, t0, E0, g0, p0, j00, e0, T0, q0, lamda0, beta0;

//start up mode (1/0): from IRI/MSIS or from results of a previous run
extern int  smod;

//altitude location where to stop calculating prod/loss/photoionization
extern int ipc;

//two stage imposed perturbation magnetic field
extern double bimp, bimpf;

//time to increase v_theta from 0 to Vimp or from Vimp to Bimpf
extern double tsl, TL;

extern double lat[4];

extern double pi, pi2, deg, rad;

//change-mass ratios, normalized speed of light
extern double qms[7], Omegae;

//a1=N1+1, a2=N2+1, a3=N3+1
extern int  a1, a2, a3, aa;

//normalized elementary charge, square of speed of light, & product of two &
//normalized rotation frequency of Earth, normalized acceleration of gravity
extern double  e, w0n, gen;

//r, theta, phi in magnetic polar coordinates
extern double *rr, *theta, *phi;

//geographic colatitude, geographic longitude
extern double *thetag, *phig;

//scale factors of curvilinear coordinates
extern double *h1, *h2, *h1h, *h2h, *h12, *h22, *hh, *dh1, *rh, *zh;

//components of earth's rotation frequency in curvilinear coordinates
extern double *Omega1, *Omega2, *Omega1h, *Omega2h;

extern double *gr;

//EUV emission flux, photo-absorption and photo-ionization cross sections & photon energies
//in given wavelength bins
extern double euvflux[37], segabs[37][5], segion[37][7], pene[37];

//for nighttime photo-absorption & ionization
extern double euvfluxn[4], segabsn[4][4], segionn[4][4];

//solar zenith
extern double *zenith;

/* photoionization rates for individual species, photoionization rate in 
   individual wavelength bins summed over all species, & total rate */
extern double ***qi, ***qib, **qit;

extern double ***Ps, ***Ls;

//electron photoelectron heating of electrons, electron cooling rates
extern double **Qe, **Ce;

//neutral EUV heating rate, exothermic heating rates for ions and neutrals
extern double ***Qeuv, ***Qch;

//neutral cooling rates
extern double **Cn;

//time step & its half and quart, & time step for E&B fields in subionosphere
extern double dt, dth, dt2;

//coefficients for collision frequencies, multiplied by n0*t0
extern double coe[5], coiO[8], coiO2[8], coiN2[8], coiH[8], coiNO[8], con[6];

/* nvt = 0/1: velocity or perturbation magnetic field as top BC; */
extern int nvt;

//velocity or perturbation magnetic field at the top boundary
extern double *btop;

extern double d;

extern double cshap;

extern double  terms[360][80][5];

#endif  /* INC PARAM */
