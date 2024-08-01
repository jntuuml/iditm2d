/************************************************************************
 header file: param_def.h
 Purpose: declair variables and arrays

 By Jiannan Tu 9/12/2013
************************************************************************/
//year, month, date, and UT in second
int    iyr, mon, idate;
    
double sec;

double f107, f107a, Ap;

int    ntot, nout, npre;

//lower and upper radial distance
double rb, ru;

//Initial magnetic local time
double LT;

//number of ion and major neutral species
int sl, sm, slm;

//# of mesh cells along x1 & x2 curvilinear coordinates
int  N1, N2;

//normalization parameters
double r0, n0, B0, v0, t0, E0, g0, p0, j00, e0, T0, q0, lamda0, beta0;

//time step in ionosphere and Earth-ionosphere waveguide
double dt, dth, dt2;

//start up mode (1/0): from IRI/MSIS or from results of a previous run
int  smod;

//altitude location where to stop calculating prod/loss/photoionization
int ipc;

//two stage imposed perturbation magnetic field
double bimp, bimpf;

//time to increase v_theta from 0 to Vimp or from Vimp to Bimpf
double tsl, TL;

//dayside & nightside latitudes (deg) of the polar cap boundary
double lat[4];

double pi, pi2, deg, rad;

double qms[7], Omegae;

int  a1, a2, a3, aa;

double  e, w0n, gen;

double *rr, *theta, *phi;

double *thetag, *phig;

double *h1, *h2, *h1h, *h2h, *h12, *h22, *hh, *dh1, *rh, *zh;

//components of earth's rotation frequency in curvilinear coordinates
double *Omega1, *Omega2, *Omega1h, *Omega2h;

double *gr;

double euvflux[37], segabs[37][5], segion[37][7], pene[37];
double euvfluxn[4], segabsn[4][4], segionn[4][4];

double *zenith;

double ***qi, ***qib, **qit;

double ***Ps, ***Ls;

double ***nue, ****nus;

double **Qe, **Ce;

double ***Qeuv, ***Qch;

double **Cn;

double coe[5], coiO[8], coiO2[8], coiN2[8], coiH[8], coiNO[8], con[6];

int nvt;

double *btop;

double d;

double cshap;

double terms[360][80][5];

