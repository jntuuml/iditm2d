//Global scope function definition
//#ifndef INC_FUNCDEF_H
//#define INC_FUNCDEF_H

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include "param.h"

int input_param(char*,AppCtx*);
int input_prev(DM,Field**,Fieldu**,char*,char*);

void array_alloc(int,int,int,int);
void array_init(int,int,int,int);
void array_dealloc(DM);

int initialize(DM,Vec,AppCtx,char*);

void const_normalize();

void curvicoords();

int dipole_magnetic(DM,Vec,char*);

void geocoords();

int output_param(ofstream&,char*,int);
int output_solution(DM,Field**,Fieldu**,double,PetscInt,char*);

int euvflux_seg(DM,char*);

void update_timedate();

void solar_zenith();

void photoionization(Fieldu**,int,int,int,int);

void prod_loss_rates(Field**,Fieldu**,int,int,int,int);

void collision_freq(Field**,Fieldu**,int,int,int,int);

void ele_heating_rate(Field**,Fieldu**,int,int,int,int);
void ele_cooling_rate(Field**,Fieldu**,int,int,int,int);
void ele_velocity(Field**,Fieldu**,int,int,int,int);

void thermal_conductivity(Field**,Fieldu**,int,int,int,int);

void neu_cooling_rate(Fieldu**,int,int,int,int);

void jparallel(Field**,Fieldu**,int,int,int,int);

void top_bc(double);

void geomag(double&, double&, double&, double&, double&, double&, int);
void sphcar(double&, double&, double&, double&, double&, double&, int);
void vec_sphcar(double, double, double, double, double, double&, double&, double&);
int sun(int, int, double, double&, double&, double&, double&);
int dayno(int, int, int);
double mag_longitude();
double cerfc(double);
double lgrg(double x[], double y[], int, double);
double expon_integral(double, int);

void filter(Field**, int, int, int, int);
void shapiro(int, double*, int);

PetscErrorCode advance_neutrals(DM,Vec,Field**,double);

double funcvec_part(int, int, int, int, int, Field**, Field**, Fieldu**);

//#endif  /* INC FUNCDEF */
