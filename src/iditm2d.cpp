
static char help[] = "Time-dependent PDE in 2d. \n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#include "param.h"
#include "param_def.h"
#include "funcdef.h"

/* Function and Jacobian calculation routine */
extern PetscErrorCode FormIFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormIJacobian(SNES,Vec,Mat,Mat,void*);
/* calls functions that calculate parameters and outputs results */
extern PetscErrorCode MyMonitor(SNES,PetscInt,PetscReal,Vec,void*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Vec            X,rhs;              /* solution, parameter,residual vectors */
  Mat            J;                    /* Jacobian matrices */
  PetscErrorCode ierr;
  DM             da;
  AppCtx         user;              /* user-defined work context */
  SNES           snes;

  PetscMPIInt    rank;
  PetscInt       nst,nt,dof=127,np,i;
  char           workdir[150], buff[250];
  char           fname[150];
  time_t         start_t, end_t;
  struct         tm *now;
  ofstream       logfstr;

  //get current time
  start_t=time(NULL);
  now = localtime( &start_t );
  user.start_t=start_t;

  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  if(getcwd(workdir, sizeof(buff)) == NULL) {
    perror("getcwd() error");
    return -1;
  };

  /* Inputs */
  if (input_param(workdir,&user) < 0) return -1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_PERIODIC,DMDA_STENCIL_STAR,
           a1,a2,1,PETSC_DECIDE,a3,1,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,NULL,NULL,NULL,NULL,NULL,&np,NULL,NULL,NULL,NULL,NULL,
          NULL,NULL);CHKERRQ(ierr);
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_PERIODIC,DMDA_STENCIL_STAR,
            a1,a2,1,np,dof,1,NULL,NULL,&user.dau);CHKERRQ(ierr);

  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);
  ierr = DMSetUp(user.dau);CHKERRQ(ierr);

  if (!rank) {
      if (np > N2) {
        cout<<"Number of processor n = "<<np<<" must be < number of grids N2 ="<<N2
            <<" on theta axis"<<endl;
        cout<<"The current version of the code limits 1 processor on radial axis "
            <<"to make it easy to process files I/O"<<endl;
        MPI_Abort(PETSC_COMM_WORLD,-1);
      }
  }

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(da,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&rhs);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&user.Xn);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&user.Xn1);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user.dau,&user.U);CHKERRQ(ierr);

  if (initialize(da,X,user,workdir)<0)
      MPI_Abort(PETSC_COMM_WORLD,-1);

  //strncpy(workdir, "/project/uml_jiannan_tu", 150);

  /* calculate background magnetic field using either dipole or eccentric 
   * dipole magnetic field (first three harmonic terms of IGRF-11). */
  dipole_magnetic(user.dau,user.U,workdir);

  strncpy(user.workdir, workdir, 150);

/* 
 * write run setup to a file
*/
  if (!rank) {
    //open the log file to write run records
    strncpy(fname, workdir, 150);
    strcat(fname, "/output/iditm2d.log");
    logfstr.open(fname,fstream::out);
    if(!logfstr) {
        cout << "Can't open file " <<fname << endl;
        MPI_Abort(PETSC_COMM_WORLD,-1);
    }
    logfstr <<"Date and Time when program run started: "
            <<now->tm_year+1900<<"-"<<now->tm_mon+1<<"-"<<now->tm_mday<< " "
            <<now->tm_hour<< ":"<<now->tm_min<<":"<<now->tm_sec<< endl<<endl;
    logfstr << "Run command and options:" << endl;
    logfstr<<argv[0]<<endl;
    for (i = 1; i < argc-1; i=i+2) {
        logfstr<<argv[i]<<" "<<argv[i+1]<<endl;
    }
    logfstr<<endl;

    output_param(logfstr,workdir,user.nt);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,rhs,FormIFunction,&user);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSetSolution(snes,X);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set Jacobian evaluation routine
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMSetMatrixPreallocateOnly(da,PETSC_TRUE);
  ierr = DMSetMatType(da,MATAIJ);CHKERRQ(ierr);
  ierr = DMCreateMatrix(da,&J);CHKERRQ(ierr);
  ierr = MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,FormIJacobian,&user);CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Sets various TS parameters from user options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  if (smod == 1) nst=npre;
  else nst=0;
  ntot=ntot+nst;

  user.ini=0;
  if (MyMonitor(snes,nst,double(nst)*dt,user.Xn,&user)<0) return -1;
  user.ini=1;

  for (nt=nst+1; nt<=ntot; nt++) {
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = SNESSolve(snes,PETSC_NULL,X);CHKERRQ(ierr);

    if (MyMonitor(snes,nt,double(nt)*dt,X,&user)<0) MPI_Abort(PETSC_COMM_WORLD,-1);
  }

  if (!rank) {
      end_t=time(NULL);
      now=localtime(&end_t);

      logfstr <<endl<<"Date and Time when program execution ended: "
            <<now->tm_year+1900<<"-"<<now->tm_mon+1<<"-"<<now->tm_mday
            <<" " <<now->tm_hour<<":"<<now->tm_min<< ":"<<now->tm_sec<<endl;
      logfstr.close();
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Free work space.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  array_dealloc(da);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  //ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&rhs);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  //ierr = DMDestroy(&dau);CHKERRQ(ierr);

  ierr = PetscFinalize();
  PetscFunctionReturn(0);
}
