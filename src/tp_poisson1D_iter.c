/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <time.h>

#define ALPHA 0
#define JAC 1
#define GS 2

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;
  double *MB;
  
  double temp, relres;

  double opt_alpha;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  /* Size of the problem */
  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double)); 
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS_iter.dat");
  write_vec(EX_SOL, &la, "EX_SOL_iter.dat");
  write_vec(X, &la, "X_grid_iter.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  /* uncomment the following to check matrix A */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_iter.dat");
  
  /********************************************/
  /* Solution (Richardson with optimal alpha) */
  printf("Richardson with optimal alpha\n\n");

  /* Computation of optimum alpha */
  opt_alpha = richardson_alpha_opt(&la);
  printf("Optimal alpha for simple Richardson iteration is : %lf\n",opt_alpha);

  /* Solve */
  double tol=1e-3;
  int maxit=1000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  //if (IMPLEM == ALPHA) {
  clock_t richardson_begin = clock();
  richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  clock_t richardson_end = clock();

  double richardson_delta = (double) (richardson_end - richardson_begin) / CLOCKS_PER_SEC;
  printf("Time RICHARDSON = %f sec\n", richardson_delta);

  write_vec(SOL, &la, "SOL_RICHARDSON.dat");
  write_vec(resvec, &nbite, "RESVEC_RICHARDSON.dat");

  //Calcul erreur relative
  temp = cblas_dnrm2(la, EX_SOL, 1);
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_dnrm2(la, EX_SOL, 1);
  relres = relres / temp;

  printf("\nThe relative forward error for RICHARDSON is relres = %e\n\n",relres);

  /* Richardson General Tridiag */

  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  kv = 1;
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(lab)*la);
  SOL=(double *) calloc(la, sizeof(double));
  resvec=(double *) calloc(maxit, sizeof(double));
  
  clock_t jacobi_begin = clock();
  extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite, 0);
  clock_t jacobi_end = clock();

  double jacobi_delta = (double) (jacobi_end - jacobi_begin) / CLOCKS_PER_SEC;
  printf("Time JAKOBI = %f sec\n", jacobi_delta);

  write_vec(SOL, &la, "SOL_JAKOBI.dat");
  write_vec(resvec, &nbite, "RESVEC_JAKOBI.dat");
  
  //Calcul erreur relative
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  temp = cblas_dnrm2(la, EX_SOL, 1);
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_dnrm2(la, EX_SOL, 1);
  relres = relres / temp;

  printf("\nThe relative forward error for JACOBI is relres = %e\n\n",relres);

  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End -----------\n");
}
