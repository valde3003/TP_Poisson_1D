/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;
  double *LU;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=15;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS_direct.dat");
  write_vec(EX_SOL, &la, "EX_SOL_direct.dat");
  write_vec(X, &la, "X_grid_direct.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

 
  AB = (double *) malloc(sizeof(double)*lab*la);
  LU = (double *) malloc(sizeof(double)*lab*la);
  double *Y = (double *) malloc(sizeof(double) * la);
  //EXO 4
  for (jj = 0; jj < la; jj++) 
  {
    Y[jj] = 0.0;
  }

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_GB_operator_colMajor_poisson1D(LU, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_direct.dat");

  //Fonction CBLAS dgbmv avec cette matrice 
  cblas_dgbmv(CblasColMajor, CblasConjTrans, la, la, kl, ku, 1.0, AB+1, lab, EX_SOL, 1, 0.0, Y, 1);

  //Méthode de validation avec l'erreur relative 
  double norm = cblas_dnrm2(la, RHS, 1);
  cblas_daxpy(la, -1.0, Y, 1, RHS, 1);
  double err = cblas_dnrm2(la, RHS, 1);
  double relative_error = err / norm;

  printf("Erreur relative : %e\n", relative_error);


  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));
  //Exo 5
  //Avec dgbtrf/dgbtrs
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  clock_t dgbtrs_begin = clock();
  info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, la, la, kl, ku, AB, lab, ipiv);
  info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
  clock_t dgbtrs_end = clock();

  write_vec(RHS, &la, "dgbtrf.dat");

  //time
  double dgbtrs_delta = (double) (dgbtrs_end - dgbtrs_begin) / CLOCKS_PER_SEC;
  printf("Time DGBTRF + DGBTRS = %f sec\n", dgbtrs_delta);

  //Avec dgbsv
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  clock_t dgbsv_begin = clock();
  info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
  clock_t dgbsv_end = clock();
  write_vec(RHS, &la, "dgbsv.dat");

  // Time
  double dgbsv_delta = (double) (dgbsv_end - dgbsv_begin) / CLOCKS_PER_SEC;
  printf("Time DGBSV = %f sec\n", dgbsv_delta);
  
  
  /*write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL_direct.dat");*/
  
  //EXO 6
  /* LU Factorization */
  clock_t LU_begin = clock();
  factorisation_LU(AB, &lab, &la, &kv);
  clock_t LU_end = clock();
  write_GB_operator_colMajor_poisson1D(LU, &lab, &la, "LU_fac.dat");

  double LU_delta = (double) (LU_end - LU_begin) / CLOCKS_PER_SEC;
  printf("Time LU factorisation = %f sec\n", LU_delta);

  // LU avec LAPACK dgbtrf
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  clock_t dgbtrf_begin = clock();
  LAPACKE_dgbtrf(LAPACK_COL_MAJOR, la, la, kl, ku, AB, lab, ipiv);
  clock_t dgbtrf_end = clock();
  write_GB_operator_colMajor_poisson1D(LU, &lab, &la, "LAPACKE_dgbtrf.dat");

  double LU_dgbtrf = (double) (LU_end - LU_begin) / CLOCKS_PER_SEC;
  printf("Time LU factorisation LAPACK = %f sec\n", LU_dgbtrf);


  //Méthode de validation 
  norm = cblas_dnrm2(la, RHS, 1);
  cblas_daxpy(la*lab, -1, AB, 1, LU, 1);
  double err2 = cblas_dnrm2(la, LU, 1);
  double relative_error2 = err2 / norm;

  printf("Erreur relative : %e\n", relative_error2);

  /* Relative forward error */
  //relres = relative_forward_error(RHS, EX_SOL, &la);
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
