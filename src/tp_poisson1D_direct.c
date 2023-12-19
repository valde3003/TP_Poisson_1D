/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <time.h>

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=10;
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
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  //EXO 4
  AB = (double *) malloc(sizeof(double)*lab*la);
  double *Y = (double *) malloc(sizeof(double) * la);

  for (jj = 0; jj < la; jj++) 
  {
    Y[jj] = 0.0;
  }

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  //Fonction CBLAS dgbmv avec cette matrice 
  cblas_dgbmv(CblasColMajor, CblasConjNoTrans, la, la, kl, ku, 1.0, AB+1, lab, EX_SOL, 1, 0.0, Y, 1);

  //Méthode de validation avec l'erreur relative 
  double norm = cblas_dnrm2(la, RHS, 1);
  cblas_daxpy(la, -1.0, Y, 1, RHS, 1);
  double err = cblas_dnrm2(la, RHS, 1);
  double relative_error = err / norm;

  printf("Erreur relative : %e\n", relative_error);

  printf("Solution with LAPACK\n");
  //EXO 5
  //Fonction dgbtrf_
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);

  top1 = clock();
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  fprintf(dgbtrf, "%f\n",((double)(clock() - top1)/CLOCKS_PER_SEC));
  
  //Fonction dgbsv
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(EX_RHS,&la,&T0,&T1);

  top2 = clock();
  dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, EX_RHS, &la, &info);
  fprintf(direct, "%f\n",((double)(clock() - top2)/CLOCKS_PER_SEC));

  //EXO 6
  /* LU Factorization */

  LU = (double *) malloc(sizeof(double)*lab*la);

  factorisation_LU(LU, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");

  //Avec dgbtrf
  info=0;
  ipiv = (int *) calloc(la, sizeof(int));

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  //Méthode de validation 
  double norm = cblas_dnrm2(la, RHS, 1);
  cblas_daxpy(la*lab, -1, AB, 1, LU, 1);
  double err2 = cblas_dnrm2(la, LU, 1);
  double relative_error2 = err2 / norm;

  printf("Erreur relative : %e\n", relative_error2);
  

  write_xy(RHS, X, &la, "SOL.dat");
  

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  // ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  
  /* Solution (Triangular) */
  /*if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }*/

  /* It can also be solved with dgbsv */
  // TODO : use dgbsv

  //write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  // TODO : Compute relative norm of the residual
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(Y);
  free(LU);
  printf("\n\n--------- End -----------\n");
}
