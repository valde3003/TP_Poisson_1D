/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la)
{
  double h = (*la) * (1.0 /(*la) + 1);
  return 4.0 * pow(sin(M_PI * h / 2.0), 2);
}

double eigmin_poisson1D(int *la)
{
  double h = 1.0 / (*la + 1);
  return 4.0 * pow(sin(M_PI * h / 2.0), 2);
}

double richardson_alpha_opt(int *la)
{
  return 2 / eigmax_poisson1D(la) + eigmin_poisson1D(la);
}
//Fonction a corriger cela n'affiche pas correctement resvec
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  double *tmp = (double *)malloc(sizeof(double) * (*la));

  cblas_dcopy(*la, RHS, 1, tmp, 1);

  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, X, 1);
  cblas_dscal(*la, -1.0, resvec, 1);
  cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);

  double relres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);

  *nbite=0;

  while (relres > (*tol) && *nbite < *maxit)
  {
    cblas_daxpy(*la, *alpha_rich, tmp, 1, X, 1);

    cblas_dcopy(*la, RHS, 1, tmp, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, X, 1);
    cblas_dscal(*la, -1.0, resvec, 1);
    cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);
 
    relres = cblas_dnrm2(*la, tmp, 1) / cblas_dnrm2(*la, RHS, 1);

    (*nbite)++; 
  }
  free(tmp);
  printf("\n\n\n iter richardson = %d\n\n\n", *nbite); 

}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{
  int offset = *kv; 

  for (int i = 0; i < *la; ++i) 
  {
    MB[i] = AB[offset + i * (*lab)]; 
  } 
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

