/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

void eig_poisson1D(double* eigval, int *la)
{
  int i;
  double val;
  for (i = 0; i < *la; i++) {
    val = (1.0 * i + 1.0) * M_PI_2 * (1.0 / (*la + 1));
    eigval[i] = sin(val);
    eigval[i] = 4 * eigval[i] * eigval[i];
  }
}

double eigmax_poisson1D(int *la)
{
  double h = (*la) * (1.0 /((*la) + 1));
  return 4.0 * pow(sin(M_PI * h / 2.0), 2);
}

double eigmin_poisson1D(int *la)
{
  double h = (1.0 /((*la) + 1));
  return 4.0 * pow(sin(M_PI * h / 2.0), 2);
}

double richardson_alpha_opt(int *la)
{
  return 2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  cblas_dgbmv(CblasColMajor, CblasConjTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);
  cblas_dscal(*la, -1.0, resvec, 1);
  cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);


  double releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);

  while (releres > (*tol) && *nbite < *maxit)
  {
    cblas_daxpy(*la, *alpha_rich, resvec, 1, X, 1);
    cblas_dgbmv(CblasColMajor, CblasConjTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);
    cblas_dscal(*la, -1.0, resvec, 1);
    cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);

    releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);

    *nbite += 1;
  }
  printf("\n\n\n iter richardson = %d\n\n\n", *nbite); 
  }

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{
  for (int j=0;j<(*la);j++)
  {

    int k = j * (*lab);
    MB[k + *kv] = 1 / AB[k + *kv];    

  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite, int resolv)
{
  cblas_dgbmv(CblasColMajor, CblasConjTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);
  cblas_dscal(*la, -1.0, resvec, 1);
  cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);

  double releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);

  while (releres > (*tol) && *nbite < *maxit)
  {
    if (resolv == 1)
    {
      cblas_dgemv(CblasColMajor, CblasConjTrans, *la, *la, 1.0, MB, *la, resvec, 1, 1.0, X, 1);   // Xk = Xk + M(gauss) * r
    }
    else {cblas_dgbmv(CblasColMajor, CblasConjTrans, *la, *la, *kl, *ku, 1.0, MB, *lab, resvec, 1, 1.0, X, 1);}   // Xk = Xk + M(jacobi) * r

    cblas_dgbmv(CblasColMajor, CblasConjTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);
    cblas_dscal(*la, -1.0, resvec, 1);
    cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);

    releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);

    *nbite += 1;

  }
  printf("\n\n\n iter richardson = %d\n\n\n", *nbite);
}

