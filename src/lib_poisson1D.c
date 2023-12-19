/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv)
{
  int i,j,k;
  int diag = *kv + *kv + 1;  

    for (j = 0; j < *la; j++) 
    {
        k = j * (*lab); 
        for (i = 0; i < *lab; i++) 
        {
            AB[k + i] = 0.0;
        }
        if (j != 0) 
        {
            AB[k] = -1.0;
        }
        AB[k + *kv] = 2.0;
        if (j != *la - 1) 
        {
            AB[k + diag - 1] = -1.0;
        }
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
    if (RHS == NULL || la == NULL || BC0 == NULL || BC1 == NULL) 
    {
        printf("Pointeur null\n");
        return;
    }

    RHS[0] = *BC0;

    for (int i = 1; i < *la - 1; i++) 
    {
        RHS[i] = 0.0;
    }

    if (*la > 1) 
    {
        RHS[*la - 1] = *BC1;
    }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1)
{
  double h = 1.0 / (*la + 1);
  double delta_T = *BC1 - *BC0;

  for (int i = 0; i < *la; i++) 
  {
      X[i] = (i + 1) * h;
      EX_SOL[i] = *BC0 + X[i] * delta_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab){
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}

void factorisation_LU(double* AB, int *lab, int *la, int *kv) 
{
    int diag = 2 * (*kv) + 1;

    int i; 

    for (i = 0; i < *la - 1; ++i) 
    {
        int k = i * (*lab);

        if (AB[k + *kv] == 0) 
        {
            printf("Erreur: Élément diagonal nul.\n");
            return;
        }
        
        int j;

        for ( j = i + 1; j < *la && j <= i + *kv; ++j) 
        {
            int kj = j * (*lab);
            double m = AB[kj + *kv - (j - i)] /= AB[k + *kv];

            for (int d = 1; d < diag - (j - i); ++d) 
            {
                AB[kj + *kv - (j - i) + d] -= m * AB[k + *kv + d];
            }
        }
    }
}


