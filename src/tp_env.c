/******************************************/
/* tp_env.c                               */
/* This file contains a main function to  */
/* test the environment of compilation    */
/******************************************/
#include "tp_env.h"

int main(int argc,char *argv[])
{
  printf("--------- Test environment of execution for Practical exercises of Numerical Algorithmics ---------\n\n");
  
  // Test des valeurs de base
  printf("=== Test des valeurs fondamentales ===\n");
  printf("The exponantial value is e = %f \n",M_E);
  printf("The maximum single precision value from values.h is maxfloat = %e \n",MAXFLOAT);
  printf("The maximum single precision value from float.h is flt_max = %e \n",FLT_MAX);
  printf("The maximum double precision value from float.h is dbl_max = %e \n",DBL_MAX);
  printf("The epsilon in single precision value from float.h is flt_epsilon = %e \n",FLT_EPSILON);
  printf("The epsilon in double precision value from float.h is dbl_epsilon = %e \n",DBL_EPSILON);

  printf("\n=== Test of BLAS/LAPACK environment ===\n");

  // Test BLAS niveau 1 (opérations vecteur-vecteur)
  double x[5], y[5];
  int ii;
  for (ii=0; ii<5; ii++){
    x[ii]=ii+1; 
    y[ii]=ii+6;
    printf("x[%d] = %lf, y[%d] = %lf\n",ii,x[ii],ii,y[ii]);
  }

  // Test DCOPY
  printf("\nTest DCOPY y <- x \n");
  cblas_dcopy(5,x,1,y,1);
  for (ii=0; ii<5; ii++){
    printf("y[%d] = %lf\n",ii,y[ii]);
  }

  // Test DDOT (produit scalaire)
  printf("\nTest DDOT (produit scalaire)\n");
  double dot = cblas_ddot(5,x,1,y,1);
  printf("x.y = %lf\n", dot);

  // Test DAXPY (y = ax + y)
  printf("\nTest DAXPY (y = 2x + y)\n");
  double alpha = 2.0;
  cblas_daxpy(5,alpha,x,1,y,1);
  for (ii=0; ii<5; ii++){
    printf("y[%d] = %lf\n",ii,y[ii]);
  }

  // Test DNRM2 (norme euclidienne)
  printf("\nTest DNRM2 (norme euclidienne)\n");
  double nrm = cblas_dnrm2(5,x,1);
  printf("||x|| = %lf\n", nrm);

  printf("\n--------- Tests completed successfully -----------\n");
  return 0;
}
