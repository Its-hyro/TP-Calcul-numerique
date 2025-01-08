/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/


#include "lib_poisson1D.h"
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2
#define DGBMV_TEST 3  
#define LU_TEST 4  

int main(int argc,char *argv[])

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
  double *X_TEST; 

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }
  
  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  X_TEST=(double *) malloc(sizeof(double)*la); 

  // Set the grid points  
  set_grid_points_1D(X, &la);
  // Set the dense RHS
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  // Set the analytical solution
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  // Allocation de la matrice AB
  AB = (double *) malloc(sizeof(double)*lab*la);

  // Set the GB operator
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  // Test de la fonction dgbmv
  if (IMPLEM == DGBMV_TEST) {
    printf("\nTest de la fonction dgbmv pour Poisson 1D\n");
    
    // Initialisation du vecteur RHS avec des 1 pour le test
    for(int i = 0; i < la; i++) {
      RHS[i] = 1.0;
      X_TEST[i] = 0.0;  // Initialisation de X_TEST à 0
    }
    
    // Utilisation du format spécifique pour DGBMV
    set_GB_operator_colMajor_poisson1D_DGBMV(AB, &lab, &la, &kv);
    
    // Affichage de la matrice AB avant dgbmv
    printf("\nMatrice AB avant dgbmv :\n");
    for(int i = 0; i < lab; i++) {
      for(int j = 0; j < la; j++) {
        printf("%f ", AB[j*lab + i]);
      }
      printf("\n");
    }
    
    // Appel de dgbmv
    dgbmv_poisson1D(AB, RHS, X_TEST, &la, &lab, &ku, &kl, &kv);
    
    // Affichage des résultats
    printf("\nRésultats :\n");
    for(int i = 0; i < la; i++) {
      double expected = 2.0 * RHS[i]; 
      if(i > 0) expected -= RHS[i-1];  
      if(i < la-1) expected -= RHS[i+1];  
      printf("X_TEST[%d] = %f (attendu : %f)\n", i, X_TEST[i], expected);
    }
    
    // Vérification du résultat
    int valid = 1;
    for(int i = 0; i < la; i++) {
      double expected = 2.0 * RHS[i];
      if(i > 0) expected -= RHS[i-1];
      if(i < la-1) expected -= RHS[i+1];
      if(fabs(X_TEST[i] - expected) > 1e-10) {
        valid = 0;
        printf("Erreur à l'indice %d : obtenu = %f, attendu = %f\n", i, X_TEST[i], expected);
      }
    }
    
    if(valid) {
      printf("\nTest dgbmv_poisson1D : SUCCÈS\n");
    } else {
      printf("\nTest dgbmv_poisson1D : ÉCHEC\n");
    }
    // Sauvegarde de la solution
    write_vec(X_TEST, &la, "DGBMV_TEST.dat");
  } else if (IMPLEM == LU_TEST) {
    printf("\nTest de la factorisation LU pour matrice tridiagonale\n");
    if(test_dgbtrftridiag()) {
      printf("Test factorisation LU : SUCCÈS\n");
    } else {
      printf("Test factorisation LU : ÉCHEC\n");
    }
  } else {
    clock_t start, end;
    double cpu_time_used;

    printf("Solution with LAPACK\n");
    ipiv = (int *) calloc(la, sizeof(int));

    /* LU Factorization */
    if (IMPLEM == TRF) {
      start = clock();
      dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
      if (info == 0) {
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      }
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("\nTemps d'exécution (DGBTRF + DGBTRS) : %f secondes\n", cpu_time_used);
    }

    /* LU for tridiagonal matrix */
    if (IMPLEM == TRI) {
      start = clock();
      dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
      if (info == 0) {
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      }
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("\nTemps d'exécution (DGBTRFTRIDIAG + DGBTRS) : %f secondes\n", cpu_time_used);
    }

    /* Solve with dgbsv */
    if (IMPLEM == SV) {
      start = clock();
      dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("\nTemps d'exécution (DGBSV) : %f secondes\n", cpu_time_used);
    }
    // Sauvegarde de la solution
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
    write_xy(RHS, X, &la, "SOL.dat");

    /* Relative forward error */
    relres = relative_forward_error(RHS, EX_SOL, &la);
    
    printf("\nThe relative forward error is relres = %e\n",relres);
  }
  // Libération des allocations
  free(RHS);
  free(EX_SOL);
  free(X);
  free(X_TEST);
  free(AB);
  if (IMPLEM == TRF || IMPLEM == TRI) {
    free(ipiv);
  }
  printf("\n\n--------- End -----------\n");
}
