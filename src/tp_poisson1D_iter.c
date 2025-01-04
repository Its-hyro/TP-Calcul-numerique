/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

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
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  /* uncomment the following to check matrix A */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  /********************************************/
  /* Solution (Richardson with optimal alpha) */

  /* Computation of optimum alpha */
  opt_alpha = richardson_alpha_opt(&la);
  printf("Optimal alpha for simple Richardson iteration is : %lf",opt_alpha); 

  /* Solve */
  double tol=1e-3;
  int maxit=1000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  if (IMPLEM == ALPHA) {
    richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    printf("\nRichardson :\n");
    printf("Nombre d'itérations : %d\n", nbite);
    printf("Résidu final : %e\n", resvec[nbite-1]);
    
    // Analyse de la convergence
    printf("\nAnalyse de la convergence :\n");
    printf("Iteration |    Résidu    | Ratio de convergence\n");
    printf("-----------------------------------------\n");
    for(int i = 0; i < nbite; i += nbite/5) {  // Affiche ~5 points
        double ratio = (i > 0) ? resvec[i]/resvec[i-1] : 0.0;
        printf("%9d | %11.1e | %19.2f\n", i, resvec[i], ratio);
    }
    printf("%9d | %11.1e | %19.2f\n", nbite-1, resvec[nbite-1], 
           resvec[nbite-1]/resvec[nbite-2]);
    
    // Affichage des solutions pour comparaison
    printf("\nComparaison des solutions :\n");
    printf("     X     |  Analytique  |  Richardson  |   Différence\n");
    printf("------------------------------------------------\n");
    for(int i = 0; i < la; i++) {
        printf("%9.6f | %11.6f | %11.6f | %11.6f\n", 
               X[i], EX_SOL[i], SOL[i], fabs(EX_SOL[i] - SOL[i]));
    }
    
    // Calcul de l'erreur par rapport à la solution analytique
    relres = relative_forward_error(SOL, EX_SOL, &la);
    printf("\nErreur relative par rapport à la solution analytique : %e\n", relres);
    
    // Sauvegarde des solutions pour comparaison
    write_vec(SOL, &la, "SOL_richardson.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
  }

  /* Richardson General Tridiag */

  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  kv = 1;
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(lab)*la);
  if (IMPLEM == JAC) {
    extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  } else if (IMPLEM == GS) {
    extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  }

  /* Solve with General Richardson */
  if (IMPLEM == JAC) {
    write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
    
    // Résolution avec Jacobi
    jacobi_tridiag(AB, RHS, SOL, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    
    printf("\nJacobi :\n");
    printf("Nombre d'itérations : %d\n", nbite);
    printf("Résidu final : %e\n", resvec[nbite-1]);
    
    // Analyse de la convergence
    printf("\nAnalyse de la convergence :\n");
    printf("Iteration |    Résidu    | Ratio de convergence\n");
    printf("-----------------------------------------\n");
    for(int i = 0; i < nbite; i += nbite/5) {  // Affiche ~5 points
        double ratio = (i > 0) ? resvec[i]/resvec[i-1] : 0.0;
        printf("%9d | %11.1e | %19.2f\n", i, resvec[i], ratio);
    }
    printf("%9d | %11.1e | %19.2f\n", nbite-1, resvec[nbite-1], 
           resvec[nbite-1]/resvec[nbite-2]);
    /* 
     // Affichage des résidus pour chaque itération
    printf("\nRésidus pour chaque itération :\n");
    for(int i = 0; i <= nbite; i++) {
        printf("%e\n", resvec[i]);
    }
    */
    // Affichage des solutions pour comparaison
    printf("\nComparaison des solutions :\n");
    printf("     X     |  Analytique  |  Jacobi  |   Différence\n");
    printf("------------------------------------------------\n");
    for(int i = 0; i < la; i++) {
        printf("%9.6f | %11.6f | %11.6f | %11.6f\n", 
               X[i], EX_SOL[i], SOL[i], fabs(EX_SOL[i] - SOL[i]));
    }
    
    // Calcul de l'erreur par rapport à la solution analytique
    relres = relative_forward_error(SOL, EX_SOL, &la);
    printf("\nErreur relative par rapport à la solution analytique : %e\n", relres);
    
    // Sauvegarde des solutions pour comparaison
    write_vec(SOL, &la, "SOL_jacobi.dat");
  }

  /* Solve with Gauss-Seidel */
  if (IMPLEM == GS) {
    write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
    
    // Résolution avec Gauss-Seidel
    gauss_seidel_tridiag(AB, RHS, SOL, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    
    printf("\nGauss-Seidel :\n");
    printf("Nombre d'itérations : %d\n", nbite);
    printf("Résidu final : %e\n", resvec[nbite-1]);
    
    // Analyse de la convergence
    printf("\nAnalyse de la convergence :\n");
    printf("Iteration |    Résidu    | Ratio de convergence\n");
    printf("-----------------------------------------\n");
    for(int i = 0; i < nbite; i += nbite/5) {  // Affiche ~5 points
        double ratio = (i > 0) ? resvec[i]/resvec[i-1] : 0.0;
        printf("%9d | %11.1e | %19.2f\n", i, resvec[i], ratio);
    }
    printf("%9d | %11.1e | %19.2f\n", nbite-1, resvec[nbite-1], 
           resvec[nbite-1]/resvec[nbite-2]);
    
    // Affichage des solutions pour comparaison
    printf("\nComparaison des solutions :\n");
    printf("     X     |  Analytique  | Gauss-Seidel |   Différence\n");
    printf("------------------------------------------------\n");
    for(int i = 0; i < la; i++) {
        printf("%9.6f | %11.6f | %11.6f | %11.6f\n", 
               X[i], EX_SOL[i], SOL[i], fabs(EX_SOL[i] - SOL[i]));
    }
    
    // Calcul de l'erreur par rapport à la solution analytique
    relres = relative_forward_error(SOL, EX_SOL, &la);
    printf("\nErreur relative par rapport à la solution analytique : %e\n", relres);
    
    // Sauvegarde des solutions pour comparaison
    write_vec(SOL, &la, "SOL_gauss_seidel.dat");
    /*
    // Sauvegarde des résidus pour le graphe
    printf("\nRésidus pour chaque itération :\n");
    for(int i = 0; i <= nbite; i++) {
        printf("%e\n", resvec[i]);
    }
    */
  }

  /* Write solution */
  write_vec(SOL, &la, "SOL.dat");

  /* Write convergence history */
  write_vec(resvec, &nbite, "RESVEC.dat");
/*
  printf("\nDonnées pour le graphe (format CSV) :\n");
  printf("Iteration,Residu\n");
  for(int i = 0; i < nbite; i++) {
      printf("%d,%e\n", i, resvec[i]);
  }
*/

  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End -----------\n");
}
