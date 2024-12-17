/******************************************/
/* tp_poisson1D_direct.c                 */
/* Programme principal pour résoudre      */
/* l'équation de Poisson 1D par          */
/* méthode directe (LU)                  */
/******************************************/
#include "lib_poisson1D.h"
#include <time.h>
#include <string.h>

int main(int argc,char *argv[])
{
  int ierr;
  int jj;
  int nbpoints, la;        // Nombre de points et taille du système
  int ku, kl, kv, lab;     // Paramètres de la matrice bande
  int *ipiv;               // Vecteur de permutation
  int info;                // Code de retour de notre factorisation
  int info2;              // Code de retour de LAPACK (pour comparaison)
  int NRHS = 1;           // Nombre de seconds membres
  double T0 = -5.0;       // Condition limite gauche
  double T1 = 5.0;        // Condition limite droite
  double *RHS, *EX_SOL, *X;  // Second membre, solution exacte et grille
  double *AB, *AB_copy;    // Matrice bande et sa copie
  clock_t start, end;      // Pour mesurer le temps
  double cpu_time_used;

  // Initialisation des paramètres
  kv = 1;                  // Nombre de sur/sous-diagonales
  ku = 1;                  // Nombre de sur-diagonales
  kl = 1;                  // Nombre de sous-diagonales
  nbpoints = 10;           // Nombre de points de la grille
  la = nbpoints-2;         // Taille du système (sans les bords)
  lab = kv+kl+ku+1;        // Largeur de bande

  // Allocation de la mémoire
  RHS = (double *) malloc(sizeof(double)*la);      // Second membre
  EX_SOL = (double *) malloc(sizeof(double)*la);   // Solution exacte
  X = (double *) malloc(sizeof(double)*la);        // Points de la grille
  AB = (double *) malloc(sizeof(double)*lab*la);   // Matrice bande
  AB_copy = (double *) malloc(sizeof(double)*lab*la);  // Copie pour comparaison
  ipiv = (int *) malloc(la * sizeof(int));         // Vecteur de permutation

  // Configuration du problème
  set_grid_points_1D(X, &la);                      // Définition de la grille
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);          // Second membre avec conditions aux limites
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);  // Construction de la matrice

  // Sauvegarde d'une copie de la matrice pour comparaison
  memcpy(AB_copy, AB, sizeof(double)*lab*la);

  printf("\n=== Test de factorisation LU ===\n");
  
  // 1. Factorisation LU avec notre implémentation
  printf("\n1. Factorisation avec dgbtrftridiag :\n");
  info = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  printf("Info dgbtrftridiag = %d\n", info);
  
  if (info == 0) {  // Si la factorisation a réussi
    // Validation de la factorisation LU
    if (!validate_LU_tridiag(AB, &lab, &la, &kv)) {
      printf("La validation de la factorisation LU a échoué\n");
      return 1;
    }
    
    // 2. Test de la solution
    printf("\n2. Test de la solution :\n");
    double *RHS_copy = (double *) malloc(sizeof(double)*la);
    memcpy(RHS_copy, RHS, sizeof(double)*la);
    
    // Résolution avec notre factorisation
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    
    // Résolution avec LAPACK pour comparaison
    dgbsv_(&la, &kl, &ku, &NRHS, AB_copy, &lab, ipiv, RHS_copy, &la, &info2);
    
    // Comparaison des deux solutions
    double max_diff = 0.0;
    for(int i = 0; i < la; i++) {
      double diff = fabs(RHS[i] - RHS_copy[i]);
      if (diff > max_diff) max_diff = diff;
    }
    printf("Différence maximale entre les solutions : %e\n", max_diff);
    
    // Calcul de l'erreur par rapport à la solution exacte
    double rel_err = relative_forward_error(RHS, EX_SOL, &la);
    printf("Erreur relative par rapport à la solution exacte : %e\n", rel_err);
    
    free(RHS_copy);
  } else {
    printf("La factorisation a échoué avec info = %d\n", info);
  }

  // Libération de la mémoire
  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(AB_copy);
  free(ipiv);

  printf("\n=====================================\n");
  return 0;
}
