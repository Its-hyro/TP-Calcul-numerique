/**********************************************/
/* lib_poisson1D.c                            */
/* Bibliothèque numérique pour résoudre       */
/* l'équation de Poisson 1D                   */
/* (Équation de la chaleur stationnaire)      */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv) {
    int ii, jj, kk;
    
    // Initialisation de la matrice à zéro
    for (ii = 0; ii < (*lab) * (*la); ii++) {
        AB[ii] = 0.0;
    }
    
    // Construction de la matrice tridiagonale
    for (ii = 0; ii < *la; ii++) {
        // Diagonale principale : -2
        AB[(*kv) + ii * (*lab)] = -2.0;
        
        // Sur-diagonale : 1 (sauf dernière colonne)
        if (ii < *la - 1) {
            AB[(*kv-1) + ii * (*lab)] = 1.0;
        }
        
        // Sous-diagonale : 1 (sauf première colonne)
        if (ii > 0) {
            AB[(*kv+1) + ii * (*lab)] = 1.0;
        }
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info) {
    int i;
    double pivot;
    *info = 0;

    // Vérification des paramètres (matrice tridiagonale uniquement)
    if (*kl != 1 || *ku != 1) {
        *info = -1;
        return *info;
    }

    // Structure de la matrice bande AB :
    // - AB[0 + i*(*lab)] : sur-diagonale (u_i)
    // - AB[1 + i*(*lab)] : diagonale principale (d_i)
    // - AB[2 + i*(*lab)] : sous-diagonale (l_i)

    // Factorisation LU sans pivotage
    for (i = 0; i < *n-1; i++) {
        // Vérification du pivot
        pivot = AB[1 + i*(*lab)];
        if (fabs(pivot) < 1e-12) {
            *info = i+1;  // Échec : pivot trop petit
            return *info;
        }

        // Calcul du multiplicateur L[i+1,i]
        AB[2 + i*(*lab)] = AB[2 + i*(*lab)] / pivot;

        // Mise à jour de l'élément diagonal suivant
        AB[1 + (i+1)*(*lab)] = AB[1 + (i+1)*(*lab)] - 
                               AB[2 + i*(*lab)] * AB[0 + i*(*lab)];
    }

    // Vérification du dernier pivot
    if (fabs(AB[1 + (*n-1)*(*lab)]) < 1e-12) {
        *info = *n;
        return *info;
    }

    // Pas de permutation (matrice tridiagonale)
    for (i = 0; i < *n; i++) {
        ipiv[i] = i+1;
    }

    return *info;
}

int validate_GB_operator_poisson1D(double* AB, int *lab, int *la, int *kv) {
    int ii;
    double eps = 1e-12;  // Tolérance pour les comparaisons
    int valid = 1;       // Flag de validation

    printf("=== Validation de la matrice de Poisson 1D ===\n");
    
    // Test 1 : Vérification de la diagonale principale (-2)
    printf("Test 1: Vérification de la diagonale principale...\n");
    for (ii = 0; ii < *la; ii++) {
        if (fabs(AB[(*kv) + ii * (*lab)] + 2.0) > eps) {
            printf("  Erreur: Diagonale principale incorrecte à la position %d: %f != -2.0\n", 
                   ii, AB[(*kv) + ii * (*lab)]);
            valid = 0;
        }
    }
    if (valid) printf("  OK: Diagonale principale correcte\n");

    // Test 2 : Vérification des diagonales secondaires (1)
    printf("\nTest 2: Vérification des diagonales secondaires...\n");
    for (ii = 0; ii < *la-1; ii++) {
        // Vérification de la sur-diagonale
        if (fabs(AB[(*kv-1) + ii * (*lab)] - 1.0) > eps) {
            printf("  Erreur: Sur-diagonale incorrecte à la position %d: %f != 1.0\n", 
                   ii, AB[(*kv-1) + ii * (*lab)]);
            valid = 0;
        }
        // Vérification de la sous-diagonale
        if (fabs(AB[(*kv+1) + (ii+1) * (*lab)] - 1.0) > eps) {
            printf("  Erreur: Sous-diagonale incorrecte à la position %d: %f != 1.0\n", 
                   ii, AB[(*kv+1) + (ii+1) * (*lab)]);
            valid = 0;
        }
    }
    if (valid) printf("  OK: Diagonales secondaires correctes\n");

    // Test 3 : Validation par produit matrice-vecteur
    printf("\nTest 3: Vérification du produit matrice-vecteur...\n");
    double *x_test = (double *) malloc(sizeof(double)*(*la));
    double *y_test = (double *) malloc(sizeof(double)*(*la));
    
    // Initialisation du vecteur test (tous les éléments à 1)
    for(ii = 0; ii < *la; ii++) {
        x_test[ii] = 1.0;
        y_test[ii] = 0.0;
    }
    
    // Calcul de y = Ax avec BLAS
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1;
    int incy = 1;
    
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, 1, 1, 
                alpha, AB, *lab, x_test, incx, beta, y_test, incy);
    
    // Vérification des propriétés du produit
    int valid_mv = 1;
    printf("  Vérification des valeurs :\n");
    for(ii = 0; ii < *la; ii++) {
        printf("    Position %d: OK (valeur: %f)\n", ii, y_test[ii]);
        
        // Vérification des propriétés attendues
        if (ii == 0 || ii == 1 || ii == *la-2 || ii == *la-1) {
            // Points aux bords : devrait être -1
            if (fabs(y_test[ii] + 1.0) > eps) {
                printf("      ATTENTION: Devrait être -1 (point de bord)\n");
                valid_mv = 0;
                valid = 0;
            }
        } else {
            // Points intérieurs : devrait être 0
            if (fabs(y_test[ii]) > eps) {
                printf("      ATTENTION: Devrait être 0 (point intérieur)\n");
                valid_mv = 0;
                valid = 0;
            }
        }
    }
    if (valid_mv) printf("  OK: Toutes les valeurs sont correctes\n");
    
    free(x_test);
    free(y_test);

    printf("\nRésultat final: %s\n", valid ? "SUCCÈS" : "ÉCHEC");
    printf("=====================================\n");
    
    return valid;
}

int validate_LU_tridiag(double *AB, int *lab, int *la, int *kv) {
    int ii;
    double eps = 1e-12;  // Tolérance pour les comparaisons
    int valid = 1;       // Flag de validation
    
    printf("\n=== Validation de la factorisation LU tridiagonale ===\n");
    
    // Test 1 : Vérification des pivots
    printf("1. Vérification des pivots :\n");
    for(ii = 0; ii < *la; ii++) {
        if(fabs(AB[1 + ii*(*lab)]) < eps) {
            printf("  ERREUR: Pivot nul ou trop petit à la position %d\n", ii);
            valid = 0;
        }
    }
    if(valid) printf("  OK: Tous les pivots sont non nuls\n");
    
    // Test 2 : Vérification des multiplicateurs de L
    printf("\n2. Vérification des multiplicateurs :\n");
    for(ii = 0; ii < *la-1; ii++) {
        if(fabs(AB[2 + ii*(*lab)]) > 1.0 + eps) {
            printf("  ERREUR: Multiplicateur > 1 à la position %d\n", ii);
            valid = 0;
        }
    }
    if(valid) printf("  OK: Tous les multiplicateurs sont <= 1\n");
    
    // Test 3 : Vérification de la structure LU
    printf("\n3. Test de la décomposition LU :\n");
    printf("  Vérification des éléments de la matrice LU :\n");
    
    // Vérification du premier élément diagonal
    if(fabs(AB[1]) != 2.0) {
        printf("  ERREUR: Premier élément diagonal incorrect : %f != 2.0\n", fabs(AB[1]));
        valid = 0;
    }
    
    // Vérification de la structure L et U
    for(ii = 1; ii < *la; ii++) {
        double u_ii = AB[1 + ii*(*lab)];      // Diagonal de U
        double l_i = AB[2 + (ii-1)*(*lab)];   // Sous-diagonal de L
        double u_i = (ii < *la-1) ? AB[0 + ii*(*lab)] : 0.0;  // Sur-diagonal de U
        
        printf("  Col %d: l = %f, u = %f, d = %f\n", 
               ii, l_i, u_i, u_ii);
        
        // Vérification des propriétés de L
        if(fabs(l_i) > 1.0 + eps) {
            printf("  ERREUR: Multiplicateur L trop grand à la position %d\n", ii);
            valid = 0;
        }
        
        // Vérification des propriétés de U
        if(ii < *la-1 && fabs(u_i - 1.0) > eps) {
            printf("  ERREUR: Sur-diagonale U incorrecte à la position %d\n", ii);
            valid = 0;
        }
    }
    
    printf("\nRésultat final de la validation : %s\n", valid ? "SUCCÈS" : "ÉCHEC");
    printf("==========================================\n");
    
    return valid;
}
