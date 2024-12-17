/**********************************************/
/* lib_poisson1D.h                            */
/* Fichier d'en-tête de la bibliothèque pour  */
/* résoudre l'équation de Poisson 1D          */
/* (Équation de la chaleur stationnaire)      */
/**********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "atlas_headers.h"

// Fonctions de construction et manipulation des matrices
void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int *la, int *kv);  // Construit la matrice de Poisson en format bande
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int *la, int *kv);  // Construit la matrice identité en format bande
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);  // Définit le second membre avec conditions aux limites
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);  // Calcule la solution analytique
void set_grid_points_1D(double* x, int* la);  // Définit les points de la grille

// Fonctions d'analyse et de mesure
double relative_forward_error(double* x, double* y, int* la);  // Calcule l'erreur relative entre deux vecteurs
int indexABCol(int i, int j, int *lab);  // Calcule l'index dans le format bande

// Fonctions de résolution directe
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info);  // Factorisation LU tridiagonale

// Fonctions de validation
int validate_GB_operator_poisson1D(double* AB, int *lab, int *la, int *kv);  // Valide la construction de la matrice
int validate_LU_tridiag(double *AB, int *lab, int *la, int *kv);  // Valide la factorisation LU

// Fonctions d'écriture et de sauvegarde
void write_GB2AIJ_operator_poisson1D(double* AB, int *la, char* filename);  // Écrit la matrice au format AIJ
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int *la, char* filename);  // Écrit la matrice en format ligne
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);  // Écrit la matrice en format colonne
void write_vec(double* vec, int* la, char* filename);  // Écrit un vecteur dans un fichier
void write_xy(double* vec, double* x, int* la, char* filename);  // Écrit les coordonnées (x,y) dans un fichier

// Fonctions pour la méthode de Richardson
void eig_poisson1D(double* eigval, int *la);  // Calcule les valeurs propres
double eigmax_poisson1D(int *la);  // Calcule la plus grande valeur propre
double eigmin_poisson1D(int *la);  // Calcule la plus petite valeur propre
double richardson_alpha_opt(int *la);  // Calcule le paramètre optimal pour Richardson
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,
                     int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite);  // Méthode de Richardson

// Fonctions pour les méthodes itératives
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv);  // Extrait la matrice de Jacobi
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv);  // Extrait la matrice de Gauss-Seidel
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
                  int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite);  // Richardson avec préconditionneur
