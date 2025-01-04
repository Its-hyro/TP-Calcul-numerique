/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, int* kv) {
    int filler = *kv;
    int nb_cols = *la;
    int nb_lines = *lab;

    printf("\nDébut set_GB_operator_colMajor_poisson1D:\n");
    printf("filler (kv) = %d\n", filler);
    printf("nb_cols (la) = %d\n", nb_cols);
    printf("nb_lines (lab) = %d\n", nb_lines);

    for (int i = 0; i < nb_cols * nb_lines; i += nb_lines) {
        for (int j = 0; j < filler; j++) {
            AB[i + j] = 0;
        }
        AB[i + filler] = -1;
        AB[i + filler + 1] = 2;
        AB[i + filler + 2] = -1;
    }

    // Make the final elem 0
    AB[nb_cols * nb_lines - 1] = 0;
    // Make the first elem of the first non kv line 0
    AB[filler] = 0;

    printf("\nMatrice construite (format bande):\n");
    for(int i = 0; i < nb_lines; i++) {
        for(int j = 0; j < nb_cols; j++) {
            printf("%6.2f ", AB[j * nb_lines + i]);
        }
        printf("\n");
    }
}

void set_GB_operator_colMajor_poisson1D_DGBMV(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  // Initialisation à 0
  for (ii = 0; ii < (*lab * *la); ii++) {
    AB[ii] = 0.0;
  }
  
  // Format pour DGBMV
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (jj > 0) AB[kk+0]=-1.0;      // sous-diagonale
    AB[kk+1]=2.0;                    // diagonale
    if (jj < (*la)-1) AB[kk+2]=-1.0; // sur-diagonale
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
    double norm_diff = 0.0;
    double norm_y = 0.0;
    int i;
    
    // Calcul de ||x-y|| et ||y||
    for(i = 0; i < *la; i++){
        norm_diff += (x[i] - y[i]) * (x[i] - y[i]);
        norm_y += y[i] * y[i];
    }
    
    // Éviter la division par zéro
    if(norm_y < 1e-15) return 0.0;
    
    // Calcul de l'erreur relative ||x-y|| / ||y||
    return sqrt(norm_diff) / sqrt(norm_y);
}

int indexABCol(int i, int j, int* lab) {
    return (j + 1) * (*lab - 1) + i - 1;
}

int dgbtrftridiag(int* la, int* n, int* kl, int* ku, double* AB, int* lab, int* ipiv, int* info) {
    int nb_lines = *lab;
    int nb_cols = *la;
    int upper_diags = *ku;
    int lower_diags = *kl;
    int result_matrix_size = *n;
    int kv = 1;  // Position de la diagonale

    printf("\nDébut dgbtrftridiag:\n");
    printf("nb_lines (lab) = %d\n", nb_lines);
    printf("nb_cols (la) = %d\n", nb_cols);
    printf("upper_diags (ku) = %d\n", upper_diags);
    printf("lower_diags (kl) = %d\n", lower_diags);
    printf("result_matrix_size (n) = %d\n", result_matrix_size);

    // Vérification des dimensions
    if (nb_lines != 4 || lower_diags != 1 || upper_diags != 1 || nb_cols < result_matrix_size) {
        printf("Erreur: dimensions incorrectes\n");
        *info = -1;
        return *info;
    }

    // Pour une matrice tridiagonale au format bande:
    // - ligne 0: zéros (non utilisé)
    // - ligne 1: sous-diagonale (l)
    // - ligne 2: diagonale (u)
    // - ligne 3: sur-diagonale

    // Factorisation LU
    for (int k = 0; k < result_matrix_size - 1; k++) {
        // Vérification du pivot
        if(fabs(AB[k * nb_lines + 2]) < 1e-10) {
            *info = k + 1;
            return *info;
        }

        // Calcul du multiplicateur
        double l = AB[k * nb_lines + 3] / AB[k * nb_lines + 2];
        
        // Stockage du multiplicateur dans la sous-diagonale
        AB[(k+1) * nb_lines + 1] = l;
        
        // Mise à jour de la diagonale suivante
        AB[(k+1) * nb_lines + 2] = AB[(k+1) * nb_lines + 2] - l * AB[k * nb_lines + 3];
    }

    // Vérification du dernier pivot
    if(fabs(AB[(result_matrix_size - 1) * nb_lines + 2]) < 1e-10) {
        *info = result_matrix_size;
        return *info;
    }

    // Mise à jour du vecteur ipiv
    for(int i = 0; i < result_matrix_size; i++) {
        ipiv[i] = i + 1;
    }

    *info = 0;
    return *info;
}

void dgbmv_poisson1D(double *AB, double *RHS, double *X, int *la, int *lab, int *ku, int *kl, int *kv) {
  // Pour une matrice tridiagonale :
  // kl = 1 (une sous-diagonale)
  // ku = 1 (une sur-diagonale)
  
  double alpha = 1.0;
  double beta = 0.0;
  int incx = 1;
  int incy = 1;
  
  // Appel à dgbmv avec les bons paramètres
  // - *la est le nombre de lignes et colonnes
  // - 1 est le nombre de sous-diagonales (kl)
  // - 1 est le nombre de sur-diagonales (ku)
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, 1, 1, alpha, 
              AB, *lab, RHS, incx, beta, X, incy);
}

int test_dgbmv_poisson1D(void) {
  // Taille du problème
  int la = 5;  // Taille de la matrice
  int ku = 1;  // Nombre de sur-diagonales
  int kl = 1;  // Nombre de sous-diagonales
  int kv = 1;  // Position de la diagonale
  int lab = kv + kl + ku + 1;
  
  // Allocation de mémoire
  double *AB = (double *)malloc(sizeof(double)*lab*la);
  double *RHS = (double *)malloc(sizeof(double)*la);
  double *X = (double *)malloc(sizeof(double)*la);
  double *X_exact = (double *)malloc(sizeof(double)*la);
  
  // Initialisation de la matrice
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  // Initialisation du vecteur RHS avec des 1
  for(int i = 0; i < la; i++) {
    RHS[i] = 1.0;
  }
  
  // Calcul avec dgbmv
  dgbmv_poisson1D(AB, RHS, X, &la, &lab, &ku, &kl, &kv);
  
  // Vérification : pour une matrice tridiagonale [-1 2 -1] et RHS = [1,...,1]
  // La solution exacte peut être calculée analytiquement
  int valid = 1;
  for(int i = 0; i < la; i++) {
    // Vérification que le résultat est cohérent (à une tolérance près)
    if(fabs(X[i] - RHS[i]*2.0 + (i > 0 ? RHS[i-1]*(-1.0) : 0.0) + 
           (i < la-1 ? RHS[i+1]*(-1.0) : 0.0)) > 1e-10) {
      valid = 0;
      break;
    }
  }
  
  // Libération de la mémoire
  free(AB);
  free(RHS);
  free(X);
  free(X_exact);
  
  return valid;
}

int test_dgbtrftridiag(void) {
    int la = 5;  // Taille de la matrice
    int ku = 1;  // Nombre de sur-diagonales
    int kl = 1;  // Nombre de sous-diagonales
    int kv = 1;  // Position de la diagonale
    int lab = 4; // kl + ku + kv + 1
    int *ipiv = (int *)malloc(sizeof(int)*la);
    int info;
    
    // Allocation et initialisation de la matrice
    double *AB = (double *)malloc(sizeof(double)*lab*la);
    double *AB_copy = (double *)malloc(sizeof(double)*lab*la);
    
    // Création de la matrice tridiagonale
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    
    // Copie de la matrice originale
    for(int i = 0; i < lab*la; i++) {
        AB_copy[i] = AB[i];
    }
    
    // Factorisation LU
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    
    // Vérification : les éléments de L et U doivent être cohérents
    int valid = 1;
    if(info == 0) {
        // Vérification des éléments diagonaux de U
        for(int j = 0; j < la; j++) {
            double u_diag = AB[j*lab + kv + 1];
            if(fabs(u_diag) < 1e-10) {
                valid = 0;
                printf("Pivot nul trouvé à la position %d\n", j);
                break;
            }
        }
        
        // Vérification que la factorisation est correcte en reconstruisant A = LU
        if(valid) {
            // Pour chaque ligne
            for(int i = 0; i < la; i++) {
                // Pour chaque colonne
                for(int j = 0; j < la; j++) {
                    if(abs(i-j) <= 1) {  // On ne vérifie que la tridiagonale
                        double sum = 0.0;
                        // Calcul de (LU)_{ij}
                        for(int k = 0; k <= i && k <= j; k++) {
                            double l_ik = (i == k) ? 1.0 : 
                                        (i == k+1) ? AB[(k+1)*lab + kv] : 0.0;
                            double u_kj = (k == j) ? AB[k*lab + kv + 1] :
                                        (k == j-1) ? AB[k*lab + kv + 2] : 0.0;
                            sum += l_ik * u_kj;
                        }
                        // Comparaison avec la matrice originale
                        double orig;
                        if(i == j) orig = 2.0;
                        else if(abs(i-j) == 1) orig = -1.0;
                        else orig = 0.0;
                        
                        if(fabs(sum - orig) > 1e-10) {
                            valid = 0;
                            printf("Erreur de reconstruction à la position (%d,%d): %f != %f\n", 
                                   i, j, sum, orig);
                            break;
                        }
                    }
                }
                if(!valid) break;
            }
        }
    } else {
        valid = 0;
        printf("La factorisation a échoué avec info = %d\n", info);
    }
    
    // Libération de la mémoire
    free(AB);
    free(AB_copy);
    free(ipiv);
    
    return valid;
}

void jacobi_tridiag(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    // Allocation des vecteurs temporaires
    double *X_new = (double *)malloc(sizeof(double)*(*la));
    double *AX = (double *)malloc(sizeof(double)*(*la));
    int kv = 1;  // Position de la diagonale dans le format GB
    
    int iter = 0;
    double resid = 1.0;
    resvec[0] = 1.0;
    
    while(iter < *maxit && resid > *tol) {
        // Calcul de AX = AB*X
        dgbmv_poisson1D(AB, X, AX, la, lab, ku, kl, &kv);
        
        // Mise à jour de X selon la méthode de Jacobi
        for(int i = 0; i < *la; i++) {
            // Pour une matrice tridiagonale au format GB :
            // AB[(*lab)*i + 0] : sous-diagonale
            // AB[(*lab)*i + 1] : diagonale
            // AB[(*lab)*i + 2] : sur-diagonale
            double diag = AB[(*lab)*i + 1];
            X_new[i] = RHS[i];
            
            if(i > 0) X_new[i] -= AB[(*lab)*i + 0] * X[i-1];
            if(i < *la-1) X_new[i] -= AB[(*lab)*i + 2] * X[i+1];
            
            X_new[i] /= diag;
        }
        
        // Calcul du résidu
        resid = 0.0;
        for(int i = 0; i < *la; i++) {
            double diff = X_new[i] - X[i];
            resid += diff * diff;
        }
        resid = sqrt(resid);
        
        // Mise à jour de X pour la prochaine itération
        for(int i = 0; i < *la; i++) {
            X[i] = X_new[i];
        }
        
        iter++;
        resvec[iter] = resid;
        
        if(iter % 100 == 0) {
            printf("Iteration %d: résidu = %e\n", iter, resid);
        }
    }
    
    *nbite = iter;
    
    // Libération de la mémoire
    free(X_new);
    free(AX);
}

void gauss_seidel_tridiag(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    // Allocation des vecteurs temporaires
    double *AX = (double *)malloc(sizeof(double)*(*la));
    int kv = 1;
    
    int iter = 0;
    double resid = 1.0;
    resvec[0] = 1.0;
    
    while(iter < *maxit && resid > *tol) {
        // Calcul de AX = AB*X pour le résidu
        dgbmv_poisson1D(AB, X, AX, la, lab, ku, kl, &kv);
        
        // Mise à jour de X selon Gauss-Seidel
        for(int i = 0; i < *la; i++) {
            double sum = RHS[i];  // bi
            
            // Soustraction des termes déjà calculés (partie E)
            if(i > 0) {
                sum -= AB[(*lab)*i + 0] * X[i-1];  // Utilise la nouvelle valeur X[i-1]
            }
            
            // Soustraction des termes non encore calculés (partie F)
            if(i < *la-1) {
                sum -= AB[(*lab)*i + 2] * X[i+1];  // Utilise l'ancienne valeur X[i+1]
            }
            
            // Division par l'élément diagonal
            X[i] = sum / AB[(*lab)*i + 1];
        }
        
        // Calcul du résidu
        resid = 0.0;
        dgbmv_poisson1D(AB, X, AX, la, lab, ku, kl, &kv);
        for(int i = 0; i < *la; i++) {
            double diff = RHS[i] - AX[i];
            resid += diff * diff;
        }
        resid = sqrt(resid);
        
        iter++;
        resvec[iter] = resid;
        
        if(iter % 100 == 0) {
            printf("Iteration %d: résidu = %e\n", iter, resid);
        }
    }
    
    *nbite = iter;
    
    // Libération de la mémoire
    free(AX);
}
