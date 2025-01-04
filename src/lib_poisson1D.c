/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  // Initialisation à 0
  for (ii = 0; ii < (*lab * *la); ii++) {
    AB[ii] = 0.0;
  }
  
  // Format pour DGBSV
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    AB[kk+0]=-1.0;  // sous-diagonale
    AB[kk+1]=2.0;   // diagonale
    AB[kk+2]=-1.0;  // sur-diagonale
  }
  // Conditions aux limites
  AB[0]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
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

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // Vérification des paramètres
  if (*kl != 1 || *ku != 1) {
    *info = -1;
    return *info;
  }

  *info = 0;
  
  // Pour une matrice tridiagonale au format GB :
  // AB[i*(*lab) + 0] : sous-diagonale (l)
  // AB[i*(*lab) + 1] : diagonale (u)
  // AB[i*(*lab) + 2] : sur-diagonale (déjà en place)
  
  for(int i = 0; i < *n-1; i++) {
    // Si pivot nul, erreur
    if(fabs(AB[i*(*lab) + 1]) < 1e-10) {
      *info = i+1;
      return *info;
    }
    
    // Calcul du multiplicateur l = a_{i+1,i}/a_{i,i}
    double l = AB[(i+1)*(*lab) + 0] / AB[i*(*lab) + 1];
    
    // Mise à jour de la sous-diagonale avec l
    AB[(i+1)*(*lab) + 0] = l;
    
    // Mise à jour de la diagonale : a_{i+1,i+1} = a_{i+1,i+1} - l*a_{i,i+1}
    AB[(i+1)*(*lab) + 1] = AB[(i+1)*(*lab) + 1] - l * AB[i*(*lab) + 2];
    
    // Stockage des permutations (pas de permutation pour une matrice tridiagonale)
    ipiv[i] = i+1;
  }
  ipiv[*n-1] = *n;
  
  // Vérification du dernier pivot
  if(fabs(AB[(*n-1)*(*lab) + 1]) < 1e-10) {
    *info = *n;
  }
  
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
  int lab = kv + kl + ku + 1;
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
  
  // Vérification : reconstruction de la matrice originale A = L*U
  int valid = 1;
  if(info == 0) {
    // Vérification des éléments L et U
    for(int j = 0; j < la; j++) {
      // Vérification de la diagonale de U
      double u_diag = AB[j*lab + 1];
      if(fabs(u_diag) < 1e-10) {
        valid = 0;
        printf("Pivot nul trouvé à la position %d\n", j);
        break;
      }
      
      // Vérification des multiplicateurs L (doivent être < 1 en valeur absolue)
      if(j < la-1) {
        double l_subdiag = AB[(j+1)*lab + 0];
        if(fabs(l_subdiag) >= 1.0) {
          valid = 0;
          printf("Multiplicateur L invalide à la position %d: %f\n", j+1, l_subdiag);
          break;
        }
      }
      
      // Vérification de la sur-diagonale de U
      if(j < la-1) {
        double u_superdiag = AB[j*lab + 2];
        if(fabs(u_superdiag + 1.0) > 1e-10) {
          valid = 0;
          printf("Sur-diagonale U incorrecte à la position %d: %f\n", j, u_superdiag);
          break;
        }
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
