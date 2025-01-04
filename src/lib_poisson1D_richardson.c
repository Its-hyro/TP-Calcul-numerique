/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){
    // Pour une matrice tridiagonale de Poisson 1D
    // λmax = 4*sin²(π/(2(n+1)))
    return 4.0 * pow(sin(M_PI/(2.0*(*la + 1))), 2);
}

double eigmin_poisson1D(int *la){
    // Pour une matrice tridiagonale de Poisson 1D
    // λmin = 4*sin²(nπ/(2(n+1)))
    return 4.0 * pow(sin(*la * M_PI/(2.0*(*la + 1))), 2);
}

double richardson_alpha_opt(int *la){
    // Le alpha optimal pour la méthode de Richardson est :
    // α = 2/(λmin + λmax)
    double lambda_max = eigmax_poisson1D(la);
    double lambda_min = eigmin_poisson1D(la);
    double alpha = 2.0/(lambda_max + lambda_min);
    
    printf("\nValeurs propres et alpha :\n");
    printf("lambda_max = %e\n", lambda_max);
    printf("lambda_min = %e\n", lambda_min);
    printf("alpha_opt calculé = %e\n", alpha);
    
    return alpha;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite){
    int i;
    int kv = 1;
    double *AX = (double *) malloc(sizeof(double)*(*la));
    double *resid = (double *) malloc(sizeof(double)*(*la));
    
    // Debug: Affichage de la matrice AB
    printf("\nMatrice AB au format GB :\n");
    for (i = 0; i < *lab; i++) {
        for (int j = 0; j < *la; j++) {
            printf("%f ", AB[j * (*lab) + i]);
        }
        printf("\n");
    }
    
    // Debug: Affichage du vecteur RHS
    printf("\nVecteur RHS :\n");
    for(i = 0; i < *la; i++) {
        printf("%f ", RHS[i]);
    }
    printf("\n");
    
    // Initialisation
    *nbite = 0;
    
    do {
        // 1. Calcul de AX
        dgbmv_poisson1D(AB, X, AX, la, lab, ku, kl, &kv);
        
        // 2. Calcul du résidu r = RHS - AX
        double norm_rhs = 0.0;
        for(i = 0; i < *la; i++){
            resid[i] = RHS[i] - AX[i];
            norm_rhs += RHS[i] * RHS[i];
        }
        
        // 3. Calcul de la norme du résidu normalisé
        double norm_res = 0.0;
        for(i = 0; i < *la; i++){
            norm_res += resid[i] * resid[i];
        }
        norm_res = sqrt(norm_res/norm_rhs);
        
        // Debug: Afficher tous les 100 itérations
        if(*nbite % 100 == 0) {
            printf("Iteration %d: résidu = %e\n", *nbite, norm_res);
        }
        
        // 4. Sauvegarde de la norme du résidu
        resvec[*nbite] = norm_res;
        
        // 5. Mise à jour de la solution : X = X + alpha * resid
        for(i = 0; i < *la; i++){
            X[i] = X[i] + (*alpha_rich) * resid[i];
        }
        
        (*nbite)++;
        
    } while (*nbite < *maxit && resvec[*nbite-1] > *tol);
    
    free(AX);
    free(resid);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

