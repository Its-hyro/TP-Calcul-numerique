# Reponses aux questions

## Question 3

En C, comment doit on déclarer et allouer une matrice pour utiliser BLAS et LAPACK ?

Pour utiliser BLAS et LAPACK en langage C, les matrices doivent être déclarées et allouées sous forme de tableaux unidimensionnels, organisés en stockage par colonnes (column-major order). Une matrice de dimensions  m \times n  sera allouée dynamiquement du style :
```c
double *A;  
A = (double *)malloc(m * n * sizeof(double));  
```

L’accès à l’élément  A_{i,j}  (ligne  i , colonne  j ) se fait alors par l’index  A[i + j * m] , conformément au stockage column-major.

2. Quelle est la signification de la constante LAPACK COL MAJOR ?

La constante LAPACK COL MAJOR est une constante qui indique que la matrice est stockée par colonnes (column-major order) donc en majorité colonne dans la mémoire. Cela signifie que les éléments de la matrice sont stockés dans un tableau unidimensionnel, où les éléments de chaque colonne sont contigus. Ce qui n 'est pas etonnant vu que ça correspond au format de stockage natif utilisé par Fortran, permettant une meilleur compatibilité avec les bibliothèques LAPACK.

3. A quoi correspond la dimension principale (leading dimension) généralement notée ld ?

La dimension principale (leading dimension) est un paramètre qui
représente la taille physique en mémoire d’une colonne d’une matrice dans un tableau unidimensionnel. Elle permet de distinguer la dimension réelle de la matrice de son espace d’allocation en mémoire. Pour une matrice de taille  m \times n , le ld doit être au moins égal au nombre de lignes  m  dans le cas du stockage par colonnes. 

4. Que fait la fonction dgbmv ? Quelle m´ethode implemente-t-elle ?

La fonction dgbmv effectue un produit matrice-vecteur pour une matrice stockée en format bande général. Elle implémente l’opération :

y = \alpha A x + \beta y

A est une matrice bande, x et y sont des vecteurs, et \alpha, \beta sont des scalaires.

5. Que fait la fonction dgbtrf ? Quelle methode implemente-t-elle ?

La fonction dgbtrf effectue la factorisation LU d'une matrice bande générale avec pivotage partiel. Elle implémente l’opération :

A = PLU

A est une matrice bande, P est une matrice de permutation, L est une matrice triangulaire inférieure et U est une matrice triangulaire supérieure.

6. Que fait la fonction dgbtrs ? Quelle m´ethode impl´emente-t-elle ?

La fonction dgbtrs effectue la résolution d'un système linéaire triangulaire apres la factorisation LU effectuée par dgbtrf. Elle utilise la décomposition  A = PLU  pour résoudre :

A X = B

X est la solution recherchée et B le second membre.

7. Que fait la fonction dgbsv ? Quelle m´ethode impl´emente-t-elle ?

La fonction dgbsv effectue la résolution d'un système linéaire avec une matrice bande générale. Elle combine les étapes de factorisation LU (par dgbtrf) et de résolution (par dgbtrs) pour fournir directement la solution  X  de :

A X = B

A est une matrice bande, X et B sont des vecteurs.

8. Comment calculer la norme du r´esidu relatif avec des appels BLAS ?

Le calcul de la norme du résidu relatif consiste à évaluer l’erreur relative entre le second membre  b  et le produit de la matrice  A  par la solution approchée  x . L’expression mathématique de cette grandeur est donnée par :

||b - Ax||₂ / ||b||₂

On peut le faire avec les appels BLAS mais globalement on fait le produit matrice-vecteur avec dgbmv, la soustraction vectorielle avec daxpy, et le calcul des normes via dnrm2. 

Calcule du produit matrice-vecteur Ax avec `dgbmv` :
```c
dgbmv("N", &la, &la, &kl, &ku, &alpha, AB, &lab, x, &incx, &beta, y, &incy);
```

Calcule du résidu b - Ax avec `daxpy` :
```c
dcopy(&la, b, &incx, r, &incy);
alpha = -1.0;
daxpy(&la, &alpha, y, &incx, r, &incy);
```

Calcule des normes avec `dnrm2` et on fait la division :
```c
double norm_r = dnrm2(&la, r, &incx);
double norm_b = dnrm2(&la, b, &incx);
double rel_res = norm_r / norm_b;
```


## Question 5

Evaluer les performances. Que dire de la complexité des méthodes appelées ?

Les résultats expérimentaux montrent les performances suivantes pour différentes tailles de matrices :
- n = 100 : 0.000094s
- n = 1000 : 0.000019s
- n = 10000 : 0.000194s

Analyse de la complexité :

1. Complexité théorique
- Factorisation LU (DGBTRF) : O(n) opérations
- Résolution (DGBTRS) : O(n) opérations
- Solution complète (DGBSV) : O(n) opérations

2. Validation expérimentale
- Le temps d'exécution croît quasi-linéairement avec n
- Ratio temps(n=10000)/temps(n=1000) ≈ 10.2
- Confirme la complexité théorique O(n)

3. Efficacité
Cette performance s'explique par une exploitation de la structure bande (3 diagonales) et un stockage optimisé en mémoire.

Cette complexité linéaire O(n) est optimale pour ce type de problème, comparée aux méthodes générales :
- LU dense : O(n³)
- QR dense : O(n³)

L'implémentation est donc parfaitement adaptée aux systèmes tridiagonaux, offrant une résolution rapide et stable même pour de grandes matrices. La complexité des méthodes appelées est en accord avec la complexité théorique.

