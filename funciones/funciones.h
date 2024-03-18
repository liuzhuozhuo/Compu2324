#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//2 pointer for the two vectors, n is the dimension of the vectors, the result is given in a1(first vector).
void suma_vec (double* a1, double* a2, int n){
    int i;
    for (i = 0; i <n; i++){
        *(a1+i) = *(a1+i)+*(a2+i);
    }
}

// generate a random number in the position i, j of the matrix of size n, m
void rand_matrix (int n, double *matrix, int i, int j){
    srand(time(NULL));
    matrix[i*n + j]= (double)rand() / RAND_MAX;
}

//return a real random number
double random_number(){ 
    double n;
    srand(time(NULL));
    n = (double)rand() / RAND_MAX;
    return n;
}