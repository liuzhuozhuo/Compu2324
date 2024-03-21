#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Update the acceleration vector, 
//r: matrix (pointer) with the position of all the planets in 3d
//a: matrix (pointer) with the acceleration of all the planets in 3d
//m: vector (pointer) with the modified mass of the planets
//N is the total number of planets
void acceleration (double* r, double* a, double* m, int N) {
    int i, j, k;
    double acc; //Dummy variable to sum all the acceleration for a planet in one component
    double R[N][N][3]; 
    double R_mod[N][N]; //Matrix NxN with the module of the distance between planets

    //Calculate the distance vector between planets. 
    for(i = 0; i < N-1; i++){
        for (j = i+1; j < N; j++){
            R_mod[i][j]= 0;
            for (k = 0; k < 3; k++){
                R[i][j][k] = *(r+i*N+k) - *(r+j*N+k); 
                R_mod[i][j] += pow(R[i][j][k], 2);
            }
        }
    }
    //Calculate the acceleration
    //iterate in the planets
     for(i = 0; i < N; i++){
        //iterate in the components
        for (k = 0; k < 3; k++){
            //iterate between planets
            acc = 0;
            for (j = 0; j < N; j++){
                if(j>i){
                    acc += (*(m+j))*R[i][j][k]/pow(R_mod[i][j], 3./2.);
                }
                if(j<i){
                    acc += -(*(m+j))*R[j][i][k]/pow(R_mod[j][i], 3./2.);
                }
            } 
            *(a+i*N+k) = acc;
        } 
    }    
}

//Function that calculate the next set of position and velocity using the verlet algorithms
void verlet_algorithm(double* r, double* v ,double* a, double* m, int N, double h){
    int i, j, k;
    double w[N][3];

    //Calculate the new positions
    for (i = 0; i < N; i++){
        for (k = 0; k < 3; k++){
            *(r+i*N+k) += h*(*(v+i*N+k)) + (pow(h, 2)/2.) * (*(a+i*N+k));
            w[N][k] = *(v+i*N+k) + (h/2.)*(*(a+i*N+k));
        }
    }
   //Using the acceleration function calculate the net acceleration for each of the planets
    acceleration(r, a, m, N);

    //Calculate the new velocity of the planets
    for (i = 0; i < N; i++){
        for (k = 0; k < 3; k++){
            *(v+i*N+k) = w[N][k] + (h/2.)*(*(a+i*N+k));
        }
    }
}
