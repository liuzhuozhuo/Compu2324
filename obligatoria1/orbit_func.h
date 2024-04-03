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
            //Initialize the value to 0, to add recursively the rest of the components
            R_mod[i][j]= 0;
            for (k = 0; k < 3; k++){
                R[i][j][k] = *(r+i*3+k) - *(r+j*3+k); 
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
                // Considering that only the upper diagonal of the matriz is calculated, and to avoid the case i=i, 
                // the following conditions are checked.
                if(j>i){
                    acc += -(*(m+j*3))*R[i][j][k]/pow(R_mod[i][j], 3./2.);
                }
                if(j<i){
                    acc += (*(m+j*3))*R[j][i][k]/pow(R_mod[j][i], 3./2.);
                }
            } 
            //Assing the acceleration checked to it's position in the acceleration matrix.
            *(a+i*3+k) = acc;
        } 
    }    
}

// Function that calculate the next set of position and velocity using the verlet algorithms
// Using the same notation as before, and the additional parameters are:
// v: matrix (pointer) with the velocity of all the planets in 3d
// h: the step in time for each iteration
// It returns the values by overwriting the previous vectors
void verlet_algorithm(double* r, double* v ,double* a, double* m, int N, double h){
    int i, j, k;
    double w[N][3]; // Angular velocity used for the calculation of the new velocities

    //Calculate the new positions, and the angular velocity.
    for (i = 0; i < N; i++){
        for (k = 0; k < 3; k++){
            *(r+i*3+k) = *(r+i*3+k) + (h*(*(v+i*3+k))) + (h * h * (*(a+i*3+k))/2);
            w[i][k] = (*(v+i*3+k)) + (h * (*(a+i*3+k))/2);
        }
    }
   //Using the acceleration function calculate the net acceleration for each of the planets
    acceleration(r, a, m, N);

    //Calculate the new velocity of the planets
    for (i = 0; i < N; i++){
        for (k = 0; k < 3; k++){
            *(v+i*3+k) = w[i][k] + (h * (*(a+i*3+k))/2);
        }
    }
}

// Function that changes the coordinates so the selected object, l,  is in the origin of coordinates
// Using the same notation as previous functions, 
void change_coord (double* r, int l,int N){
    //Coordinates of the object l
    double r_l[3] = {*(r+l*N+0), *(r+l*N+1), *(r+l*N+2)};
    int i, k;

    for(i = 0; i < N; i++){
        for (k = 0; k < N; k++){
            *(r+i*N+k) -= r_l[k];
        }
    }
}

//Function that calculate the total energy of the system to check if the energy is conserved.
// double energy (double* r, double*v, double* m, int N){
//     int i, j, k;
//     double T = 0, V= 0;
//     double v_mod_sq[N];



//     for (i = 0; i < N; i++){
//         for (k = 0; k < 3; k++){

//         }
//         T += (1./2.)*(*(m+i))*v_mod_sq[i];
//     }
    
    
// }

// Function that calculate the angle, from the position of a planet.
// In this case r_i indicate a vector of size 3, NOT A MATRIX
double angle(double* r_i){
    return atan(*(r_i+1)/(*(r_i)));
}

//Function that returns the number of turns that the planet has done.
double turn_count(){
    return 0;
}