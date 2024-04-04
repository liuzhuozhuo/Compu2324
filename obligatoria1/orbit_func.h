#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Update the acceleration vector, 
//r: matrix (pointer) with the position of all the planets in 3d
//a: matrix (pointer) with the acceleration of all the planets in 3d
//m: vector (pointer) with the modified mass of the planets
//N is the total number of planets
//D: the dimentions of the simulation
void acceleration (double* r, double* a, double* m, int N, int D) {
    int i, j, k;
    double acc; //Dummy variable to sum all the acceleration for a planet in one component
    double R[N][N][D]; 
    double R_mod[N][N]; //Matrix NxN with the module of the distance between planets

    //Calculate the distance vector between planets. 
    for(i = 0; i < N-1; i++){
        for (j = i+1; j < N; j++){
            //Initialize the value to 0, to add recursively the rest of the components
            R_mod[i][j]= 0;
            for (k = 0; k < D; k++){
                R[i][j][k] = *(r+i*D+k) - *(r+j*D+k); 
                R[j][i][k] = -R[i][j][k];
                R_mod[i][j] += pow(R[i][j][k], 2);
            }
            R_mod[j][i] = R_mod[i][j];
        }
    }
    //Calculate the acceleration
    //iterate in the planets
     for(i = 0; i < N; i++){
        //iterate in the components
        for (k = 0; k < D; k++){
            //iterate between planets
            *(a+i*D+k) = 0;
            for (j = 0; j < N; j++){
                if(i!=j){
                    *(a+i*D+k) += -(*(m+j))*R[i][j][k]/pow(R_mod[i][j], 3./2.);
                }
            } 
            
        } 
    }    
}

// Function that calculate the next set of position and velocity using the verlet algorithms
// Using the same notation as before, and the additional parameters are:
// v: matrix (pointer) with the velocity of all the planets in 3d
// h: the step in time for each iteration
// It returns the values by overwriting the previous vectors
void verlet_algorithm(double* r, double* v ,double* a, double* m, int N, int D, double h){
    int i, j, k;
    double w[N][D]; // Angular velocity used for the calculation of the new velocities

    //Calculate the new positions, and the angular velocity.
    for (i = 0; i < N; i++){
        for (k = 0; k < D; k++){
            *(r+i*D+k) = *(r+i*D+k) + (h*(*(v+i*D+k))) + (h * h * (*(a+i*D+k))/2);
            w[i][k] = (*(v+i*D+k)) + (h * (*(a+i*D+k))/2);
        }
    }
   //Using the acceleration function calculate the net acceleration for each of the planets
    acceleration(r, a, m, N, D);

    //Calculate the new velocity of the planets
    for (i = 0; i < N; i++){
        for (k = 0; k < D; k++){
            *(v+i*D+k) = w[i][k] + (h * (*(a+i*D+k))/2);
        }
    }
}

// Function that changes the coordinates so the selected object, l,  is in the origin of coordinates
// Using the same notation as previous functions, 
// r_exp: the positions with the coordinates changed
void change_coord (double* r, double* r_exp, int l, int N, int D){
    //Coordinates of the object l
    double r_l3[3] = {*(r+l*D+0), *(r+l*D+1), *(r+l*D+2)};
    double r_l2[3] = {*(r+l*D+0), *(r+l*D+1)};
    int i, k;

    if (D == 3){
        for(i = 0; i < N; i++){
            for (k = 0; k < D; k++){
                *(r_exp+i*D+k) =  r_l3[k] - *(r+i*D+k);
            }
        }
    }else if (D == 2){
        for(i = 0; i < N; i++){
            for (k = 0; k < D; k++){
                *(r_exp+i*D+k) = r_l2[k] - *(r+i*D+k);
            }
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