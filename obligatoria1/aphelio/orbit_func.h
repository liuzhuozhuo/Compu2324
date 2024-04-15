#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

//Update the acceleration vector, 
//r: matrix (pointer) with the position of all the planets in 3d
//a: matrix (pointer) with the acceleration of all the planets in 3d
//m: vector (pointer) with the modified mass of the planets
//R_mod: metrix (pointer) with the distance between planets, calculated by the function
//N is the total number of planets
//D: the dimentions of the simulation
void acceleration (double* r, double* a, double* m, double* R_mod, int N, int D) {
    int i, j, k;
    double acc; //Dummy variable to sum all the acceleration for a planet in one component
    double R[N][N][D]; 

    //Calculate the distance vector between planets. 
    for(i = 0; i < N-1; i++){
        for (j = i+1; j < N; j++){
            //Initialize the value to 0, to add recursively the rest of the components
            *(R_mod+i*N+j)= 0;
            for (k = 0; k < D; k++){
                R[i][j][k] = *(r+i*D+k) - *(r+j*D+k); 
                R[j][i][k] = -R[i][j][k];
               *(R_mod+i*N+j) += pow(R[i][j][k], 2);
            }
            *(R_mod+j*N+i)= *(R_mod+i*N +j);
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
                    *(a+i*D+k) += -(*(m+j))*R[i][j][k]/pow(*(R_mod+i*N+j), 3./2.);
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
void verlet_algorithm(double* r, double* v ,double* a, double* m, double* R_mod, int N, int D, double h){
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
    acceleration(r, a, m, R_mod, N, D);

    //Calculate the new velocity of the planets
    for (i = 0; i < N; i++){
        for (k = 0; k < D; k++){
            *(v+i*D+k) = w[i][k] + (h * (*(a+i*D+k))/2);
        }
    }
}

// Function that rescales the values of distance, velocity and mass
void rescale (double* r, double* v ,double* a, double* m, 
                    int N, int D, double t_prime, double M_s, double c){

    int i, k;
    for (i=0; i<N; i++){
        m[i] = m[i]/M_s;
        for (k=0; k<D; k++){
            *(r+i*D+k) = *(r+i*D+k)/c;
            *(v+i*D+k) = *(v+i*D+k)/(c*t_prime);
        }
    }
}

// Function that derescales the values of distance, velocity and mass
void derescale(double* r, double* v ,double* a, double* m, double t, 
                    int N, int D, double h, double t_prime, double M_s, double c){

    int i, k;
    for (i=0; i<N; i++){
        m[i] = m[i]*M_s;
        for (k=0; k<D; k++){
            *(r+i*D+k) = *(r+i*D+k)*c;
            *(v+i*D+k) = *(v+i*D+k)*(c*t_prime);
        }
    }
}

// Function that returns the value of the energy in units of international system
double derescaled_energy (double E, double t_prime, double M_s, double c){
    return E * M_s *pow(c, 2) /pow(t_prime, 2);
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
                *(r_exp+i*D+k) =  -r_l3[k] + *(r+i*D+k);
            }
        }
    }else if (D == 2){
        for(i = 0; i < N; i++){
            for (k = 0; k < D; k++){
                *(r_exp+i*D+k) = -r_l2[k] + *(r+i*D+k);
            }
        }
    }
    
}

//Function that calculate the energy for each of the planets of the system to check if the energy is conserved.
// T: vector (pointer) with the kinetic energy of the planets
// V: vector (pointer) with the potential energy of the planets
void energy (double*v, double* m, double* R_mod, double* T, double* V,int N, int D){
    int i, j, k;
    for (i = 0; i < N; i++){
        *(T+i) = 0;
        *(V+i) = 0;
        for (k = 0; k < D; k++){
            *(T+i) += *(m+i)*(*(v+i*D+k))*(*(v+i*D+k))/2.;
        }
        for (j = 0; j < N; j++){
            if(i != j) {
                *(V+i) -= (*(m+i))*(*(m+j))/pow(*(R_mod+i*N+j), 0.5);
            }
        }
    }
}

//Returns the total kinetic energy of the system
double k_energy (double*v, double* m, int N, int D){
    int i, j, k;
    double T = 0;

    for (i = 0; i < N; i++){
        for (k = 0; k < D; k++){
            T += *(m+i)*pow(*(v+i*D+k), 2)/2.;
        }        
    }
    return T;
}

//Returns the total potential energy of the system
double p_energy (double* R_mod, double* m, int N, int D){
    int i, j, k;
    double V = 0;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if(i != j) {
                V -= (*(m+i))*(*(m+j))/pow(*(R_mod+i*N+j), 0.5);
            }
        } 
    }
    return V;
}

// Function that calculate the angle, from the position of a planet, NOT THE SUN. (Only valid on 2D)
void angle(double* r, double* curr_angle, int N){
    int i;
    for (i = 1; i < N; i++){
        if (*(r+i*2) < 0){
            *(curr_angle+i) = -atan(*(r+i*2+1)/(*(r+i*2)));
        }else{
            *(curr_angle+i) = atan(*(r+i*2+1)/(*(r+i*2)));
        }
    }   
}

//Function that checks if the planet has completed an turn.
bool turn_count(double prev_angle, double curr_angle, double init_angle){
    if (curr_angle > init_angle && prev_angle < init_angle){
        return true;
    }else{
        return false;
    }
}