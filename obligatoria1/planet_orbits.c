#include <stdio.h>
#include "orbit_func.h"

int main(){
        // Define t of the simulation, and the step h
    double t = 1.;
    double h= 0.01;

    // Rescale the magnitudes and define some constants
    double G = 6.67408e-11; //N m2/kg2, Gravitational constant
    double c = 1.4966e8; //km, astronomical units
    double M_s = 1.99e30; //kg, solar mass
    double t_prime = pow(G*M_s/pow(c, 3.),0.5) ;//
    int N = 5; //Number of planets

    int i, j, k; // Used for the loop

    // Define the matrix for position, velociy and acceleration, and it's pointers
    long double r[N][3];
    long double v[N][3];
    long double a[N][3];
    double m[N];

    long double* r_p = r[0];
    long double* v_p = v[0];
    long double* a_p = a[0];
    double* m_p = m;

    // Open folder where the initial conditions are stored, and the folder to write the data.
    // The initial condition file has format, where each line --> mass   x   y   z   v_x   v_y   v_z, separated by tab of a different planet
    // Data taken from 21/03/2024 at 00:00 (https://ssd.jpl.nasa.gov/horizons/app.html#/)

    FILE *f_init, *f_exp; //pointer to the files

    f_init = fopen("initial_condition.csv", "r"); // open the init_cond file to read
    f_exp = fopen("planets_data.txt", "w"); // open the export file to write

    //Read the file and store the data in r and v
    i = 0;
/*     while (fscanf(f_init, "%lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf",&m[i] , &r[i][0], &r[i][1], &r[i][2], &v[i][0], &v[i][1], &v[i][2])!= EOF){
        i++;
    } */
    fscanf(f_init, "");

    
    // From those data rescale them
    for (i=0; i<N; i++){
        for (k=0; k<3; k++){
            r[i][k] = r[i][k]/c;
            v[i][k] = v[i][k]*t_prime/c;
            printf("%Lf, ", r[i][k]);
            printf("%Lf, ", v[i][k]);
        }
        m[i] = m[i]/M_s;
        printf("%lf\n", m[i]);
    }

    //Calculate the first set of acceleration used for the Verlet algorithm
    acceleration(r_p, a_p, m_p, N);

    for (i = 0; i < (int)(t/h); i++){
        verlet_algorithm(r_p, v_p, a_p, m_p, N, h);
        for (j = 0; j < N; j++){
            //fprintf(f_exp, "%Lf, %Lf, %Lf \n", r[j][0], r[j][1], r[j][2]);
            fprintf(f_exp, "%Lf, %Lf\n", r[j][0], r[j][1]);
        }
        fprintf(f_exp, "\n");
    }
    fclose(f_init);
    fclose(f_exp);

    return(0);
}