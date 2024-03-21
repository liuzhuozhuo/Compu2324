#include <stdio.h>
#include "orbit_func.h"

int main(){
        // Define t of the simulation, and the step h
    double t = 1.;
    double h= 0.1;

    // Rescale the magnitudes and define some constants
    double G = 6.67408e-11; //N m2/kg2, Gravitational constant
    double c = 1.4966e11; //m, astronomical units
    double M_s = 1.99e30; //kg, solar mass
    int N = 3; //Number of planets

    int i, j, k; // Used for the loop

    // Define the matrix for position, velociy and acceleration, and it's pointers
    double r[N][3];
    double v[N][3];
    double a[N][3];
    double m[3] = {1, 1, 1};

    double* r_p = r[0];
    double* v_p = v[0];
    double* a_p = a[0];
    double* m_p = m;

    // Open folder where the initial conditions are stored, and the folder to write the data.
    FILE *f_init, *f_exp; //pointer to the files

    f_init = fopen("initial_cond.txt", "r"); // open the init_cond file to read
    f_exp = fopen("data_orbit.txt", "w"); // open the export file to write

    //Read the file and store the data in r and v
    i = 0;
    while (fscanf(f_init, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &r[i][0], &r[i][1], &r[i][2], &v[i][0], &v[i][1], &v[i][2])!= EOF){
        i++;
    }
    
    // From those data rescale them


    //Calculate the first set of acceleration used for the Verlet algorithm
    acceleration(r_p, a_p, m_p, N);
    for (i=0; i<N; i++){
        printf("%lf, ", a[i][0]);
    }

    for (i = 0; i < (int)(t/h); i++){
        verlet_algorithm(r_p, v_p, a_p, m_p, N, h);
        printf("Iteracion %i: \n", i);
        for (j = 0; j < N; j++){
            printf("Planeta %i: %lf, %lf, %lf \n",j , r[j][0], r[j][1], r[j][2]);
        }

    }
    fclose(f_init);
    fclose(f_exp);

    return(0);
}