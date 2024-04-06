#include <stdio.h>
#include "orbit_func.h"

int main(){
    // Modified arguments for different simulations
    double total_time =  10;// in years
    double interval = 1; // in days

    int N = 10; //Number of planets
    int D = 2; // Dimensions of the simulation

    // Open folder where the initial conditions are stored, and the folder to write the data.
    // The initial condition file has format, where each line --> mass   x   y   z   v_x   v_y   v_z, 
    // separated by tab of the different planets

    FILE *f_init, *f_exp, *f_geo; //pointer to the files

    f_init = fopen("data/init_cond.txt", "r"); // open the init_cond file to read
    f_exp = fopen("data/planets_data.txt", "w"); // open the export file to write
    f_geo = fopen("data/geocentric_data.txt", "w"); // open the export file to write

    // Define t of the simulation, and the step h
    double t = total_time * 3.154e7; //s
    double h= interval * 8.64e4; //s

    // Rescale the magnitudes and define some constants
    double G = 6.67408e-11; //N km2/kg2, Gravitational constant
    double c = 1.4966e11; //km, astronomical units
    double M_s = 1.99e30; //kg, solar mass
    double t_prime = pow(G*M_s/pow(c, 3.),0.5) ;// Rescalation constant

    int i, j, k; // Used for the loop

    // Define the matrix for position, velociy and acceleration, and it's pointers
    double r[N][D];
    double r_geo[N][D];
    double v[N][D];
    double a[N][D];
    double m[N];

    double* r_p = r[0];
    double* r_pg = r_geo[0];
    double* v_p = v[0];
    double* a_p = a[0];
    double* m_p = m;

    //Read the file and store the data in r and v
    i = 0;
    while (fscanf(f_init, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
            &m[i] , &r[i][0], &r[i][1], &r[i][2], &v[i][0], &v[i][1], &v[i][2])!= EOF){
        i++;
    }
    
    // From those data rescale them
    h = t_prime*h;
    t = t*t_prime;
    for (i=0; i<N; i++){
        m[i] = m[i]/M_s;
        for (k=0; k<D; k++){
            r[i][k] = r[i][k]/c;
            v[i][k] = v[i][k]/(c*t_prime);
        }
    }

    //Calculate the first set of acceleration used for the Verlet algorithm
    acceleration(r_p, a_p, m_p, N, D);

    for (i = 0; i < (int)(t/h); i++){
        // Using the Verlet algorithms get the new position for the planets
        verlet_algorithm(r_p, v_p, a_p, m_p, N, D, h);
        // Using the function for the change of coordinates, get the position of the planets in terms of the Earth
        change_coord(r_p, r_pg, 3, N, D);
        if (D == 2){
            for (j = 0; j < N; j++){
            // Print the positions in the file
            fprintf(f_exp, "%lf, %lf\n", r[j][0], r[j][1]); // For 2D
            fprintf(f_geo, "%lf, %lf\n", r_geo[j][0], r_geo[j][1]); // For 2D
            }
        }else if (D == 3){
            // Print the positions in the file
            fprintf(f_exp, "%lf, %lf, %lf \n", r[j][0], r[j][1], r[j][2]); // For 3D
            fprintf(f_geo, "%lf, %lf, %lf \n", r_geo[j][0], r_geo[j][1], r_geo[j][2]); // For 2D
        }

        fprintf(f_exp, "\n");
        fprintf(f_geo, "\n");
    }

    fclose(f_init);
    fclose(f_exp);
    fclose(f_geo);

    return(0);
}