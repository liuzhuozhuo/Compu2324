#include <stdio.h>
#include "orbit_func.h"
#include <time.h>

int main(){
    clock_t begin = clock();

    // Modified arguments for different simulations
    double total_time =  250;// in years
    double interval = 3; // in days

    int N = 10; //Number of planets
    int D = 2; // Dimensions of the simulation
    int q = 5; // Every q frames the data gets logged.

    // Open folder where the initial conditions are stored, and the folder to write the data.
    // The initial condition file has format, where each line --> mass   x   y   z   v_x   v_y   v_z, 
    // separated by tab of the different planets

    FILE *f_init, *f_exp, *f_geo, *f_energy_t, *f_energy, *f_period; //pointer to the files

    f_init = fopen("data_time_joel/init_cond.txt", "r"); // open the init_cond file to read
    f_exp = fopen("data_time_joel/planets_data_.txt", "w"); // open the export file to write the planet orbit
    f_geo = fopen("data_time_joel/geocentric_data_.txt", "w"); // open the export file to write the planet orbit for geocentric
    f_energy_t = fopen("data_time_joel/energy_total.txt", "w"); // open the export file to write the total energy of the system
    f_energy = fopen("data_time_joel/energy.txt", "w"); // open the export file to write the total energy of each planet
    f_period = fopen("data_time_joel/period.txt", "w"); // open the export file to write the period for each planet

    // Define t of the simulation, and the step h
    double t = total_time * 3.154e7; //s
    double h= interval * 8.64e4; //s

    // Rescale the magnitudes and define some constants
    double G = 6.67408e-11; //N m2/kg2, Gravitational constant
    double c = 1.4966e11; //m, astronomical units
    double M_s = 1.99e30; //kg, solar mass
    double t_prime = pow(G*M_s/pow(c, 3.),0.5) ;// Rescalation constant

    int i, j, k; // Used for the loop

    // Define the matrix for position, velocity, acceleration, and other variables, and it's pointers
    double r[21][D];
    double r_geo[N][D];
    double v[21][D];
    double a[N][D];
    double R_mod[N][N];
    double m[N];
    double T_total, V_total;
    double T[N], V[N];

    double curr_angle[N], prev_angle[N], init_angle[N];
    double n_turns[N];
    double period[N];

    double* r_p = r[0];
    double* r_pg = r_geo[0];
    double* v_p = v[0];
    double* a_p = a[0];
    double* R_mod_p = R_mod[0];
    double* m_p = m;
    double* T_p = T;
    double* V_p = V;

    double* c_angle_p = curr_angle;
    double* p_angle_p = prev_angle;

    //Read the file and store the data in r and v
    i = 0;
    while (fscanf(f_init, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
            &m[i] , &r[i][0], &r[i][1], &r[i][2], &v[i][0], &v[i][1], &v[i][2])!= EOF){
        i++;
    }
    
    // From those data rescale them
    h = t_prime*h;
    t = t*t_prime;
    rescale(r_p, v_p, a_p, m_p, N, D, t_prime, M_s, c);

    //Calculate the current angle of the different planets:
    angle(r_p, c_angle_p, N);

    // Asign the calculated angle to the array of initial angles
    for (i = 1; i < N; i++){
        n_turns[i] = 0;
        init_angle[i] = curr_angle[i];
    }
    
    //Calculate the first set of acceleration used for the Verlet algorithm
    acceleration(r_p, a_p, m_p, R_mod_p, N, D);

    //Iterate in time 
    for (i = 0; i < (int)(t/h); i++){
        // Using the Verlet algorithms get the new position for the planets
        verlet_algorithm(r_p, v_p, a_p, m_p, R_mod_p, N, D, h);

        // Using the function for the change of coordinates, get the position of the planets in terms of the Earth
        change_coord(r_p, r_pg, 3, N, D);
        
        // Depending on the dimention used, the data is written accordingly
        if (D == 2 && i%q == 0){
            for (j = 0; j < N; j++){
                // Print the positions in the file
                fprintf(f_exp, "%lf, %lf\n", r[j][0], r[j][1]); // For 2D
                fprintf(f_geo, "%lf, %lf\n", r_geo[j][0], r_geo[j][1]); // For 2D
            }
            fprintf(f_exp, "\n");
            fprintf(f_geo, "\n");

        }else if (D == 3 && i%q == 0){
            for (j = 0; j < N; j++){
                // Print the positions in the file
                fprintf(f_exp, "%lf, %lf, %lf \n", r[j][0], r[j][1], r[j][2]); // For 3D
                fprintf(f_geo, "%lf, %lf, %lf \n", r_geo[j][0], r_geo[j][1], r_geo[j][2]); // For 2D
            }
            fprintf(f_exp, "\n");
            fprintf(f_geo, "\n");
        }

        // Asign the values for the angles of the previous iteration
        for (j = 1; j < N; j++){
            prev_angle[j] = curr_angle[j];
        }
        // Calculate the angles in the new positions
        angle(r_p, c_angle_p, N);

        for (j = 1; j < N; j++){
            // Check if a planet has completed a full turn, if so, calculate the period the orbit in years
            if (turn_count(prev_angle[j], curr_angle[j], init_angle[j])){
                n_turns[j] ++;
                period[j] = i*h/((n_turns[j])*t_prime*3.154e7);
            }
        }

        // Calculate the total kinetic and potential energy of the system
        T_total = k_energy(v_p, m_p, N, D);
        V_total = p_energy(R_mod_p, m_p, N, D);

        // Calculate the kinetic and potential energy of each planet
        energy(v_p, m_p, R_mod_p, T_p ,V_p, N, D);
        for (j = 0; j < N; j++){
            T[j] = derescaled_energy(T[j], t_prime, M_s, c);
            V[j] = derescaled_energy(V[j], t_prime, M_s, c);
            fprintf(f_energy, "%lf, %lf, %lf, %lf\n", i*h, T[j], V[j], T[j]+V[j]);
        }

        if(i != (int)(t/h)-1){
            fprintf(f_energy, "\n");
        }
        

        fprintf(f_energy_t, "%lf, %lf, %lf, %lf\n", i*h, T_total, V_total, T_total+V_total);
        
    }

    for(i = 0; i < N; i++){
        if(n_turns[i] != 0){
            fprintf(f_period, "%lf, %lf\n", n_turns[i], period[i]);
        }else{
            fprintf(f_period, "%lf, %lf\n", n_turns[i], 0.);
        }
    }

    fclose(f_init);
    fclose(f_exp);
    fclose(f_geo);
    fclose(f_energy);
    fclose(f_energy_t);
    fclose(f_period);


    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("%lf\n", time_spent);

    return(0);
}