#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics/graphics.h"
#include <pthread.h>



// Global variables
const double epsilon = 0.001;

// Struct descibing a particle
typedef struct particle {
    double x;
    double y;
    double mass;
    double x_vel;
    double y_vel;
    double brightness;
} part;

// Calculates the force in X-direction of two particles
inline double force_x(part * part1, part * part2, double G) {

    const double r = sqrt(pow((part1->x )- (part2->x), 2 ) + pow((part1->y) - (part2->y), 2));
    double f;
    f = (((part1->x) - (part2->x))) * (G * (part1->mass) * (part2->mass))/pow((r+epsilon),3); 
    return f;
}

// Calculates the force in Y-direction of two particles
inline double force_y(part * part1, part * part2, double G) {
    
    const double r = sqrt(pow((part1->x )- (part2->x), 2 ) + pow((part1->y) - (part2->y), 2));
    double f;
     
    f = (((part1->y) - (part2->y))) * (G * (part1->mass)*(part2->mass))/pow((r+epsilon),3);

    return f;
}

// Calculates the total force in X-direction on a particle
inline double force_tot_x(part ** particles, const int N, const int particleNr, const double G) {

    double f_tot_x = 0;
    
    for(int i = 0; i < N; i++)
    {
        if(particleNr != i) 
        {
        f_tot_x += force_x(particles[particleNr], particles[i], G);
        }
    }
    
    return -f_tot_x;
}

// Calculates the total force in Y-directin on a particle
inline double force_tot_y(part ** particles, int N, int particleNr, double G) {
    
    double f_tot_y = 0;

    for(int i = 0; i < N; i++)
    {
        if(particleNr != i) 
        {
            f_tot_y += force_y(particles[particleNr], particles[i], G);
        }
    }

    return -f_tot_y;
}


inline double acc_x (part ** particles, int N, int particleNr, double G) {
    return  force_tot_x(particles, N, particleNr, G)/(particles[particleNr]->mass);
}

inline double acc_y (part ** particles, int N, int particleNr, double G) {
    return force_tot_y(particles, N, particleNr, G)/(particles[particleNr]->mass);
}


inline void vel(part ** particles, int N, int particleNr, double G, double delta_t) {
    
    particles[particleNr]->x_vel = particles[particleNr]->x_vel + delta_t* acc_x(particles,N,particleNr,G); 
    particles[particleNr]->y_vel = particles[particleNr]->y_vel + delta_t* acc_y(particles,N,particleNr,G);
}

inline void pos(part ** particles, int particleNr, double delta_t) {
    particles[particleNr]->x += delta_t* (particles[particleNr]->x_vel); 
    particles[particleNr]->y += delta_t* (particles[particleNr]->y_vel);
}





int main(int argc, char* argv[]) {

    if (argc != 6) {
        printf("Wrong number of inputs");
        return -1;
    }
    // Read user input
    const int N = atoi(argv[1]);
    const char* filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const int graphics = atoi(argv[5]);

    // Creates array to store particles
    part **array = (part**) malloc(N*sizeof(part*));
    for(int i = 0; i < N; i++)
    {
        array[i] = (part*) malloc(sizeof(part));
    }
    

    // Read file
    
    FILE * fp = fopen(filename, "r");

    for(int i = 0; i < N; i++)
    {

        fread(&(array[i]->x), sizeof(double), 1, fp);
        fread(&(array[i]->y), sizeof(double), 1, fp);
        fread(&(array[i]->mass), sizeof(double), 1, fp);
        fread(&(array[i]->x_vel), sizeof(double), 1, fp);
        fread(&(array[i]->y_vel), sizeof(double), 1, fp);
        fread(&(array[i]->brightness), sizeof(double), 1, fp);

    }


    fclose(fp);


    // Prepare for the loop
    double t = 0;
    const double G = (double) 100/N;
    const float W = 1; 
    const float L = 1;

    // Loop with graphics
    if(graphics==1) {
        InitializeGraphics(argv[0], 800, 800);
        SetCAxes(0,1);

        while(t <= nsteps){
            ClearScreen();
            
            for(int i = 0; i < N; i++)
            {
                DrawCircle(array[i]->x, array[i]->y, W, L, 0.01, 0);
                vel(array, N, i, G, delta_t);
                pos(array, i, delta_t);

            }
            

            Refresh();
            usleep(3000);

            t +=1;
        }
        FlushDisplay();
        CloseDisplay();
    }

    // Loop without graphics
    else {
        while(t <= nsteps) {
            for(int i = 0; i < N; i++)  
            {
                vel(array, N, i, G, delta_t);
                pos(array, i, delta_t);
            }
        t += 1;                                
        }

    }

    // Writing to new file
    FILE * pw = fopen("result.gal", "w");

    for(int i = 0; i < N; i++)
    {
        fwrite(&(array[i]->x), sizeof(double), 1, pw);
        fwrite(&(array[i]->y), sizeof(double), 1, pw);
        fwrite(&(array[i]->mass), sizeof(double), 1, pw);
        fwrite(&(array[i]->x_vel), sizeof(double), 1, pw);
        fwrite(&(array[i]->y_vel), sizeof(double), 1, pw);
        fwrite(&(array[i]->brightness), sizeof(double), 1, pw);

    }
    
    fclose(pw);

    for(int i = 0; i < N; i++)
    {
        free(array[i]);
    }
    free(array);
return 0;
}