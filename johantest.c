
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics/graphics.h"
#include <pthread.h>



// Global variables
const double epsilon = 0.001;


// Struct descibing a particle
typedef struct particle {
    double x_force;
    double y_force;
    double x;
    double y;
    double mass;
    double x_vel;
    double y_vel;
    double brightness;
} part;

inline void force_x(part *p1, part *p2, double r, double G) {
 

    double force = G * (p1->mass) * (p2->mass) * ((p1->x) - (p2->x)) /
    pow(r + epsilon, 3);

    p1->x_force = p1->x_force - force;
    p2->x_force = p2->x_force + force;

}

inline void force_y(part *p1, part *p2, double r, double G) {

    double force = G * (p1->mass) * (p2->mass) * ((p1->y) - (p2->y)) /
    pow(r + epsilon, 3);

    p1->y_force = p1->y_force - force;
    p2->y_force = p2->y_force + force;
}

void vel_update(part** particles, int N, double t) {

    for(int i = 0; i < N; i++) {
        particles[i]->x_vel = particles[i]->x_vel + ((particles[i]->x_force)/(particles[i]->mass)) * t;
        particles[i]->y_vel = particles[i]->y_vel + ((particles[i]->y_force)/(particles[i]->mass)) * t;    
    }
}

void pos_update(part** particles, int N, double t) {

        for(int i = 0; i < N; i++) {
            particles[i]->x = particles[i]->x + (particles[i]->x_vel) * t;
            particles[i]->y = particles[i]->y + (particles[i]->y_vel) * t;
            //printf("X-force = %lf \n", particles[i]->x);
            //printf("Y-force = %lf \n", particles[i]->y);
            particles[i]->y_force = 0;
            particles[i]->x_force = 0;

    }

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
    int t = 0;
    const double G = (double) 100/N;
    const float W = 1; 
    const float L = 1;

    // Loop with graphics
    if(graphics==1) {
        InitializeGraphics(argv[0], 800, 800);
        SetCAxes(0,1);

        while(t <= nsteps){
            ClearScreen();
            printf("Loop nr %d \n", t);

            for(int i = 0; i < N; i++) 
            {
                array[i]->y_force = 0;
                array[i]->x_force = 0;
            }

            for(int i = 0; i < N; i++)  
            {
                
                DrawCircle(array[i]->x, array[i]->y, W, L, 0.01, 0);

                for(int j = i+1; j < N; j++)
                {
                    double r = sqrt(pow((array[i]->x)-(array[j]->x), 2) + pow((array[i]->y)-(array[j]->y), 2));
                    force_x(array[i], array[j], r, G );
                    force_y(array[i], array[j], r, G);

                }
            }


            vel_update(array, N, delta_t);
            pos_update(array, N, delta_t);
            
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
                for(int j = i+1; j < N; j++)
                {

                    double r = sqrt(pow((array[i]->x)-(array[j]->x), 2) + pow((array[i]->y)-(array[j]->y), 2));
                    force_x(array[i], array[j], r, G );
                    //force_y(array[i], array[j], r, G);
                }
            }
            vel_update(array, N, delta_t);
            pos_update(array, N, delta_t);
            t +=1;
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