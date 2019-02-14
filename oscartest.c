#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics/graphics.h"
#include <pthread.h>



// Global variables
const double epsilon = 0.001;


inline double computex(double x1, double x2, double constant) {
    return (x1 - x2) * constant;
    
}

inline double computey(double y1, double y2, double constant) {
    return (y1-y2) * constant;
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


    // Create arrays to store values
    double *force = malloc(2*N*sizeof(double));
    double *position = malloc(2*N*sizeof(double));    
    double *velocity = malloc(2* N*sizeof(double));
    double *mass = malloc(N*sizeof(double));
    double *brightness = malloc(N*sizeof(double));

    // Read file
    
    FILE * fp = fopen(filename, "r");

    for(int i = 0; i < N; i++)
    {

        
        fread(&(position[2*i]), sizeof(double), 1, fp);
        fread(&(position[2*i+1]), sizeof(double), 1, fp);
        fread(&(mass[i]), sizeof(double), 1, fp);
        fread(&(velocity[2*i]), sizeof(double), 1, fp);
        fread(&(velocity[2*i+1]), sizeof(double), 1, fp);
        fread(&(brightness[i]), sizeof(double), 1, fp);
        

    }
    fclose(fp);
    int t=0;
    const double G = (double) 100/N;

    if(graphics==1) {
        InitializeGraphics(argv[0], 800, 800);
        SetCAxes(0,1);
        const float W = 1;
        const float L = 1;
        
        while (t <= nsteps) {
            ClearScreen();
            //double *force = malloc(2*N*sizeof(double));
            for (int i = 0; i < N; i++) {
                int p = 2*i;
                DrawCircle(position[2*i], position[2*i+1], W, L, 0.01, 0);
                for(int j = i+1; j < N; j++) {
                    double constant =  -G * (mass[i]) * (mass[j]) / (
                    (sqrt((position[p]-position[2*j])*(position[p]-position[2*j]) + (position[p+1]-position[2*j+1])*(position[p+1]-position[2*j+1])) + epsilon) * 
                    (sqrt((position[p]-position[2*j])*(position[p]-position[2*j]) + (position[p+1]-position[2*j+1])*(position[p+1]-position[2*j+1])) + epsilon) *
                    (sqrt((position[p]-position[2*j])*(position[p]-position[2*j]) + (position[p+1]-position[2*j+1])*(position[p+1]-position[2*j+1])) + epsilon));
                    //printf("Constat = %lf\n", constant); 
                    force[2*i] += computex(position[2*i], position[2*j], constant); 
                    force[2*j] -= force[2*i];
                    force[2*i+1] += computey(position[2*i+1], position[2*j+1], constant);
                    force[2*j+1] -= force[2*i+1];
                }

                //printf("particle: %d x-force =%lf \n",i,  force[2*i]);
                //printf("particle: %d y-force =%lf \n",i , force[2*i+1]);
            
                velocity[2*i] += delta_t*force[2*i]/mass[i];
                position[2*i] += delta_t*velocity[2*i];

                velocity[2*i+1] += delta_t*force[2*i+1]/mass[i];
                position[2*i+1] += delta_t*velocity[2*i+1];

                //printf("particle: %d x-pos =%lf \n",i,  position[2*i]);
                //printf("particle: %d y-pos =%lf \n",i , position[2*i+1]);

                force[2*i] = 0;
                force[2*i+1] = 0;

            }

            //free(force);
            Refresh();
            usleep(3000);
        t++;
        }
        FlushDisplay();
        CloseDisplay();
    } else {

       while (t <= nsteps) {

           //printf("===== %d ===== \n", t);
            //force = (double *)malloc(2*N*sizeof(double));
            //printf("Force after malloc %lf \n", force[3]);
            for (int i = 0; i < N; i++) {
                int p = 2*i;
                for(int j = i+1; j < N; j++) {
                    double constant =  - G * (mass[i]) * (mass[j]) / (
                    (sqrt((position[p]-position[2*j])*(position[p]-position[2*j]) + (position[p+1]-position[2*j+1])*(position[p+1]-position[2*j+1])) + epsilon) * 
                    (sqrt((position[p]-position[2*j])*(position[p]-position[2*j]) + (position[p+1]-position[2*j+1])*(position[p+1]-position[2*j+1])) + epsilon) *
                    (sqrt((position[p]-position[2*j])*(position[p]-position[2*j]) + (position[p+1]-position[2*j+1])*(position[p+1]-position[2*j+1])) + epsilon)); 
                    force[p] += computex(position[p], position[2*j], constant); 
                    force[2*j] -= force[2*i];
                    force[p+1] += computey(position[p+1], position[2*j+1], constant);
                    force[2*j+1] -= force[p+1];

                }

                //printf("particle: %d x-force =%lf \n",i,  force[2*i]);
                //printf("particle: %d y-force =%lf \n",i , force[2*i+1]);

            
                velocity[p] += delta_t*force[p]/mass[i];
                position[p] += delta_t*velocity[p];

                velocity[p+1] += delta_t*force[p]/mass[i];
                position[p+1] += delta_t*velocity[p];

                force[p] = 0;
                force[p] = 0;


            }

            //free(temp);
            //printf("Force after free %lf \n", force[3]);
            
        t++;
        } 
    }
    

    


    // Writing to new file
    FILE * pw = fopen("result.gal", "w");

    for(int i = 0; i < N; i++)
    {
        fwrite(&position[2*i], sizeof(double), 1, pw);
        fwrite(&position[2*i+1], sizeof(double), 1, pw);
        fwrite(&(mass[i]), sizeof(double), 1, pw);
        fwrite(&(velocity[2*i]), sizeof(double), 1, pw);
        fwrite(&(velocity[2*i+1]), sizeof(double), 1, pw);
        fwrite(&(brightness[i]), sizeof(double), 1, pw);

    }
    
    fclose(pw);
return 0;
}