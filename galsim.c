#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "graphics.h"
 

const float circleRadius=0.0025, circleColor=0;
const int windowWidth=800;
const float L=1, W=1; //opt

typedef struct particle{	
	double posx;
	double posy;
    	double newposx;
	double newposy;
	double velx;
	double vely;
	double forcex;
	double forcey;
	double mass;
	double bri;
}particle;

double *distance(double x1, double y1, double x2, double y2){
	double *r = (double*)malloc(2*sizeof(double));
	if (!r){
		return NULL;
	}
	r[0] = (x1-x2);
	r[1] = (y1-y2);
	return r;
}

void force (double G, int N, particle** allPart, int ind){
	int i;
	double eps = 0.001;
	particle *p1 = allPart[ind];
	double m = p1->mass;
	for (i=ind+1; i<N; i++){
		particle *p2 = allPart[i];
		double *r_b = distance (p1->posx, p1->posy, p2->posx, p2->posy);
		double  r = sqrt((r_b[0]*r_b[0]+ r_b[1]*r_b[1])); //opt
		double temp = G*m*(p2->mass)/((r+eps)*(r+eps)*(r+eps));  //opt
		p1->forcex -= temp*r_b[0];
		p1->forcey -= temp*r_b[1];
		p2->forcex += temp*r_b[0];
		p2->forcey += temp*r_b[1];
		free(r_b);
	}
}

void add(double a[6], particle **p){
	(*p)->posx = a[0];
	(*p)->posy = a[1];
	(*p)->mass = a[2];
	(*p)->velx = a[3];
	(*p)->vely = a[4];
	(*p)->bri = a[5];
	(*p)->forcex = 0;
	(*p)->forcey = 0;
}

void writeOut(double* a, particle *p){
	a[0] = p->posx;
	a[1] = p->posy;
	a[2] = p->mass;
	a[3] = p->velx;
	a[4] = p->vely;
	a[5] = p->bri;
}

int main(int argc, char *argv[])
{
	if (argc != 6){
		printf("\nUsage %s N filename nsteps delta_t graphics\n", argv[0]);
		return -1;
	}

	//Store the input
	int N = atoi(argv[1]);
	char const *filename = argv[2];
	int nsteps = atoi(argv[3]);  //opt: put const on all
	double delta_t = atof(argv[4]);
	int graph = atoi(argv[5]);
	double G = (double)100/N;

	//create pointer to particle pointers
	particle **allPart = (particle **)malloc(N*sizeof(particle*));
	if(!allPart){
		printf("Couldn't allocate allPart\n");
		return -1;
	}
	int i, k;

	FILE *input = fopen(filename, "r");
	for (i =0; i<N; i++){
		double attr[6] = {0}; //compiler will set all values to zero
		fread(attr,sizeof(double),6,input);
		allPart[i] = (particle*)malloc(sizeof(particle));
		if (!allPart[i]){
			printf("Couldn't allocate allPart[%d]", i);
			return -1;
		}
		add(attr, &allPart[i]);
	}
		fclose(input);

    if(graph){
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    
    	for (k=0; k<nsteps; k++){
    		//printf("k: %d", k);
            ClearScreen();
	    	for(i=0; i<N; i++){
	    		particle *p = allPart[i]; //opt restricted?
	    		double m = (p->mass);
			force(G, N, allPart, i);
	    		p->velx = p->velx + delta_t*(p->forcex)/m;
	    		p->vely = p->vely + delta_t*(p->forcey)/m;
	    		p->newposx = p->posx + delta_t*(p->velx);
	    		p->newposy = p->posy + delta_t*(p->vely);
                DrawCircle(p->posx,p->posy , L, W, circleRadius, circleColor);
            }
            for (i=0; i<N; i++){
                particle *p = allPart[i];
                p->posx = p->newposx;
                p->posy = p->newposy;
		p->forcex = 0;
		p->forcey = 0;
            }
            Refresh();
            usleep(30000);
        }
        FlushDisplay();
        CloseDisplay();
    
    } else {
        for (k=0; k<nsteps; k++){
    		//printf("k: %d \n", k);
    		for(i=0; i<N; i++){
    			particle *p = allPart[i]; //opt restricted?
    			double m = (p->mass);
    			force(G, N, allPart, i);
	    		p->velx = p->velx + delta_t*(p->forcex)/m;
	    		p->vely = p->vely + delta_t*(p->forcey)/m;
    			p->newposx = p->posx + delta_t*(p->velx);
    			p->newposy = p->posy + delta_t*(p->vely);
    		}
              for (i=0; i<N; i++){
                particle *p = allPart[i];
                p->posx = p->newposx;
                p->posy = p->newposy;
		p->forcex = 0;
		p->forcey = 0;
            }
    	}
    }



	FILE *output = fopen("result.gal","wb");
	assert(output);
	for(i=0; i<N; i++){
		double attr[6] = {0};
		writeOut(attr, allPart[i]);
		fwrite(attr, sizeof(double), 6, output);
	}

	for (i=0; i<N; i++){
		free(allPart[i]);
	} free(allPart);

    return 0;
}
