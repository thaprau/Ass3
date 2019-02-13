#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct particle{
double posx;
double posy;
double mass;
double vx;
double vy;
double brightness;
}particle;


double distancex(particle *p1, particle *p2){
double d;
d=(p1->posx)-(p2->posx);
//printf("\ndx%lf posx%lf",d,p1->posx);
return d;
}

double distancey(particle *p1, particle *p2){
double d;;
d= (p1->posy)-(p2->posy);
return d;
}

double distance(particle *p1, particle *p2){
double d;
d= sqrt(pow(distancex(p1,p2),2)+pow(distancey(p1,p2),2));
printf("%lf \n", pow(distancey(p1,p2),2));
return d;
}

double forcex(particle** particles, int N, int G, int index){
    double F=0;
    int e0=pow(10,-3);
    for (int i=0; i<N; i++){
        if (i!=index){  
            F += ((particles[i]->mass)/
            pow((distance(particles[index],particles[i])+e0),3))*distancex(particles[index],particles[i]);
        }
    }
    F=F*(-G)*(particles[index]->mass);
    return F;
}




double forcey(particle** particles, int N, int G, int index){
double F=0;
int e0=pow(10,-3);
for (int i=0; i<N; i++){
    if (i!=index){
        F=F+((particles[i]->mass)/pow((distance(particles[index],particles[i])+e0),3))*distancey(particles[index],particles[i]);
    }
}
F=F*(-G)*(particles[index]->mass);
return F;
}


int main(int argc,const char* argv[]){
//check if correct number of input arguments
if (argc!=6){
printf("\nUsage:%s N (int), filename (char), nsteps (int),timestep (double), graphics (int=0,1)",argv[0]);
return -1;
}

//initializations
int N=atoi(argv[1]);
int nsteps=atoi(argv[3]);
int graphics=atoi(argv[5]);
int i,t;

double G=100/N; 
double ax,ay;
double timestep=atof(argv[4]);
double *attributes=malloc(6*sizeof(double));

char const *filename=argv[2]; 
particle **particles=(particle**)malloc(N*sizeof(particle*));

//read file and create N particles
FILE *filein=fopen(filename,"r");

for(i=0; i<N;i++){
fread(attributes,sizeof(double),6,filein);

particle *new = malloc(sizeof(particle));
new->posx=attributes[0];
new->posy=attributes[1];
new->mass=attributes[2];
new->vx=attributes[3];
new->vy=attributes[4];
new->brightness=attributes[5];

particles[i]=new;
}

//printf("\n%lf %lf ",particles[0]->posx,particles[0]->posy);
//printf("\n%lf %lf ",particles[1]->posx,particles[1]->posy);




for(t=0; t<1;t++){
    for(i=0; i<N;i++){
        ax=forcex(particles,N,G,i)/(particles[i]->mass);
        ay=forcey(particles,N,G,i)/(particles[i]->mass);

        //printf("Före ax %lf \n", ax);
        //printf("Före ay %lf \n", ay);
        
        particles[i]->vx = (particles[i]->vx)+ timestep*ax;
        particles[i]->vy = particles[i]->vy+ timestep*ay;

        particles[i]->posx=particles[i]->posx+timestep*particles[i]->vx;
        //printf("\n%lf",particles[i]->posx);
        particles[i]->posy=particles[i]->posy+timestep*particles[i]->vy;



/*particles[i]->posx=particles[i]->posx+xx;
particles[i]->posy=particles[i]->posy+xy;
particles[i]->vx=particles[i]->vx+ux;
particles[i]->vy=particles[i]->vy+uy;*/
}
}

//BARA FÖR ATT PRINTA OCH SE OM DET ÄR RÄTT
//printf("\n%lf %lf %lf %lf %lf %lf",particles[0]->posx,particles[0]->posy,particles[0]->mass,particles[0]->vx,particles[0]->vy,particles[0]->brightness);

//printf("\n%lf %lf %lf %lf %lf %lf",particles[1]->posx,particles[1]->posy,particles[1]->mass,particles[1]->vx,particles[1]->vy,particles[1]->brightness);

char const *filename2="ellipse_N_00010_after200steps.gal";
particle **particles2=(particle**)malloc(N*sizeof(particle*));
FILE *filein2=fopen(filename2,"r");
for(i=0; i<N;i++){
fread(attributes,sizeof(double),6,filein2);
particle *new = malloc(sizeof(particle));
new->posx=attributes[0];
new->posy=attributes[1];
new->mass=attributes[2];
new->vx=attributes[3];
new->vy=attributes[4];
new->brightness=attributes[5];
particles2[i]=new;
}
fclose(filein);

//printf("\n%lf %lf %lf %lf %lf %lf",particles2[0]->posx,particles2[0]->posy,particles2[0]->mass,particles2[0]->vx,particles2[0]->vy,particles2[0]->brightness);

//printf("\n%lf %lf %lf %lf %lf %lf\n",particles2[1]->posx,particles2[1]->posy,particles2[1]->mass,particles2[1]->vx,particles2[1]->vy,particles2[1]->brightness);


return 0;
}


