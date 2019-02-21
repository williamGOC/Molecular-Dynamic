#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

#define L 54.77225575
#define NC 30
#define RCUT 2.5
#define RHO 0.3
#define EMPTY -1
#define N 900
#define square(A) (A)*(A)


typedef double vector_R[2];
typedef vector_R *vector_RN;
typedef void(* force)(vector_RN, vector_RN);

typedef struct particles{
    force f;
    vector_RN position, momentum, aceleration;
    double pass_time;
    int accountant;
} class_system;


class_system *builder_system(double pass_time, force f);
double energy_K(class_system *sys);
double energy_P(class_system *sys);
double histogram_system(vector_RN v);
double dist_ij(vector_R ri, vector_R rj);
double minimal_image(double x_1, double x_2);
void f(vector_RN pos, vector_RN fu);
void calculate_fij(vector_R ri, vector_R rj, vector_R fij);
void verlet_iteration(class_system *sys);
void periodic_condition(class_system *sys);
void printer_system(class_system *sys, FILE *gnuplot_pipe);
void destroy_system(class_system *sys);
void inversor();


void main(){
    int i;
    double K, U;
    class_system *sys;
    
    FILE *feop = fopen("Energy","w");
    FILE *pipe = popen("gnuplot", "w");
    assert(pipe);

    sys = builder_system(0.005, f);
 
    for (i=0;i<2000;i++){ 
        verlet_iteration(sys);
        periodic_condition(sys);

        K = energy_K(sys);
        U = energy_P(sys);
        fprintf(feop,"%g\t%lf\t%lf\t%lf\n", i * 0.005, K, U, K + U);
        histogram_system(sys->momentum);
        printer_system(sys, pipe);
    }
  
    fclose(pipe);
    fclose(feop);
    destroy_system(sys);

    //inversor();
}


/*
Construye un objeto de tipo class_system * y lo inicializa adecuadamente.
Keyword Arguments:
pass_time: double (paso de tiempo del iterador).
f: force (fuerza de interacción del sistema).
*/
class_system *builder_system(double pass_time, force f){
    int i;

    class_system *sys = (class_system *)malloc(sizeof(class_system));        
    assert(sys);

    sys -> position = (vector_RN)malloc(N * sizeof(vector_R));
    assert(sys -> position);
    sys -> momentum = (vector_RN)malloc(N * sizeof(vector_R));
    assert(sys -> momentum);
    sys -> aceleration = (vector_RN)malloc(N * sizeof(vector_R));
    assert(sys -> aceleration);
  
    sys -> f = f;
    sys -> pass_time = pass_time;
    sys -> accountant = 0;
  
    for (i = 0; i < N; i++){
        sys->position[i][0] = (i % NC) * L / (NC + 1);
        sys->position[i][1] = (i / NC) * L / (NC + 1);
        sys->momentum[i][0] = ((rand() % 2) ? 1.1 : -1.1);
        sys->momentum[i][1] = 0;
    }

    f(sys -> position,sys -> aceleration);
    return sys;
}


/*
Fuerza de interacción de las partículas del sistema.
Keyword Arguments:
pos: vector_RN (posiciones instantaneas de las partículas del sistema).
fu: vector_RN (almacenamiento de las fuerzas de interacción).
*/
void f(vector_RN pos, vector_RN fu){
    int i,j,k;
    vector_R fij;
  
    for(i = 0; i < N; i++)
        fu[i][0] = fu[i][1] = 0;

    for (i = 0; i < N; i++){
        for (j = i + 1; j < N; j++){
            calculate_fij(pos[i], pos[j], fij);
            for (k = 0; k < 2; k++){
                fu[i][k] += fij[k];
                fu[j][k] -= fij[k];
            }
        }
    }
}


/*
Algoritmo de Verlet de velocidades, evolución del sistema.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
void verlet_iteration(class_system *sys){
    int i,j;

    for(j = 0; j < 2; j++){
        for (i = 0; i < N; i++){
            sys -> momentum[i][j] += 0.5 * sys -> aceleration[i][j] * sys -> pass_time;
            sys -> position[i][j] += sys -> momentum[i][j] * sys -> pass_time;
        }
    }
  
    sys -> f(sys -> position, sys -> aceleration);

    for(j = 0; j < 2; j++){
        for(i = 0; i < N; i++){
            sys -> momentum[i][j] += 0.5 * sys -> aceleration[i][j] * sys -> pass_time;
        }
    }
  
    sys -> accountant++;
}


/*
Condiciones de contorno periódicas.
Keyword Arguments:
sys: class_system * (sistema de partícula).
*/
void periodic_condition(class_system *sys){
    int i;
    double size = 1.0 / L;

    for(i = 0; i < N; i++){
        sys -> position[i][0] -= floor(sys -> position[i][0] * size) * L;
        sys -> position[i][1] -= floor(sys -> position[i][1] * size) * L;
    }
}


/*
Cálculo de la fuerza entre dos partículas en las posiciones ri y rj.
Keyword Arguments:
ri: vector_R (posición de la primera partícula).
rj: vector_R (posición de la segunda partícula).
fij: vector_R (componentes de la fuerza).
*/
void calculate_fij(vector_R ri, vector_R rj, vector_R fij){
    double dij = dist_ij(ri, rj);
    double x = (pow(dij,-14) - 0.5 * pow(dij,-8));
  
    fij[0] = 48 * minimal_image(ri[0], rj[0]) * x;                                 
    fij[1] = 48 * minimal_image(ri[1], rj[1]) * x;
}


/*
Cálculo de la distancia entre dos partículas en las posiciones ri y rj.
Keyword Arguments:
ri: vector_R (posición de la primera partícula).
rj: vector_R (posición de la segunda partícula).
*/
double dist_ij(vector_R ri, vector_R rj){
    return sqrt(square(minimal_image(ri[0], rj[0])) + square(minimal_image(ri[1], rj[1])));
}


/*
Implementación de método de imagen mínima.
Keyword Arguments:
x_1: double (proyección de la primera partícula).
x_2: double (proyección de la segunda partícula).
*/
double minimal_image(double x_1, double x_2){
    double x_ij = x_1 - x_2;
    return x_ij - L * round(x_ij / L); 
}


/*
Cálculo de la energía cinética instantánea del sistema.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
double energy_K(class_system *sys){ 
    int i;
    double K = 0;
  
    for(i = 0; i < N; i++)
        K += 0.5 * (square(sys -> momentum[i][0])+ square(sys -> momentum[i][1]));

    return K;
}


/*
Calculo de la energía potencial instantánea del sistema.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
double energy_P(class_system *sys){
    int i,j;
    double U = 0, dij, tem;
   
    for(i = 0; i < N; i++){
        for (j = i + 1; j < N; j++){
            dij = dist_ij(sys -> position[i], sys -> position[j]);
            tem = pow(dij,-6);
            U += 4 * tem * (tem - 1);
       }
    }
    
    return U;
}


/*
Generación de los histogramas de módulos de velocidades. Retorna la entropía
del sistema de partículas.
Keyword Arguments:
v: vector_RN (velocidades instantaneas del sistema).
*/
double histogram_system(vector_RN v){
    FILE *histogram_1 = fopen("histogr_1", "w");
    FILE *histogram_2 = fopen("histogr_2", "w");
    
    int index, i, NCajas = 20;
    
    double dx, Vmax, Vmin, Entropia = 0;
    double histo[NCajas];
    double tem, aux[900];

    
 
    for(i = 0; i < N; i++)
        aux[i] = sqrt(square(v[i][0]) + square(v[i][1]));

    for(i = 0; i < NCajas; i++)
        histo[i] = 0;
  
    Vmax = Vmin = aux[0];

    for(i = 0; i < N; i++){
        if(Vmax < aux[i])  Vmax = aux[i];
        if(Vmin > aux[i])  Vmin = aux[i];
    }
    
    dx = (Vmax - Vmin) / NCajas;

    for(i = 0; i < N; i++){
        fprintf(histogram_1, "%lg\n", aux[i]);
        index = (int)((aux[i] - Vmin) / (1.0 * dx));
        histo[index]++;
    }
  
    for(i = 0; i < NCajas; i++){ 
        fprintf(histogram_2, "%lg  %lg\n", Vmin + i * dx * 0.5, histo[i] / N);
        if(histo[i] != 0){
           tem = histo[i]/N;
           Entropia -= dx * tem * log(tem);
        }
    }    
  
    fclose(histogram_1);
    fclose(histogram_2);

    return Entropia;
}

/*
Impresión que permite visualizar la evolución del sistema.
Keyword Arguments:
sys: class_system * (sistema de particulas).
gnuplot_pipe: ----> variable de entorno (Window).
              ----> objeto FILE * (Linux).
*/
void printer_system(class_system *sys, FILE *gnuplot_pipe){
    int i;
   
    fprintf(gnuplot_pipe, "set title '{/=20 Modelo de Gas pass=%d}'\n", sys -> accountant);
    fprintf(gnuplot_pipe, "unset tics \n");   
    fprintf(gnuplot_pipe, "set xlabel  '{/=15 X}'\n");
    fprintf(gnuplot_pipe, "set ylabel  '{/=15 Y}'\n");
    fprintf(gnuplot_pipe, "plot '-' with p pointtype 6 t ''\n");
  
    for(i = 0; i < N; i++)
        fprintf(gnuplot_pipe, "%lf\t%lf\n", sys -> position[i][0], sys -> position[i][1]);
    
    fprintf(gnuplot_pipe, "e\n");
    fflush(gnuplot_pipe);
}


/*
Funcion que demuestra la simetria de inversión del sistema.
*/
void inversor(){
    int i, j;
    class_system *sys = builder_system(0.005, f);

    FILE *feop_1 = fopen("veloc_1","w");
    FILE *feop_2 = fopen("veloc_2","w");

    for(j = 0; j < N; j++)
        fprintf(feop_1, "%lf\t%lf\n", sys -> momentum[j][0], sys -> momentum[j][1]);

    for(i = 0; i < 4000; i++){
        verlet_iteration(sys);
        periodic_condition(sys);

        if(i == 1999){
            for(j = 0; j < N; j++){
                sys -> momentum[j][0] = - sys -> momentum[j][0];
                sys -> momentum[j][1] = - sys -> momentum[j][1];
            }
        }

    }

    for(j = 0; j < N; j++)
    fprintf(feop_2, "%lf\t%lf\n", sys -> momentum[j][0], sys -> momentum[j][1]);

    fclose(feop_1);
    fclose(feop_2);
    destroy_system(sys);
}


/*
Liberación de la memoria solicitada para crear un objeto class_system *.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
void destroy_system(class_system *sys){
    free(sys -> position);
    free(sys -> momentum);
    free(sys -> aceleration);
    free(sys);
}
