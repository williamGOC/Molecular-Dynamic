#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include "generation.c"


#define NC 10
#define RCUT 2.5
#define RHO 0.8
#define DEPS 0.5
#define EMPTY -1
#define FIRST_N 13
#define N 4*(NC)*(NC)*(NC)
#define square(A) (A)*(A)


typedef double vector_R[3];
typedef vector_R *vector_RN;


typedef struct particles{
    vector_RN position, momentum, aceleration;
    double pass_time;
    double temperature;
    double * etta;
    int ** map;
    int accountant;
    long seed;
} class_system;


class_system *builder(double pass_time, double T_0);
double crystalline_parameter(class_system *sys);
double dist_ij(vector_R ri, vector_R rj);
double minimal_image(double x_1, double x_2);
void calculate_fij(vector_R ri, vector_R rj, vector_R fij);
void difusion();
void radial_distri(class_system *sys);
double T(class_system *sys);
void periodic_condition(class_system *sys);
void nose_hoover_iteration(class_system *sys);
void printer_system(class_system *sys, FILE *gnuplot_pipe);
void linked_list(vector_RN r, vector_RN fu, int** map);
void destroy_system(class_system *sys);
int index_cell(int i_x, int i_y, int i_z);
int **neihgbors_map(class_system *sys);



void main(){
    int i, j;
    double par, t;
    int n_scal = 50;
    int NEQ = 1000;
    int NRUM = 1000;
double T_0 = 1.0;
    //class_system *sys = builder(0.005, T_0);
    
    
    srand(time(NULL));
   /* FILE *feop = fopen("parametro.txt","w");
    FILE *pipe = popen("gnuplot", "w");
    assert(pipe);
    
    for (i=0;i<NEQ + NRUM;i++){ 
        nose_hoover_iteration(sys);
        periodic_condition(sys);
        //if(i % n_scal == 0 && i < NEQ) thermals(sys);
        //par = crystalline_parameter(sys);
        //t = T(sys);
        //fprintf(feop, "%lf\t%lf\t%lf\n", i * sys -> pass_time, t, par);
        
        printer_system(sys, pipe);
    }*/
   	//radial_distri(sys);
    difusion();
    //fclose(pipe);
    //fclose(feop);
    //destroy_system(sys);
}


/*
Constructor de un objeto de tipo class_system *. Inicialización de las
posiciones en una estructura FCC y las velocidades con una distribución gaussiana.
Keyword Arguments:
pass_time: double (paso de tiempo del iterador).
f: force (fuerza de interación entre las partículas del sistema).
*/
class_system *builder(double pass_time, double T_0){
    int n_x, n_y, n_z, i, k;
    double L = pow(N / RHO, 1.0 / 3.0);
    double a = L / NC;
    double V_CM[3] = {0.0, 0.0, 0.0};
    
    class_system *sys = (class_system *)malloc(sizeof(class_system));        
    assert(sys);

    sys -> position = (vector_RN)malloc(N * sizeof(vector_R));
    assert(sys -> position);
    sys -> momentum = (vector_RN)malloc(N * sizeof(vector_R));
    assert(sys -> momentum);
    sys -> aceleration = (vector_RN)malloc(N * sizeof(vector_R));
    assert(sys -> aceleration);
    sys -> etta = (double *)malloc(sizeof(double));
    assert(sys -> etta);
  
    sys -> pass_time = pass_time;
    sys -> temperature = T_0;
    sys -> accountant = 0;
    sys -> seed = time(NULL);
    sys -> map = neihgbors_map(sys);
    *sys -> etta = 0.0;

    for(i = 0; i < N; i += 4){
        n_x = (i / 4) % NC;
        n_y = (i / (4 * NC)) % NC;
        n_z = (i / (4 * NC * NC)) % NC;

        sys -> position[i][0] = n_x * a;
        sys -> position[i][1] = n_y * a;
        sys -> position[i][2] = n_z * a;

        sys -> position[i + 1][0] = n_x * a;
        sys -> position[i + 1][1] = (n_y + DEPS) * a;
        sys -> position[i + 1][2] = (n_z + DEPS) * a;

        sys -> position[i + 2][0] = (n_x + DEPS) * a;
        sys -> position[i + 2][1] = n_y * a;
        sys -> position[i + 2][2] = (n_z + DEPS) * a;

        sys -> position[i + 3][0] = (n_x + DEPS) * a;
        sys -> position[i + 3][1] = (n_y + DEPS) * a;
        sys -> position[i + 3][2] = n_z * a;
    }

    for(i = 0; i < N; i++){
        for(k = 0; k < 3; k++){
            sys -> momentum[i][k] = sqrt(T_0) * gasdev(&sys -> seed);
            V_CM[k] += sys -> momentum[i][k] / (N);
        }
    }

    for(i = 0; i < N; i++){
        for(k = 0; k < 3; k++)
            sys -> momentum[i][k] -= V_CM[k];
        
    }

    linked_list(sys -> position, sys -> aceleration, sys -> map);
   
    return sys;
}


/*
Cálculo de la temperatura instantánea del sistema.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
double T(class_system *sys){
    int i;
    double K = 0;
  
    for(i = 0; i < N; i++){
        K += 0.5 * (square(sys -> momentum[i][0])+ square(sys -> momentum[i][1]) 
                + square(sys -> momentum[i][2]));
    }
    return K / (3 * N);
}


/*
Algoritmo de Verlet (Nosé-Hoover) de velocidades, evolución del sistema.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
void nose_hoover_iteration(class_system *sys){
    int i,j;
    double g = 0, aux;
    int NF = 3 * N - 3;
    double tau = (double)(1 / sqrt(200));

    for(j = 0; j < 3; j++){
        for (i = 0; i < N; i++){
            sys -> momentum[i][j] += 0.5 * sys -> pass_time * (sys -> aceleration[i][j] - (*sys -> etta) * sys -> momentum[i][j]);
            sys -> position[i][j] += sys -> momentum[i][j] * sys -> pass_time;
        }
    }
    

    linked_list(sys -> position, sys -> aceleration, sys -> map);
    for(i = 0; i < N; i++){
        for(j = 0; j < 3; j++)
            g += (sys -> momentum[i][j] * sys -> momentum[i][j]) / (NF * sys -> temperature);
    }

    g = (g - 1) / (tau * tau); 
    *sys -> etta += sys -> pass_time * g;

    for(j = 0; j < 3; j++){
        for(i = 0; i < N; i++){
            aux = 1.0 + 0.5 * sys -> pass_time * (*sys -> etta); 
            sys -> momentum[i][j] = (0.5 * sys -> aceleration[i][j] * sys -> pass_time + sys -> momentum[i][j]) / aux;
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
    double L = pow(N / RHO, 1.0 / 3.0);
    double size = 1.0 / L;

    for(i = 0; i < N; i++){
        sys -> position[i][0] -= floor(sys -> position[i][0] * size) * L;
        sys -> position[i][1] -= floor(sys -> position[i][1] * size) * L;
        sys -> position[i][2] -= floor(sys -> position[i][2] * size) * L;
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
    fij[2] = 48 * minimal_image(ri[2], rj[2]) * x;
}


/*
Cálculo de la distancia entre dos partículas en las posiciones ri y rj.
Keyword Arguments:
ri: vector_R (posición de la primera partícula).
rj: vector_R (posición de la segunda partícula).
*/
double dist_ij(vector_R ri, vector_R rj){
    return sqrt(square(minimal_image(ri[0], rj[0])) + square(minimal_image(ri[1], rj[1])) 
            + square(minimal_image(ri[2], rj[2])));
}


/*
Implementación de método de imagen mínima.
Keyword Arguments:
x_1: double (proyección de la primera partícula).
x_2: double (proyección de la segunda partícula).
*/
double minimal_image(double x_1, double x_2){
    double x_ij = x_1 - x_2;
    double L = pow(N / RHO, 1.0 / 3.0);
    return x_ij - L * round(x_ij / L); 
}


/*
Parámetro cristalino del sistema de partículas.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
double crystalline_parameter(class_system *sys){
    double aux, x_1 = 0, x_2 = 0;
    int i;
    double L = pow(N / RHO, 1.0 / 3.0);
    double a = L / NC, value;
    
    aux = (2*M_PI)/a;
    
    for(i = 0; i < N; i++){
        x_1 += sin((aux * (sys -> position[i][1] - sys -> position[i][0]- sys -> position[i][2])));
        x_2 += cos((aux * (sys -> position[i][1] - sys -> position[i][0]- sys -> position[i][2])));
    }
    value = (x_1*x_1 + x_2*x_2)/(N*N);
    return value;
}

/*
Parámetro de difusión del sistema de partículas
*/
void difusion(){
    int n_scal = 50;
    int NEQ = 1000;
    int NRUM = 1000;
    int i, j, k;
    double aux, cont = 0;

    class_system *sys;
    sys = builder(0.005, 1.0);
    
    FILE * folder = fopen("difusion.txt","w");
    
    vector_RN temp = (vector_RN)malloc(N * sizeof(vector_R));
    assert(temp);

    for (i = 0; i < NEQ + NRUM; i++){ 
        nose_hoover_iteration(sys);
        periodic_condition(sys);
        //if(i % n_scal == 0 && i < NEQ) thermals(sys);
        
        if(i == 0){
            for(k = 0; k < N; k++){
                for(j = 0; j < 3; j++)
                    temp[k][j] = sys->position[k][j];
            }
        }

        if(i > 0){
            for(k = 0; k < N; k++){
                aux = dist_ij(sys->position[k], temp[k]);
                cont += (aux * aux);
            }
        }
       
        cont = cont / 6;
        fprintf(folder, "%lf\t%lf\n", i*sys->pass_time, cont);
        cont = 0;
    }
    free(temp);
    fclose(folder);
}



/*
Función de distribución radial en equilibrio. Genera un histograma promediado en el
tiempo de las distancias entre partículas del sistema y lo normaliza respecto al del
gas ideal.
Keyword Arguments:
sys: class_system * (sistema de partículas).
*/
void radial_distri(class_system *sys){
    int i, ii, j, k, n_scal = 50;
    double cont = 0;
    double dij, delt_r = 0.1;
    double L = pow(N / RHO, 1.0/3.0);
    int CAJAS = (int)(L / delt_r);
    int NUM = 2000;
    
    FILE * folder = fopen("radial.txt", "w");
    double **generator = (double **)malloc(NUM * sizeof(double *));
    assert(generator);

    for(i = 0; i < NUM; i++){
        generator[i] = (double *)malloc(CAJAS * sizeof(double));
        assert(generator[i]);
    }

    for(i =0 ; i < NUM; i ++){
        for(j = 0; j < CAJAS; j++)
            generator[i][j] = 0.0;
    }
    
    for(i = 0; i < NUM; i++){
        nose_hoover_iteration(sys);
        periodic_condition(sys);
        //if(i % n_scal == 0 && i < NUM / 2) thermals(sys);
        
        if(i > NUM / 2){
            for(ii = 0; ii < N; ii ++){
                for(j = ii + 1; j < N; j ++){
                    dij = dist_ij(sys -> position[ii], sys -> position[j]);
                    k = floor(dij / delt_r);
                    generator[i][k] += 2;
                }
            }
        }
    }
    
    for(i =0 ; i < CAJAS; i++){
        for(j = 0; j < NUM; j++){
            cont += (2 * generator[j][i]) / (N * NUM * delt_r);
        }
        fprintf(folder, "%lf\t%lf\n", (i + 1) * delt_r, (double)(cont/(4 * M_PI * RHO * delt_r *delt_r * (i+1)*(1+i))));
        cont = 0;
    }

    for(i = 0; i < NUM; i++)
        free(generator[i]);
    free(generator);
    fclose(folder);
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
   
    fprintf(gnuplot_pipe, "set title '{/=20 Modelo de Gas, pass %d}'\n", sys->accountant);
    fprintf(gnuplot_pipe, "set xlabel  '{/=15 X}'\n");  
    fprintf(gnuplot_pipe, "set zlabel  '{/=15 Z}'\n");
    fprintf(gnuplot_pipe, "set ylabel  '{/=15 Y}'\n");
    fprintf(gnuplot_pipe, "splot '-' with p pointtype 6 t ''\n");
  
    for(i = 0; i < N; i++)
        fprintf(gnuplot_pipe, "%lf\t%lf\t%lf\n", sys -> position[i][0], 
                                    sys -> position[i][1], sys ->position[i][2]);
    
    fprintf(gnuplot_pipe, "e\n");
    fflush(gnuplot_pipe);
}


/*
Equivalencia entre índices de partículas del sistema y lista de vecinos.
Keyword Arguments:
i_x: int (indice del eje Ox).
i_y: int (indice del eje Oy).
i_z: int (indice del eje Oz).
*/
int index_cell(int i_x, int i_y, int i_z){
    double L = pow(N / RHO, 1.0 / 3.0);
    int M = (int)(L / RCUT);

    return (i_x + M) % M + ((i_y + M) % M) * M + ((i_z + M) % M) * M *M;
}


/*
Creación de la lista de vecinos.
Keyword Arguments:
sys: class_system * (sistema de particulas).
*/
int **neihgbors_map(class_system *sys){
    double L = pow(N / RHO, 1.0 / 3.0);
    int M = (int)(L / RCUT);
    int NCEL = M * M * M;
    double L_cell = L / M;
    int i, i_x, i_y, i_z;

    int ** map = (int **)malloc(NCEL * sizeof(int *));
    assert(map);

    for(i = 0; i < NCEL; i++){
        map[i] = (int *)malloc(FIRST_N * sizeof(int));
        assert(map[i]);
    }

    for(i_x = 0; i_x < M; i_x++){
        for(i_y = 0; i_y < M; i_y++){
            for(i_z = 0; i_z <M; i_z++){
                i = index_cell(i_x, i_y, i_z);
                map[i][0] = index_cell(i_x + 1, i_y, i_z);
                map[i][1] = index_cell(i_x, i_y + 1, i_z);
                map[i][2] = index_cell(i_x + 1, i_y + 1, i_z);
                map[i][3] = index_cell(i_x + 1, i_y - 1, i_z);
                map[i][4] = index_cell(i_x, i_y, i_z + 1);
                map[i][5] = index_cell(i_x + 1, i_y, i_z + 1);
                map[i][6] = index_cell(i_x + 1, i_y, i_z - 1);
                map[i][7] = index_cell(i_x, i_y + 1, i_z + 1);
                map[i][8] = index_cell(i_x, i_y + 1, i_z - 1);
                map[i][9] = index_cell(i_x + 1, i_y + 1, i_z + 1);
                map[i][10] = index_cell(i_x + 1, i_y - 1, i_z + 1);
                map[i][11] = index_cell(i_x + 1, i_y + 1, i_z - 1);
                map[i][12] = index_cell(i_x + 1, i_y - 1, i_z - 1);
            }
        }
    }

    return map;
}


/*
Cálculo de las fuerzas entre las partículas del sistema mediante el algoritmo linked-cell.
Keyword Arguments:
sys: class_system * (sistema de particulas).
map: int ** (mapa de vecinos).
*/
void linked_list(vector_RN r, vector_RN fu, int** map){
    int  i, icell, ivec, j, k, jcell;
    double L = pow(N / RHO, 1.0 / 3.0);
    int M = (int)(L / RCUT);
    int NCEL = M * M * M;
    double L_cell = L / M;
    double x, y, z, rij;
    
    vector_R fij;

    int head[NCEL];
    int list[N];

    for(icell = 0; icell < NCEL; icell++)
        head[icell] = EMPTY;
    
    for(i = 0; i < N; i++){
        x = r[i][0] - L * floor(r[i][0] / L);
        y = r[i][1] - L * floor(r[i][1] / L);
        z = r[i][2] - L * floor(r[i][2] / L);
        icell = (int)(x / L_cell) + ((int)(y / L_cell)) * M + ((int)(z / L_cell)) * M *M;
        list[i] = head[icell];
        head[icell] = i;
    }

    for(i = 0; i < N; i++)
        fu[i][0] = fu[i][1] = fu[i][2] = 0;

    for(icell = 0; icell < NCEL; icell++){
        i = head[icell];
        while(i >= 0){
            j = list[i];
            while(j >= 0){
                rij = dist_ij(r[i], r[j]);
                if(square(rij) < square(RCUT)){
                    calculate_fij(r[i], r[j], fij);
                    for (k = 0; k < 3; k++){
                        fu[i][k] += fij[k];
                        fu[j][k] -= fij[k];
                    }
                }
                
                j = list[j];
            }    
            
            for(ivec = 0; ivec < FIRST_N; ivec++){
                jcell = map[icell][ivec];
                j = head[jcell];

                while(j >= 0){
                    rij = dist_ij(r[i], r[j]);
                    if(square(rij) < square(RCUT)){
                        calculate_fij(r[i], r[j], fij);
                        for (k = 0; k < 3; k++){
                            fu[i][k] += fij[k];
                            fu[j][k] -= fij[k];
                        }
                    }
                            
                    j = list[j];
                }
            }
            
            i = list[i];
        }
    }
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
    free(sys -> etta);
    free(sys);
}
