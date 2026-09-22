#include "particles.h"

classParticles *initSystem(double dt, force f) {

    classParticles *sys = (classParticles *)malloc(sizeof(classParticles));        
    assert(sys);
    
    sys -> memoryX = DIM * N * sizeof(double);
    sys -> memoryV = DIM * N * sizeof(double);
    sys -> memoryA = DIM * N * sizeof(double);

    sys -> X = (double *)malloc(sys -> memoryX);
    assert(sys -> X != NULL);

    sys -> V = (double *)malloc(sys -> memoryV);
    assert(sys -> V != NULL);

    sys -> A = (double *)malloc(sys -> memoryA);
    assert(sys -> A != NULL);

    sys -> f = f;
    sys -> dt = dt;
    sys -> t = 0;

    for (int i = 0; i < N; i++){

        sys -> X[DIM * i + 0] = (i % NC) * L / (NC + 1);
        sys -> X[DIM * i + 1] = (i / NC) * L / (NC + 1);
        
        sys -> V[DIM * i + 1] = ((rand() % 2) ? 1.1 : -1.1);
        sys -> V[DIM * i + 1] = 0.0;
    }

    f(sys -> X, sys -> A);

    return sys;
}


void f(double *X, double *F){
    
    double *Fij = (double *)malloc(DIM * sizeof(double));
    assert(Fij != NULL);

    for(int i = 0; i < N; i++) {

        F[DIM * i + 0] = 0.0;
        F[DIM * i + 1] = 0.0;
    }

    for (int i = 0; i < N; i++){
        for (int j = i + 1; j < N; j++){

            computeFij(X, i, j, Fij);

            for (int dir = 0; dir < DIM; dir++){
                F[DIM * i + dir] += Fij[dir];
                F[DIM * j + dir] -= Fij[dir];
            }
        }
    }

    free(Fij);
}


void computeFij(double *X, int i, int j, double *Fij){
    
    // Compute the distance between particles i and j using the 
    // minimal image principe

    double xi_0 = X[DIM * i + 0];
    double xi_1 = X[DIM * i + 1];

    double xj_0 = X[DIM * j + 0];
    double xj_1 = X[DIM * j + 1];

    double minXij_0 = minImage(xi_0, xj_0);
    double minXij_1 = minImage(xi_1, xj_1);

    double dij = sqrt(minXij_0 * minXij_0 + minXij_1 * minXij_1);

    double x = pow(dij,-14) - 0.5 * pow(dij,-8);
  
    Fij[0] = 48 * minXij_0 * x;                                 
    Fij[1] = 48 * minXij_1 * x;
}


double minImage(double x_0, double x_1){
    double x_ij = x_0 - x_1;
    return x_ij - L * round(x_ij / L); 
}