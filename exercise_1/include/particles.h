#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#include "config.h"

typedef void (* force)(double *, double *);

typedef struct particles {

    size_t memoryX;
    size_t memoryV;
    size_t memoryA;

    double *X;
    double *V;
    double *A;
    
    double dt;
    double t;
    
    force f;

} classParticles;


classParticles *initSystem(double, force);
void freeSystem(classParticles *);



#endif // __PARTICLES_H__