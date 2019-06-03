//To compile should add "-lfftw3 -lgsl"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>
#include <fftw3.h>
#include <complex.h>
using namespace std;

#include "Grid.h"
#include "Function.h"

int main( int argc, char *argv[] ){
  //================Simulation Constants
  int weightFunction = 2;   //0/1/2 : NGP/CIC/TSC
  int dim = 2;        
  double L = 10.0;        //Length of box (from -L/2 ~ L/2)
  int Nx = 11;        //Number of grid in x direction. (should be odd number)
  int NParticle=2;//Number of particles used in simulation
  double massParticle=1.0;

  //================Random number generator.
    //To use : d=gsl_rng_uniform(rng);
    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,123456);//The seed is 123456.

  //================Initialize Grid Parameter=========== 
    struct grid2D grid;
    grid.L = L;   
    grid.Nx = Nx;   
    grid.Ny = grid.Nx;
    grid.N = grid.Nx * grid.Ny;
    grid.dx = grid.L / (grid.Nx-1); //-1 is because boundary (make box closed)
    grid.dy = grid.L / (grid.Ny-1);
    grid.density = (double*) malloc(grid.N*sizeof(double)); //Density on the grid
    grid.phi = (double*) malloc(grid.N*sizeof(double));   //Potential on the grid
    grid.Fx = (double*) malloc(grid.N*sizeof(double));      //Force on grid comes from potential calculation
    grid.Fy = (double*) malloc(grid.N*sizeof(double));

  //================Initialize Particles Parameter======
    struct particle2D myParticle;
    myParticle.number = NParticle;
    myParticle.mass = (double*)malloc(myParticle.number*sizeof(double));
    myParticle.x = (double*)malloc(myParticle.number*sizeof(double));
    myParticle.y = (double*)malloc(myParticle.number*sizeof(double));
    myParticle.Fx = (double*)malloc(myParticle.number*sizeof(double));
    myParticle.Fy = (double*)malloc(myParticle.number*sizeof(double));
    myParticle.vx = (double*)malloc(myParticle.number*sizeof(double));
    myParticle.vy = (double*)malloc(myParticle.number*sizeof(double));
  
  //================Initialize density field=============
    for(int i=0;i<grid.N;i++){
      grid.density[i]=0.0;    
    }

  //Initialize Initial Position of Particle.
    for (int i = 0; i < myParticle.number; ++i){
      myParticle.mass[i]=1.0;
      myParticle.x[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
      myParticle.y[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
      printf("At (%f,%f) \n",myParticle.x[i],myParticle.y[i]);
    }

  //Initialize Force field value
    for(int i=0;i<grid.N;i++){
      grid.Fx[i]=0.0;
      grid.Fy[i]=0.0;   
    }

  //Deposit Particles to grid
  Weight2d(&grid,&myParticle,weightFunction);
  
  //Print out Density field
    for(int i=0;i<grid.N;i++){
      if(i % grid.Nx == 0){
        printf("\n");
      }
      printf("%.2f\t",grid.density[i]);
    }
    printf("\n");

  //Use Fourier Transform to calculate potential and force.
  poisson_solver_fft_force_2d(dim, &grid);

  //Remap the force to particle.
  WeightForce2d(&grid,&myParticle,weightFunction);

  
  //Print out the force on a particle.
    for(int i=0;i<myParticle.number;i++){
      printf("(Fx,Fy)=(%.2f,%.2f)\n",myParticle.Fx[i],myParticle.Fy[i]);
    }

  return 0;
}
