#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string>
#include <iostream>
#include <time.h>
#include <fftw3.h>
#include <complex.h>
#include <gsl/gsl_rng.h>

#include "Grid.h"
#include "Function.h"

using namespace std;


int main(){
  int weightFunction = 1;  //0/1/2 : NGP/CIC/TSC
  int dim = 3;
  //================Random number generator.
  //To use : d=gsl_rng_uniform(rng);
  gsl_rng *rng;
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,123456);//The seed is 123456.

  //================Grid 
  struct grid3D grid;
  grid.L = 10.0;    //Length of box (-L/2 ~ L/2)
  grid.Nx = 11;   //Number of grid in x direction. (Should be an odd number)
  grid.Ny = grid.Nx;
  grid.Nz = grid.Nx;
  grid.N = grid.Nx * grid.Ny * grid.Nz;
  grid.dx = grid.L / (grid.Nx-1); //-1 is because boundary (make box closed)
  grid.dy = grid.L / (grid.Ny-1);
  grid.dz = grid.L / (grid.Nz-1);
  grid.density = (double*) malloc(grid.N*sizeof(double)); //Density on the grid
  grid.phi = (double*) malloc(grid.N*sizeof(double));
  grid.Fx = (double*) malloc(grid.N*sizeof(double));      //Force on grid comes from potential calculation
  grid.Fy = (double*) malloc(grid.N*sizeof(double));
  grid.Fz = (double*) malloc(grid.N*sizeof(double));

  

  //================Particles
  int NParticle=10;//Number of particles used in simulation
  double massParticle=1.0;
    
  struct particle3D myParticle;
  myParticle.number = NParticle;
  myParticle.mass = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.x = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.y = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.z = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.Fx = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.Fy = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.Fz = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.vx = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.vy = (double*)malloc(myParticle.number*sizeof(double));
  myParticle.vz = (double*)malloc(myParticle.number*sizeof(double));
  

  //Initialize density grid
  for(int i=0;i<grid.N;i++){
    grid.density[i]=0.0;    
  }

  //Initialize Initial Position of Particle.
  for (int i = 0; i < myParticle.number; ++i){
    myParticle.mass[i]=1.0;
    myParticle.x[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
    myParticle.y[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
    myParticle.z[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
    printf("At (%f,%f) \n",myParticle.x[i],myParticle.y[i]);
  }

  // //Initialize Force field (for test use)
  // for(int i=0;i<grid.N;i++){
  //   grid.Fx[i]=0.0;
  //   grid.Fy[i]=0.0;
  //   grid.Fz[i]=0.0;
    
  // }
  // grid.Fx[3+6*grid.Nx] = -1.0;
  // grid.Fy[3+6*grid.Nx] = 2.0;

  Weight3d(&grid,&myParticle,weightFunction);

  // for(int i=0;i < grid.Nx * grid.Ny;i++){
  //   if(i % grid.Nx == 0){
  //     printf("\n");
  //   }
  //   printf("%.2f\t",grid.density[i]);
  // }
  // printf("\n");


  //Use Fourier Transform to calculate potential and force.
  poisson_solver_fft_force_3d(dim, &grid);

  //Remap the force to particle.
  WeightForce3d(&grid,&myParticle,weightFunction);

  
  //Print out the force on a particle.
  for(int i=0;i<myParticle.number;i++){
    printf("(Fx,Fy,Fz)=(%.2f,%.2f,%.2f)\n",myParticle.Fx[i],myParticle.Fy[i],myParticle.Fz[i]);
  }



  return 0;
}