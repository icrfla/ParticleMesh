
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string>
#include <iostream>
#include <time.h>
#include <fftw3.h>
#include <complex.h>

#include "Grid.h"
#include "Function.h"

using namespace std;

void poisson_solver_fft_force_3d(int const dim, struct grid3D *grid){

  // double G_const = 6.67408e-8; // #g^-1 s^-2 cm^3 
  double G_const = 1.; // #g^-1 s^-2 cm^3 

  double dNx = (double) grid->Nx;
  double dNy = (double) grid->Ny;
  double dNz = (double) grid->Nz;

  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  int const Nz = grid->Nz;

  // int total_n = Nx*Ny*Nz;

  fftw_complex *sigma_a, *fftsigma_a;
  fftw_complex *phika, *phia; 
  fftw_plan p, q;

  sigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  fftsigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  phika = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  phia = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);

  /////////// initialization of sigma_a ///////////
  int ii, jj, kk, index;
  for (ii=0; ii<Nx; ii+=1) {
    for (jj=0; jj<Ny; jj+=1) {
      for (kk=0; kk<Nz; kk+=1){
        index = ii+ jj*Nx + kk*Nx*Ny;
        sigma_a[ index ][0] = grid->density[ index ];
        sigma_a[ index ][1] = 0; 
      }
    }
  }


  printf("fft\n");
  /////////// fft ///////////
  p = fftw_plan_dft_3d(Nx, Ny, Nz, sigma_a, fftsigma_a, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  double kxx, kyy, kzz;
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      for (kk=0; kk<Nz; kk+=1){
        kxx = pow((double)(ii+1)*2.*M_PI/dNx, 2.);
        kyy = pow((double)(jj+1)*2.*M_PI/dNy, 2.);
        kzz = pow((double)(kk+1)*2.*M_PI/dNz, 2.);
        index = ii+ jj*Nx + kk*Nx*Ny;
        phika[ index ][0] = -4. * M_PI * G_const * fftsigma_a[ index ][0] / ((kxx+kyy+kzz));
        phika[ index ][1] = -4. * M_PI * G_const * fftsigma_a[ index ][1] / ((kxx+kyy+kzz));
      }
    }
  }

  printf("inverse fft\n" );
  /////////// inverse fft ///////////
  q = fftw_plan_dft_3d(Nx, Ny, Nz, phika, phia, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(q);
  fftw_destroy_plan(q);
  /////////// normalization ///////////
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      for (kk=0; kk<Nz; kk+=1){
        index = ii+ jj*Nx + kk*Nx*Ny;
        phia[ index ][0] /= grid->N;
        phia[ index ][1] /= grid->N;
      }
    }
  }

  // for (ii=0; ii<Nx; ii+=1){
  //   for (jj=0; jj<Ny; jj+=1){
  //     for (kk=0; kk<Nz; kk+=1){
  //       index = ii*Nz*Ny + jj*Nz + kk;
  //       printf("recover: %3d %3d %3d %+.5e %+.5e I vs. %+.5e %+.5e I \n",
  //           ii, jj, kk, sigma_a[ index ][0], sigma_a[ index ][1], 
  //           phia[ index ][0], phia[ index ][1]);      
  //     }
  //   } 
  // }

  printf("phi magnitude\n");
  double *phi_total = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  ////////// magnitude of potential //////////
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){ 
      for (kk=0; kk<Nz; kk+=1){
        index = ii+ jj*Nx + kk*Nx*Ny;
        printf("index: %d, phia[0]: %.3e, phia[1]: %.3e\n", index, phia[index][0], phia[index][1]);
        grid->phi[ index ] = sqrt( pow(phia[ index ][0],2) + pow(phia[ index ][1],2) );        
      }
    }
  } 

  printf("calculate force\n");
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      for (kk=0; kk<Nz; kk+=1){
        _2nd_order_diff_3d( grid, ii, jj, kk );
      }
    }
  }
    
  fftw_cleanup();

}
