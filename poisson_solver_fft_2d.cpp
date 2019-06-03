
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

void poisson_solver_fft_force_2d(int const dim, struct grid2D *grid) {

  // double G_const = 6.67408e-8; // #g^-1 s^-2 cm^3 
  double G_const = 1.; // #g^-1 s^-2 cm^3 

  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  // int total_n = Nx*Ny;
  double dNx = (double) (grid->Nx), dNy = (double) (grid->Ny); // default Nx=Ny
  int ii, jj, index;

  fftw_complex *sigma_a, *fftsigma_a;
  fftw_complex *phika, *phia; 
  fftw_plan p, q;

  sigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
  fftsigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
  phika = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
  phia = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);

  for (ii=0; ii<Nx; ii+=1) {
    for (jj=0; jj<Ny; jj+=1) {
      index = ii + jj*Nx;
      sigma_a[ index ][0] = grid->density[ index ]; 
      sigma_a[ index ][1] = 0.;
    }
  }

  /////////// fft ///////////
  p = fftw_plan_dft_2d(Nx, Ny, sigma_a, fftsigma_a, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  double kxx, kyy;
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      index = ii + jj*Nx;
      kxx = pow((double)(ii+1)*2.*M_PI/dNx, 2.);
      kyy = pow((double)(jj+1)*2.*M_PI/dNy, 2.);
      phika[index][0] = -4. * M_PI * G_const * fftsigma_a[index][0] / ((kxx+kyy));
      phika[index][1] = -4. * M_PI * G_const * fftsigma_a[index][1] / ((kxx+kyy));
    }
  }
  /////////// inverse fft ///////////
  q = fftw_plan_dft_2d(Nx, Ny, phika, phia, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(q);
  fftw_destroy_plan(q);
  /////////// normalization ///////////
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      index = ii + jj*Nx;
      phia[index][0] /= grid->N;
      phia[index][1] /= grid->N;
    }
  }

  double *phi_total = (double *) malloc(Nx*Ny * sizeof(double));
  ////////// magnitude of potential //////////
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){ 
      index = ii + jj*Nx;
      grid->phi[index] = sqrt( pow(phia[index][0],2) + pow(phia[index][1],2) );
    }
  }  

  ////////// calculate force based on potential //////////
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){ 
      _2nd_order_diff_2d(grid, ii, jj);
    }
  }

  fftw_cleanup();

}

