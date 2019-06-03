// g++-9 poisson_solver_fft_2d.cpp -L/usr/local/Cellar/fftw/3.3.8_1/lib -I/usr/local/Cellar/fftw/3.3.8_1/include -lfftw3

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

// void poisson_solver_fft_force_2d(int const dim, struct grid2D *grid);

// void _2nd_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

// void _4th_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

// void _6th_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

// int main( int argc, char *argv[] ){

//   int canvas_size = 20, Nx = 20, Ny = 20;
//   int dim = 2;
//   int thread_num = 2;
//   double dx = 1., dy = 1.;


//   double * sigma_a = (double *) malloc(Nx*Ny * sizeof(double));
//   /////////// initialization of sigma_a ///////////
//   double r2;
//   int ii, jj, lb = -1, rb = 1;
//   for (ii=lb; ii<=rb; ii+=1) {
//     for (jj=lb; jj<=rb; jj+=1) {
//       r2 = (double) (pow(ii,2)+pow(jj,2));
//       printf("ii:%d jj:%d\n", ii, jj);
//       sigma_a[ (Nx/2+ii)*Ny + (Ny/2+jj) ] = 1. * sqrt(1. - r2/1600.);
//     }
//   }

//   double *Fx = (double *) malloc(Nx*Ny * sizeof(double));
//   double *Fy = (double *) malloc(Nx*Ny * sizeof(double));

//   poisson_solver_force_2d(dim, grid);

//   return 0;
// }


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

  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){ 
      _2nd_order_diff_2d(grid, ii, jj);
    }
  }

  // for (ii=0; ii<Nx; ii+=1) {
  //   for (jj=0; jj<Ny; jj+=1){
  //     index = ii + jj*Nx;
  //     printf("recover: %3d %3d %+.5e %+.5e I vs. %+.5e %+.5e I \n",
  //         ii, jj, sigma_a[index][0], sigma_a[index][1], 
  //         phia[index][0], phia[index][1]);    
  //   }
  // }

  // for (ii=0; ii<Nx; ii+=1) {
  //   for (jj=0; jj<Ny; jj+=1){
  //     index = ii + jj*Nx;
  //     printf("ii:%3d, jj:%3d, Fx:%.5e Fy:%.5e\n",
  //         ii, jj, Fx[index], Fy[index]);    
  //   }
  // }


  fftw_cleanup();

}

// void _2nd_order_diff_2d(struct grid2D *grid, int const ii, int const jj ) {

//   double factor1 = -1./(2.*grid->dx);
//   double factor2 = -1./(2.*grid->dy);
//   int const Nx = grid->Nx;
//   int const Ny = grid->Ny;
//   int index = ii + jj*Nx;

//   grid->Fx[ index ] = factor1*( grid->phi[ ( (Nx+ii+1)%Nx ) + jj*Nx ] - grid->phi[ ( (Nx+ii-1)%Nx ) + jj*Nx ] );
//   grid->Fy[ index ] = factor2*( grid->phi[ ii + ((Ny+jj+1)%Ny)*Nx ] - grid->phi[ ii + ((Ny+jj-1)%Ny)*Nx ] );  

// }

// void _4th_order_diff_2d(struct grid2D *grid, int const ii, int const jj ) {

//   double factor1 = -1./(12.*grid->dx);
//   double factor2 = -1./(12.*grid->dy);
//   int const Nx = grid->Nx;
//   int const Ny = grid->Ny;
//   int index = ii + jj*Nx;

//   grid->Fx[ index ] = factor1*( 
//     +1.*grid->phi[ ( (Nx+ii-2)%Nx ) + jj*Nx ] 
//     -8.*grid->phi[ ( (Nx+ii-1)%Nx ) + jj*Nx ] 
//     +8.*grid->phi[ ( (Nx+ii+1)%Nx ) + jj*Nx ]
//     -1.*grid->phi[ ( (Nx+ii+2)%Nx ) + jj*Nx ]
//     );
  

//   grid->Fy[ index ] = factor2*( 
//     +1.*grid->phi[ ii + ((Ny+jj-2)%Ny)*Nx ] 
//     -8.*grid->phi[ ii + ((Ny+jj-1)%Ny)*Nx ] 
//     +8.*grid->phi[ ii + ((Ny+jj+1)%Ny)*Nx ] 
//     -1.*grid->phi[ ii + ((Ny+jj+2)%Ny)*Nx ] 
//     );

// }


// void _6th_order_diff_2d(struct grid2D *grid, int const ii, int const jj ) {

//   double factor1 = -1./(60.*grid->dx);
//   double factor2 = -1./(60.*grid->dy);
//   int const Nx = grid->Nx;
//   int const Ny = grid->Ny;  
//   int index = ii + jj*Nx;

//   grid->Fx[ index ] = factor1*( 
//     -1.* grid->phi[ ( (Nx+ii-3)%Nx )*Nx + jj ]
//     +9.* grid->phi[ ( (Nx+ii-2)%Nx )*Nx + jj ] 
//     -45.*grid->phi[ ( (Nx+ii-1)%Nx )*Nx + jj ] 
//     +45.*grid->phi[ ( (Nx+ii+1)%Nx )*Nx + jj ]
//     -9.* grid->phi[ ( (Nx+ii+2)%Nx )*Nx + jj ]
//     +1.* grid->phi[ ( (Nx+ii+3)%Nx )*Nx + jj ]
//     );
  

//   grid->Fy[ index ] = factor2*( 
//     -1.* grid->phi[ ii + ((Ny+jj-3)%Ny)*Nx ] 
//     +9.* grid->phi[ ii + ((Ny+jj-2)%Ny)*Nx ] 
//     -45.*grid->phi[ ii + ((Ny+jj-1)%Ny)*Nx ] 
//     +45.*grid->phi[ ii + ((Ny+jj+1)%Ny)*Nx ] 
//     -9.* grid->phi[ ii + ((Ny+jj+2)%Ny)*Nx ] 
//     +1.* grid->phi[ ii + ((Ny+jj+3)%Ny)*Nx ] 
//     );

// }
