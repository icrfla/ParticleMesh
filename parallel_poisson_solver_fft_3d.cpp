// g++-9 parallel_poisson_solver_fft_3d.cpp -L/usr/local/Cellar/fftw/3.3.8_1/lib -I/usr/local/Cellar/fftw/3.3.8_1/include -lfftw3 -fopenmp

#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include <string>
#include <iostream>
#include <time.h>
#include <fftw3.h>
#include <complex.h>
using namespace std;

double _fftw_3d_(int const thread_num, int const dim, int const Nx, int const Ny, int const Nz, 
  double const ori_sigma_a[], double* Fx, double *Fy, double *Fz, 
  double const dx, double const dy, double const dz);

double _2nd_order_diff_(int const Nx, int const Ny, int const Nz, double const phia[], 
  double* Fx, double* Fy, double *Fz, double const dx, double const dy, double const dz,
  int const ii, int const jj, int const kk );

int main( int argc, char *argv[] ){

  int canvas_size = 20, Nx = 20, Ny = 20, Nz = 20;
  int dim = 3;
  double dx = 1., dy = 1., dz = 1.;
  int thread_num = 2;

  double * sigma_a = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  /////////// initialization of sigma_a ///////////
  double r2;
  int ii, jj, kk, lb = -1, rb = 1, index;
  for (ii=lb; ii<=rb; ii+=1) {
    for (jj=lb; jj<=rb; jj+=1) {
      for (kk=lb; kk<=rb; kk+=1){
        r2 = (double) (pow(ii,2)+pow(jj,2)+pow(kk,2));
        // printf("ii:%d jj:%d kk:%d\n", ii, jj, kk);
        index = (Nx/2+ii)*Nz*Ny + (Ny/2+jj)*Nz + kk;
        sigma_a[ index ] = 1. * sqrt(1. - r2/1600.);
      }
    }
  }

  double *Fx = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  double *Fy = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  double *Fz = (double *) malloc(Nx*Ny*Nz * sizeof(double));

  _fftw_3d_(thread_num, dim, Nx, Ny, Nz, sigma_a, Fx, Fy, Fz, dx, dy, dz);

  return 0;
}

double _2nd_order_diff_(int const Nx, int const Ny, int const Nz, double const phia[], 
  double* Fx, double* Fy, double* Fz, double const dx, double const dy, double const dz, 
  int const ii, int const jj, int const kk ) {

  double factor1 = -1.*(1.+1.)/dx;
  double factor2 = -1.*(1.+1.)/dy;
  double factor3 = -1.*(1.+1.)/dz;

  Fx[ ii*Nx + jj ] = factor1*( phia[ ( (Nx+ii+1)%Nx )*Ny*Nz + jj*Nz + kk ] 
                             - phia[ ( (Nx+ii-1)%Nx )*Nx + jj*Nz + kk ] );
  Fy[ ii*Nx + jj ] = factor2*( phia[ ii*Ny*Nz + ( (Ny+jj+1)%Ny )*Nz + kk ] 
                             - phia[ ii*Ny*Nz + ( (Ny+jj-1)%Ny )*Nz + kk] );  
  Fz[ ii*Nx + jj ] = factor3*( phia[ ii*Ny*Nz + jj*Nz + (Nz+kk+1)%Nz ] )
                             - phia[ ii*Ny*Nz + jj*Nz + (Nz+kk-1)%Nz ];


  return 0;

}



double _fftw_3d_(int const thread_num, int const dim, int const Nx, int const Ny, int const Nz, 
  double const ori_sigma_a[], double* Fx, double *Fy, double *Fz, 
  double const dx, double const dy, double const dz){

  omp_set_num_threads( thread_num );

  double G_const = 6.67408e-8; // #g^-1 s^-2 cm^3 
  
  double canvas_cal = 20.;
  int total_n = Nx*Ny*Nz;

  fftw_complex *sigma_a, *fftsigma_a;
  fftw_complex *phika, *phia; 
  fftw_plan p, q;

  sigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  fftsigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  phika = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  phia = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);

  /////////// initialization of sigma_a ///////////
  double r2;
  int ii, jj, kk, lb = -1, rb = 1, index;

  # pragma parallel shared( sigma_a ) private( ii, jj, kk, index )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1) {
      for (jj=0; jj<Ny; jj+=1) {
        for (kk=0; kk<Nz; kk+=1){
          index = (Nx/2+ii)*Nz*Ny + (Ny/2+jj)*Nz + kk;
          sigma_a[ index ][0] = ori_sigma_a[ index ];
          sigma_a[ index ][1] = 0; 
        }
      }
    }
  }


  /////////// fft ///////////
  p = fftw_plan_dft_3d(Nx, Ny, Nz, sigma_a, fftsigma_a, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  double kxx, kyy, kzz;

  # pragma parallel shared( phika ) private( ii, jj, kk, index, kxx, kyy, kzz )
  {
    # pragma for  
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){
        for (kk=0; kk<Nz; kk+=1){
          kxx = pow((double)(ii+1)*2.*M_PI/canvas_cal, 2.);
          kyy = pow((double)(jj+1)*2.*M_PI/canvas_cal, 2.);
          kzz = pow((double)(kk+1)*2.*M_PI/canvas_cal, 2.);
          index = ii*Nz*Ny + jj*Nz + kk;
          phika[ index ][0] = -4. * M_PI * G_const * fftsigma_a[ index ][0] / ((kxx+kyy+kzz));
          phika[ index ][1] = -4. * M_PI * G_const * fftsigma_a[ index ][1] / ((kxx+kyy+kzz));
        }
      }
    }
  }

  /////////// inverse fft ///////////
  q = fftw_plan_dft_3d(Nx, Ny, Nz, phika, phia, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(q);
  fftw_destroy_plan(q);
  /////////// normalization ///////////
  # pragma parallel shared( phia ) private( ii, jj, kk, index )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){
        for (kk=0; kk<Nz; kk+=1){
          index = ii*Nz*Ny + jj*Nz + kk;
          phia[ index ][0] /= total_n;
          phia[ index ][1] /= total_n;
        }
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

  double *phi_total = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  ////////// magnitude of potential //////////
  # pragma parallel shared( sigma_a ) private( ii, jj, kk, index )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){ 
        for (kk=0; kk<Nz; kk+=1)
        {
          index = ii*Nz*Ny + jj*Nz + kk;
          phi_total[ index ] = sqrt( pow(phia[ index ][0],2) + pow(phia[ index ][1],2) );        
        }
      }
    } 
  }

  # pragma parallel shared( Fx, Fy, Fz, phi_total ) private( ii, jj, kk )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){
        for (kk=0; kk<Nz; kk+=1){
          _2nd_order_diff_( Nx, Ny, Nz, phi_total, Fx, Fy, Fz, dx, dy, dz, ii, jj, kk );
        }
      }
    }
  }
    
  fftw_cleanup();

  return 0;
}
