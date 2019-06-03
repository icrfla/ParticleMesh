// g++-9 parallel_poisson_solver_fft_2d.cpp -L/usr/local/Cellar/fftw/3.3.8_1/lib -I/usr/local/Cellar/fftw/3.3.8_1/include -lfftw3 -fopenmp

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

double parallel_force_2d(int const thread_num, int const dim, int const Nx, int const Ny, 
  double const ori_sigma_a[], double* Fx, double *Fy, 
  double const dx, double const dy);

double _2nd_order_diff_2d(int const Nx, int const Ny, double const phi[], 
  double* Fx, double* Fy, double const dx, double const dy, 
  int const ii, int const jj );

double _4th_order_diff_2d(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy, int const ii, int const jj );

double _6th_order_diff_2d(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy, int const ii, int const jj );

int main( int argc, char *argv[] ){

  int canvas_size = 20, Nx = 20, Ny = 20;
  int dim = 2;
  int thread_num = 2;
  double dx = 1., dy = 1.;


  double * sigma_a = (double *) malloc(Nx*Ny * sizeof(double));
  /////////// initialization of sigma_a ///////////
  double r2;
  int ii, jj, lb = -1, rb = 1;
  for (ii=lb; ii<=rb; ii+=1) {
    for (jj=lb; jj<=rb; jj+=1) {
      r2 = (double) (pow(ii,2)+pow(jj,2));
      printf("ii:%d jj:%d\n", ii, jj);
      sigma_a[ (Nx/2+ii)*Ny + (Ny/2+jj) ] = 1. * sqrt(1. - r2/1600.);
    }
  }

  double *Fx = (double *) malloc(Nx*Ny * sizeof(double));
  double *Fy = (double *) malloc(Nx*Ny * sizeof(double));

  parallel_poisson_solver_2d(thread_num, dim, Nx, Ny, sigma_a, Fx, Fy, dx, dy);

  return 0;
}

double parallel_force_2d(int const thread_num, int const dim, int const Nx, int const Ny, 
  double const ori_sigma_a[], double* Fx, double *Fy, double const dx, double const dy) {

  omp_set_num_threads( thread_num );

  double G_const = 6.67408e-8; // #g^-1 s^-2 cm^3 
  int total_n = Nx*Ny;
  double dNx = (double) (Nx), dNy = (double) (Ny); // default Nx=Ny
  int ii, jj, index;

  fftw_complex *sigma_a, *fftsigma_a;
  fftw_complex *phika, *phia; 
  fftw_plan p, q;

  sigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  fftsigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  phika = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  phia = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);


  # pragma parallel shared( sigma_a ) private( ii, jj )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1) {
      for (jj=0; jj<Ny; jj+=1) {
        index = ii + jj*Nx;
        sigma_a[ index ][0] = ori_sigma_a[ index ]; 
        sigma_a[ index ][1] = 0.;
      }
    }
  }

  /////////// fft ///////////
  p = fftw_plan_dft_2d(Nx, Ny, sigma_a, fftsigma_a, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  double kxx, kyy;
  # pragma parallel shared( phika ) private( ii, jj, kxx, kyy )
  { 
    # pragma for
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){
        index = ii + jj*Nx;
        kxx = pow((double)(ii+1)*2.*M_PI/dNx, 2.);
        kyy = pow((double)(jj+1)*2.*M_PI/dNy, 2.);
        phika[index][0] = -4. * M_PI * G_const * fftsigma_a[index][0] / ((kxx+kyy));
        phika[index][1] = -4. * M_PI * G_const * fftsigma_a[index][1] / ((kxx+kyy));
      }
    }
  }

  /////////// inverse fft ///////////
  q = fftw_plan_dft_2d(Nx, Ny, phika, phia, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(q);
  fftw_destroy_plan(q);
  /////////// normalization ///////////
  # pragma parallel shared( phia ) private( ii, jj, total_n )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){
        index = ii + jj*Nx;
        phia[index][0] /= total_n;
        phia[index][1] /= total_n;
      }
    }
  }

  double *phi_total = (double *) malloc(Nx*Ny * sizeof(double));
  ////////// magnitude of potential //////////
  # pragma parallel shared( phi_total ) private( ii, jj )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){
        index = ii + jj*Nx;
        phi_total[index] = sqrt( pow(phia[index][0],2) + pow(phia[index][1],2) );
      }
    }  
  }

  # pragma parallel shared( Fx, Fy ) private( ii, jj )
  {
    # pragma for
    for (ii=0; ii<Nx; ii+=1){
      for (jj=0; jj<Ny; jj+=1){     
        _2nd_order_diff_2d(Nx, Ny, phi_total, Fx, Fy, dx, dy, ii, jj);
      }
    }
  }

  for (ii=0; ii<Nx; ii+=1) {
    for (jj=0; jj<Ny; jj+=1){
      index = ii + jj*Nx;
      printf("recover: %3d %3d %+.5e %+.5e I vs. %+.5e %+.5e I \n",
          ii, jj, sigma_a[index][0], sigma_a[index][1], 
          phia[index][0], phia[index][1]);    
    }
  }

  for (ii=0; ii<Nx; ii+=1) {
    for (jj=0; jj<Ny; jj+=1){
      index = ii + jj*Nx;
      printf("ii:%3d, jj:%3d, Fx:%.5e Fy:%.5e\n",
          ii, jj, Fx[index], Fy[index]);    
    }
  }


  fftw_cleanup();

  return 0;
}

double _2nd_order_diff_2d(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy, int const ii, int const jj ) {

  double factor1 = -1./(2.*dx);
  double factor2 = -1./(2.*dy);
  int index = ii + jj*Nx;

  Fx[ index ] = factor1*( phia[ ( (Nx+ii+1)%Nx ) + jj*Nx ] - phia[ ( (Nx+ii-1)%Nx ) + jj*Nx ] );
  Fy[ index ] = factor2*( phia[ ii + ((Ny+jj+1)%Ny)*Nx ] - phia[ ii + ((Ny+jj-1)%Ny)*Nx ] );  

  return 0;

}

double _4th_order_diff_2d(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy, int const ii, int const jj ) {

  double factor1 = -1./(12.*dx);
  double factor2 = -1./(12.*dy);
  int index = ii + jj*Nx;

  Fx[ index ] = factor1*( 
    +1.*phia[ ( (Nx+ii-2)%Nx ) + jj*Nx ] 
    -8.*phia[ ( (Nx+ii-1)%Nx ) + jj*Nx ] 
    +8.*phia[ ( (Nx+ii+1)%Nx ) + jj*Nx ]
    -1.*phia[ ( (Nx+ii+2)%Nx ) + jj*Nx ]
    );
  

  Fy[ index ] = factor2*( 
    +1.*phia[ ii + ((Ny+jj-2)%Ny)*Nx ] 
    -8.*phia[ ii + ((Ny+jj-1)%Ny)*Nx ] 
    +8.*phia[ ii + ((Ny+jj+1)%Ny)*Nx ] 
    -1.*phia[ ii + ((Ny+jj+2)%Ny)*Nx ] 
    );


  return 0;

}


double _6th_order_diff_2d(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy, int const ii, int const jj ) {

  double factor1 = -1./(60.*dx);
  double factor2 = -1./(60.*dy);
  int index = ii + jj*Nx;

  Fx[ index ] = factor1*( 
    -1.* phia[ ( (Nx+ii-3)%Nx )*Nx + jj ]
    +9.* phia[ ( (Nx+ii-2)%Nx )*Nx + jj ] 
    -45.*phia[ ( (Nx+ii-1)%Nx )*Nx + jj ] 
    +45.*phia[ ( (Nx+ii+1)%Nx )*Nx + jj ]
    -9.* phia[ ( (Nx+ii+2)%Nx )*Nx + jj ]
    +1.* phia[ ( (Nx+ii+3)%Nx )*Nx + jj ]
    );
  

  Fy[ index ] = factor2*( 
    -1.* phia[ ii + ((Ny+jj-3)%Ny)*Nx ] 
    +9.* phia[ ii + ((Ny+jj-2)%Ny)*Nx ] 
    -45.*phia[ ii + ((Ny+jj-1)%Ny)*Nx ] 
    +45.*phia[ ii + ((Ny+jj+1)%Ny)*Nx ] 
    -9.* phia[ ii + ((Ny+jj+2)%Ny)*Nx ] 
    +1.* phia[ ii + ((Ny+jj+3)%Ny)*Nx ] 
    );


  return 0;

}
