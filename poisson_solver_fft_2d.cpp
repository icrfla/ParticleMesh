// g++-9 poisson_solver_fft_2d.cpp -L/usr/local/Cellar/fftw/3.3.8_1/lib -I/usr/local/Cellar/fftw/3.3.8_1/include -lfftw3

#include <cstdio>
#include <cstdlib>
// #include <omp.h>
#include <math.h>
#include <string>
#include <iostream>
#include <time.h>
#include <fftw3.h>
#include <complex.h>
using namespace std;

double _fftw_2d_(int const dim, int const Nx, int const Ny, 
  double const ori_sigma_a[], double* Fx, double *Fy, 
  double const dx, double const dy);

double _2nd_order_diff_(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy, int const ii, int const jj );

double _4th_order_diff_(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy );

double _6th_order_diff_(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy );

int main( int argc, char *argv[] ){

  int canvas_size = 20, Nx = 20, Ny = 20;
  int dim = 2;
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

  _fftw_2d_(dim, Nx, Ny, sigma_a, Fx, Fy, dx, dy);

  return 0;
}

double _fftw_2d_(int const dim, int const Nx, int const Ny, 
  double const ori_sigma_a[], double* Fx, double *Fy, 
  double const dx, double const dy) {

  double G_const = 6.67408e-8; // #g^-1 s^-2 cm^3 
  int total_n = Nx*Ny;
  double dNx = (double) (Nx), dNy = (double) (Ny); // default Nx=Ny
  int ii, jj;

  fftw_complex *sigma_a, *fftsigma_a;
  fftw_complex *phika, *phia; 
  fftw_plan p, q;

  sigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  fftsigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  phika = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  phia = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

  for (ii=0; ii<Nx; ii+=1) {
    for (jj=0; jj<Ny; jj+=1) {
      sigma_a[ ii*Ny + jj ][0] = ori_sigma_a[ ii*Ny + jj ]; 
      sigma_a[ ii*Ny + jj ][1] = 0.;
    }
  }


  /////////// fft ///////////
  p = fftw_plan_dft_2d(Nx, Ny, sigma_a, fftsigma_a, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  double kxx, kyy;
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      kxx = pow((double)(ii+1)*2.*M_PI/dNx, 2.);
      kyy = pow((double)(jj+1)*2.*M_PI/dNy, 2.);
      phika[ii*Nx + jj][0] = -4. * M_PI * G_const * fftsigma_a[ii*Nx + jj][0] / ((kxx+kyy));
      phika[ii*Nx + jj][1] = -4. * M_PI * G_const * fftsigma_a[ii*Nx + jj][1] / ((kxx+kyy));
    }
  }
  /////////// inverse fft ///////////
  q = fftw_plan_dft_2d(Nx, Ny, phika, phia, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(q);
  fftw_destroy_plan(q);
  /////////// normalization ///////////
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      phia[ii*Nx + jj][0] /= total_n;
      phia[ii*Nx + jj][1] /= total_n;
    }
  }

  double *phi_total = (double *) malloc(Nx*Ny * sizeof(double));
  ////////// magnitude of potential //////////
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){ 
      phi_total[ii*Nx + jj] = sqrt( pow(phia[ii*Nx + jj][0],2) + pow(phia[ii*Nx + jj][1],2) );
    }
  }  

  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){     
      _2nd_order_diff_(Nx, Ny, phi_total, Fx, Fy, dx, dy, ii, jj);
    }
  }

  for (ii=0; ii<Nx; ii+=1) {
    for (jj=0; jj<Ny; jj+=1){
      printf("recover: %3d %3d %+.5e %+.5e I vs. %+.5e %+.5e I \n",
          ii, jj, sigma_a[ii*Nx + jj][0], sigma_a[ii*Nx + jj][1], 
          phia[ii*Nx + jj][0], phia[ii*Nx + jj][1]);    
    }
  }

  for (ii=0; ii<Nx; ii+=1) {
    for (jj=0; jj<Ny; jj+=1){
      printf("ii:%3d, jj:%3d, Fx:%.5e Fy:%.5e\n",
          ii, jj, Fx[ii*Nx + jj], Fy[ii*Nx + jj]);    
    }
  }


  fftw_cleanup();

  return 0;
}

double _2nd_order_diff_(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy, int const ii, int const jj ) {

  double factor1 = -1.*(1.+1.)/dx;
  double factor2 = -1.*(1.+1.)/dy;

  Fx[ ii*Nx + jj ] = factor1*( phia[ ( (Nx+ii+1)%Nx )*Ny + jj ] - phia[ ( (Nx+ii-1)%Nx )*Ny + jj ] );
  Fy[ ii*Nx + jj ] = factor2*( phia[ ii*Ny + (Ny+jj+1)%Ny ] - phia[ ii*Ny + (Ny+jj-1)%Ny ] );  

  return 0;

}

double _4th_order_diff_(int const Nx, int const Ny, double const phia[], 
  double* Fx, double* Fy, double const dx, double const dy ) {

  double factor1 = -1.*(1.+1.)/dx;
  double factor2 = -1.*(1.+1.)/dy;

  int ii, jj;
  for (ii=0; ii<Nx; ii+=1){
    for (jj=0; jj<Ny; jj+=1){
      Fx[ ii*Nx + jj ] = factor1*( phia[ ( (Nx+ii+1)%Nx )*Nx + jj ] - phia[ ( (Nx+ii-1)%Nx )*Nx + jj ] );
      Fy[ ii*Nx + jj ] = factor2*( phia[ ii*Nx + (Ny+jj+1)%Ny ] - phia[ ii*Nx + (Ny+jj-1)%Ny ] );  
    } // for jj
  } // for ii

  return 0;

}
