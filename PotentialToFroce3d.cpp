
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

double _2nd_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk ) {

  double factor1 = -1./(2.*dx);
  double factor2 = -1./(2.*dy);
  double factor3 = -1./(2.*dz);

  grid->Fx[ ii*Nx + jj ] = factor1*( grid->phi[ ( (Nx+ii+1)%Nx )*Ny*Nz + jj*Nz + kk ] 
                             - grid->phi[ ( (Nx+ii-1)%Nx )*Ny*Nz + jj*Nz + kk ] );
  grid->Fy[ ii*Nx + jj ] = factor2*( grid->phi[ ii*Ny*Nz + ( (Ny+jj+1)%Ny )*Nz + kk ] 
                             - grid->phi[ ii*Ny*Nz + ( (Ny+jj-1)%Ny )*Nz + kk] );  
  grid->Fz[ ii*Nx + jj ] = factor3*( grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk+1)%Nz ] )
                             - grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk-1)%Nz ];


  return 0;

}

double _4th_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk ) {

  double factor1 = -1.*(1.+1.)/dx;
  double factor2 = -1.*(1.+1.)/dy;
  double factor3 = -1.*(1.+1.)/dz;

  grid->Fx[ ii*Nx + jj ] = factor1*( 
                      +1.*grid->phi[ ( (Nx+ii-2)%Nx )*Ny*Nz + jj*Nz + kk ] 
                      -8.*grid->phi[ ( (Nx+ii-1)%Nx )*Ny*Nz + jj*Nz + kk ] 
                      +8.*grid->phi[ ( (Nx+ii+1)%Nx )*Ny*Nz + jj*Nz + kk ]
                      -1.*grid->phi[ ( (Nx+ii+2)%Nx )*Ny*Nz + jj*Nz + kk ]
                             );
  
  grid->Fy[ ii*Nx + jj ] = factor2*( 
                      +1.*grid->phi[ ii*Ny*Nz + ( (Ny+jj-2)%Ny )*Nz + kk ] 
                      -8.*grid->phi[ ii*Ny*Nz + ( (Ny+jj-1)%Ny )*Nz + kk ] 
                      +8.*grid->phi[ ii*Ny*Nz + ( (Ny+jj+1)%Ny )*Nz + kk ] 
                      -1.*grid->phi[ ii*Ny*Nz + ( (Ny+jj+2)%Ny )*Nz + kk ] 
                             );
  
  grid->Fz[ ii*Nx + jj ] = factor3*( 
                      +1.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk-2)%Nz ] 
                      -8.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk-1)%Nz ] 
                      +8.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk+1)%Nz ] 
                      -1.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk+2)%Nz ] 
                             );

  return 0;

}

double _6th_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk ) {

  double factor1 = -1.*(1.+1.)/dx;
  double factor2 = -1.*(1.+1.)/dy;
  double factor3 = -1.*(1.+1.)/dz;

  grid->Fx[ ii*Nx + jj ] = factor1*( 
                      -1.*grid->phi[ ( (Nx+ii-3)%Nx )*Ny*Nz + jj*Nz + kk ]
                      +9.*grid->phi[ ( (Nx+ii-2)%Nx )*Ny*Nz + jj*Nz + kk ] 
                      -45.*grid->phi[ ( (Nx+ii-1)%Nx )*Ny*Nz + jj*Nz + kk ] 
                      +45.*grid->phi[ ( (Nx+ii+1)%Nx )*Ny*Nz + jj*Nz + kk ]
                      -9.*grid->phi[ ( (Nx+ii+2)%Nx )*Ny*Nz + jj*Nz + kk ]
                      +1.*grid->phi[ ( (Nx+ii+3)%Nx )*Ny*Nz + jj*Nz + kk ]
                             );
  
  grid->Fy[ ii*Nx + jj ] = factor2*( 
                      -1.*grid->phi[ ii*Ny*Nz + ( (Ny+jj-3)%Ny )*Nz + kk ]
                      +9.*grid->phi[ ii*Ny*Nz + ( (Ny+jj-2)%Ny )*Nz + kk ] 
                      -45.*grid->phi[ ii*Ny*Nz + ( (Ny+jj-1)%Ny )*Nz + kk ] 
                      +45.*grid->phi[ ii*Ny*Nz + ( (Ny+jj+1)%Ny )*Nz + kk ] 
                      -9.*grid->phi[ ii*Ny*Nz + ( (Ny+jj+2)%Ny )*Nz + kk ] 
                      +1.*grid->phi[ ii*Ny*Nz + ( (Ny+jj+3)%Ny )*Nz + kk ] 
                             );
  
  grid->Fz[ ii*Nx + jj ] = factor3*( 
                      -1.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk-3)%Nz ] 
                      +9.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk-2)%Nz ] 
                      -45.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk-1)%Nz ] 
                      +45.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk+1)%Nz ] 
                      -9.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk+2)%Nz ] 
                      +1.*grid->phi[ ii*Ny*Nz + jj*Nz + (Nz+kk+3)%Nz ] 
                             );

  return 0;

}