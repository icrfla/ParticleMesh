
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

void _2nd_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk ) {

  double factor1 = -1./(2.*grid->dx);
  double factor2 = -1./(2.*grid->dy);
  double factor3 = -1./(2.*grid->dz);
  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  int const Nz = grid->Nz;
  int index = ii + jj*Nx + kk*Nx*Ny;

  grid->Fx[ index ] = factor1*( grid->phi[ ( (Nx+ii+1)%Nx ) + jj*Nx + kk*Nx*Ny ] 
                              - grid->phi[ ( (Nx+ii-1)%Nx ) + jj*Nx + kk*Nx*Ny ] );

  grid->Fy[ index ] = factor2*( grid->phi[ ii + ( (Ny+jj+1)%Ny )*Nx + kk*Nx*Ny ] 
                              - grid->phi[ ii + ( (Ny+jj-1)%Ny )*Nx + kk*Nx*Ny ] );  

  grid->Fz[ index ] = factor3*( grid->phi[ ii + jj*Nx + ((Nz+kk+1)%Nz)*Nx*Ny ] )
                              - grid->phi[ ii + jj*Nx + ((Nz+kk-1)%Nz)*Nx*Ny ];

}

void _4th_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk ) {

  double factor1 = -1./(12.*grid->dx);
  double factor2 = -1./(12.*grid->dy);
  double factor3 = -1./(12.*grid->dz);
  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  int const Nz = grid->Nz;  
  int index = ii + jj*Nx + kk*Nx*Ny;

  grid->Fx[ index ] = factor1*( 
                      +1.*grid->phi[ ( (Nx+ii-2)%Nx ) + jj*Nx + kk*Nx*Ny ] 
                      -8.*grid->phi[ ( (Nx+ii-1)%Nx ) + jj*Nx + kk*Nx*Ny ] 
                      +8.*grid->phi[ ( (Nx+ii+1)%Nx ) + jj*Nx + kk*Nx*Ny ]
                      -1.*grid->phi[ ( (Nx+ii+2)%Nx ) + jj*Nx + kk*Nx*Ny ]
                             );
  
  grid->Fy[ index ] = factor2*( 
                      +1.*grid->phi[ ii + ( (Ny+jj-2)%Ny )*Nx + kk*Nx*Ny ] 
                      -8.*grid->phi[ ii + ( (Ny+jj-1)%Ny )*Nx + kk*Nx*Ny ] 
                      +8.*grid->phi[ ii + ( (Ny+jj+1)%Ny )*Nx + kk*Nx*Ny ] 
                      -1.*grid->phi[ ii + ( (Ny+jj+2)%Ny )*Nx + kk*Nx*Ny ] 
                             );
  
  grid->Fz[ index ] = factor3*( 
                      +1.*grid->phi[ ii + jj*Nx + ((Nz+kk-2)%Nz)*Nx*Ny ] 
                      -8.*grid->phi[ ii + jj*Nx + ((Nz+kk-1)%Nz)*Nx*Ny ] 
                      +8.*grid->phi[ ii + jj*Nx + ((Nz+kk+1)%Nz)*Nx*Ny ] 
                      -1.*grid->phi[ ii + jj*Nx + ((Nz+kk+2)%Nz)*Nx*Ny ] 
                             );

}

void _6th_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk ) {

  double factor1 = -1./(60.*grid->dx);
  double factor2 = -1./(60.*grid->dy);
  double factor3 = -1./(60.*grid->dz);
  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  int const Nz = grid->Nz;  
  int index = ii + jj*Nx + kk*Nx*Ny;

  grid->Fx[ index ] = factor1*( 
                      -1.*grid->phi[ ( (Nx+ii-3)%Nx ) + jj*Nx + kk*Nx*Ny ]
                      +9.*grid->phi[ ( (Nx+ii-2)%Nx ) + jj*Nx + kk*Nx*Ny ] 
                      -45.*grid->phi[ ( (Nx+ii-1)%Nx ) + jj*Nx + kk*Nx*Ny ] 
                      +45.*grid->phi[ ( (Nx+ii+1)%Nx ) + jj*Nx + kk*Nx*Ny ]
                      -9.*grid->phi[ ( (Nx+ii+2)%Nx ) + jj*Nx + kk*Nx*Ny ]
                      +1.*grid->phi[ ( (Nx+ii+3)%Nx ) + jj*Nx + kk*Nx*Ny ]
                             );
  
  grid->Fy[ index ] = factor2*( 
                      -1.*grid->phi[ ii + ( (Ny+jj-3)%Ny )*Nx + kk*Nx*Ny ]
                      +9.*grid->phi[ ii + ( (Ny+jj-2)%Ny )*Nx + kk*Nx*Ny ] 
                      -45.*grid->phi[ ii + ( (Ny+jj-1)%Ny )*Nx + kk*Nx*Ny ] 
                      +45.*grid->phi[ ii + ( (Ny+jj+1)%Ny )*Nx + kk*Nx*Ny ] 
                      -9.*grid->phi[ ii + ( (Ny+jj+2)%Ny )*Nx + kk*Nx*Ny ] 
                      +1.*grid->phi[ ii + ( (Ny+jj+3)%Ny )*Nx + kk*Nx*Ny ] 
                             );
  
  grid->Fz[ index ] = factor3*( 
                      -1.*grid->phi[ ii + jj*Nx + ((Nz+kk-3)%Nz)*Nx*Ny ] 
                      +9.*grid->phi[ ii + jj*Nx + ((Nz+kk-2)%Nz)*Nx*Ny ] 
                      -45.*grid->phi[ ii + jj*Nx + ((Nz+kk-1)%Nz)*Nx*Ny ] 
                      +45.*grid->phi[ ii + jj*Nx + ((Nz+kk+1)%Nz)*Nx*Ny ] 
                      -9.*grid->phi[ ii + jj*Nx + ((Nz+kk+2)%Nz)*Nx*Ny ] 
                      +1.*grid->phi[ ii + jj*Nx + ((Nz+kk+3)%Nz)*Nx*Ny ] 
                             );

}