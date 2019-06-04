
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

void _2nd_order_diff_2d(struct grid2D *grid, int const ii, int const jj ) {

  double factor1 = -1./(2.*grid->dx);
  double factor2 = -1./(2.*grid->dy);
  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  int index = ii + jj*Nx;

  grid->Fx[ index ] = factor1*( grid->phi[ ( (Nx+ii+1)%Nx ) + jj*Nx ] - grid->phi[ ( (Nx+ii-1)%Nx ) + jj*Nx ] );
  grid->Fy[ index ] = factor2*( grid->phi[ ii + ((Ny+jj+1)%Ny)*Nx ] - grid->phi[ ii + ((Ny+jj-1)%Ny)*Nx ] );  

}

void _4th_order_diff_2d(struct grid2D *grid, int const ii, int const jj ) {

  double factor1 = -1./(12.*grid->dx);
  double factor2 = -1./(12.*grid->dy);
  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  int index = ii + jj*Nx;

  grid->Fx[ index ] = factor1*( 
    +1.*grid->phi[ ( (Nx+ii-2)%Nx ) + jj*Nx ] 
    -8.*grid->phi[ ( (Nx+ii-1)%Nx ) + jj*Nx ] 
    +8.*grid->phi[ ( (Nx+ii+1)%Nx ) + jj*Nx ]
    -1.*grid->phi[ ( (Nx+ii+2)%Nx ) + jj*Nx ]
    );
  

  grid->Fy[ index ] = factor2*( 
    +1.*grid->phi[ ii + ((Ny+jj-2)%Ny)*Nx ] 
    -8.*grid->phi[ ii + ((Ny+jj-1)%Ny)*Nx ] 
    +8.*grid->phi[ ii + ((Ny+jj+1)%Ny)*Nx ] 
    -1.*grid->phi[ ii + ((Ny+jj+2)%Ny)*Nx ] 
    );

}


void _6th_order_diff_2d(struct grid2D *grid, int const ii, int const jj ) {

  double factor1 = -1./(60.*grid->dx);
  double factor2 = -1./(60.*grid->dy);
  int const Nx = grid->Nx;
  int const Ny = grid->Ny;  
  int index = ii + jj*Nx;

  grid->Fx[ index ] = factor1*( 
    -1.* grid->phi[ ( (Nx+ii-3)%Nx )*Nx + jj ]
    +9.* grid->phi[ ( (Nx+ii-2)%Nx )*Nx + jj ] 
    -45.*grid->phi[ ( (Nx+ii-1)%Nx )*Nx + jj ] 
    +45.*grid->phi[ ( (Nx+ii+1)%Nx )*Nx + jj ]
    -9.* grid->phi[ ( (Nx+ii+2)%Nx )*Nx + jj ]
    +1.* grid->phi[ ( (Nx+ii+3)%Nx )*Nx + jj ]
    );
  

  grid->Fy[ index ] = factor2*( 
    -1.* grid->phi[ ii + ((Ny+jj-3)%Ny)*Nx ] 
    +9.* grid->phi[ ii + ((Ny+jj-2)%Ny)*Nx ] 
    -45.*grid->phi[ ii + ((Ny+jj-1)%Ny)*Nx ] 
    +45.*grid->phi[ ii + ((Ny+jj+1)%Ny)*Nx ] 
    -9.* grid->phi[ ii + ((Ny+jj+2)%Ny)*Nx ] 
    +1.* grid->phi[ ii + ((Ny+jj+3)%Ny)*Nx ] 
    );

}