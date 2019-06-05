
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>
#include <fftw3.h>
#include <complex.h>

#include "Grid.h"
#include "Function.h"

using namespace std;

void Weight2d(struct grid2D *grid, struct particle2D *particle,int type){
	//type = 0/1/2  => NGP/CIC/TSC
	for(int i=0;i < particle->number;i++){
		int lx,ly,sx,sy;
		double shift = -grid->L/2;	//make (0,0) to be in the center of grid.
		lx = (particle->x[i]-shift)/grid->dx;
		ly = (particle->y[i]-shift)/grid->dy;

		if(type == 0){
			//NGP
			sx = particle->x[i]-shift-lx * grid->dx - 0.5*grid->dx + 1;
			sy = particle->y[i]-shift-ly * grid->dy - 0.5*grid->dy + 1;
			grid->density[(lx+sx)+(ly+sy) * grid->Nx] += particle->mass[i];
		}else if(type == 1){
			//CIC
			//xx
			//ox
			for(int j=0;j<4;j++){
				int p = j / 2;
				int q = j % 2;
				double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
				grid->density[(lx+p)+(ly+q)*grid->Nx] += particle->mass[i] * wFactor;
			}
		}else if(type == 2){
			//TSC
			//xxx
			//xox
			//xxx
			lx = (particle->x[i]-shift+0.5*grid->dx)/grid->dx;	//Find the nearest point in lattice index.
			ly = (particle->y[i]-shift+0.5*grid->dy)/grid->dy;
			double weightX[3];	//Weight factor in x direction for 3 affected points
			double weightY[3];
			//Construct weighting factor
			for(int xx=-1;xx<2;xx++){
				double ddx = fabs(particle->x[i]-shift-(lx+xx)*grid->dx);
				if(ddx <= grid->dx/2){
					weightX[xx+1] = 3.0/4 - pow(ddx / grid->dx,2);
				}else if(ddx<= grid->dx/2*3.0){
					weightX[xx+1] = 0.5*pow(1.5-ddx / grid->dx,2);
				}else{
					printf("Should not be here");
					weightX[xx+1]=0.0;
				}
				
			}
			for(int yy=-1;yy<2;yy++){
				double ddy = fabs(particle->y[i]-shift-(ly+yy)*grid->dy);
				if(ddy <= grid->dy/2){
					weightY[yy+1] = 3.0/4 - pow(ddy / grid->dy,2);
				}else if(ddy <= grid->dy/2*3.0){
					weightY[yy+1] = 0.5*pow(1.5-ddy / grid->dy,2);
				}else{
					weightY[yy+1]=0.0;
				}
			}
			//Weight mass into density
			int indx,indy;
			for(int xx=-1;xx<2;xx++){
				for(int yy=-1;yy<2;yy++){
					//Account for periodic boundary
					indx = ((lx+xx)+grid->Nx)%grid->Nx;
					indy = ((ly+yy)+grid->Ny)%grid->Ny;
					grid->density[indx+indy*grid->Nx]+=weightX[xx+1]*weightY[yy+1]*particle->mass[i];
				}
			}

		}
		
	}
}
void WeightForce2d(struct grid2D *grid,struct particle2D *particle,int type){
	//type = 0/1/2  => NGP/CIC/TSC
	for(int i=0;i < particle->number;i++){
		int lx,ly,sx,sy;
		double shift = -grid->L/2;	//make (0,0) to be in the center of grid.
		lx = (particle->x[i]-shift)/grid->dx;
		ly = (particle->y[i]-shift)/grid->dy;
		particle->Fx[i]=0.0;
		particle->Fy[i]=0.0;

		if(type == 0){
			sx = particle->x[i]-shift-lx * grid->dx - 0.5*grid->dx + 1;
			sy = particle->y[i]-shift-ly * grid->dy - 0.5*grid->dy + 1;
			particle->Fx[i]=grid->Fx[(lx+sx)+(ly+sy)*grid->Nx];
			particle->Fy[i]=grid->Fy[(lx+sx)+(ly+sy)*grid->Nx];
		}else if(type == 1){
			for(int j=0;j<4;j++){
				int p = j / 2;
				int q = j % 2;
				double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
				particle->Fx[i] += grid->Fx[(lx+p)+(ly+q)*grid->Nx] * wFactor;
				particle->Fy[i] += grid->Fy[(lx+p)+(ly+q)*grid->Nx] * wFactor;
			}
		}else if(type == 2){
			//TSC
			//xxx
			//xox
			//xxx
			lx = (particle->x[i]-shift+0.5*grid->dx)/grid->dx;	//Find the nearest point in lattice index.
			ly = (particle->y[i]-shift+0.5*grid->dy)/grid->dy;
			double weightX[3];	//Weight factor in x direction for 3 affected points
			double weightY[3];
			//Construct weighting factor
			for(int xx=-1;xx<2;xx++){
				double ddx = fabs(particle->x[i]-shift-(lx+xx)*grid->dx);
				if(ddx <= grid->dx/2){
					weightX[xx+1] = 3.0/4 - pow(ddx / grid->dx,2);
				}else if(ddx<= grid->dx/2*3.0){
					weightX[xx+1] = 0.5*pow(1.5-ddx / grid->dx,2);
				}else{
					printf("Should not be here");
					weightX[xx+1]=0.0;
				}
				
			}
			for(int yy=-1;yy<2;yy++){
				double ddy = fabs(particle->y[i]-shift-(ly+yy)*grid->dy);
				if(ddy <= grid->dy/2){
					weightY[yy+1] = 3.0/4 - pow(ddy / grid->dy,2);
				}else if(ddy <= grid->dy/2*3.0){
					weightY[yy+1] = 0.5*pow(1.5-ddy / grid->dy,2);
				}else{
					weightY[yy+1]=0.0;
				}
			}
			//Weight mass into density
			int indx,indy;
			for(int xx=-1;xx<2;xx++){
				for(int yy=-1;yy<2;yy++){
					//Account for periodic boundary
					indx = ((lx+xx)+grid->Nx)%grid->Nx;
					indy = ((ly+yy)+grid->Ny)%grid->Ny;
					particle->Fx[i] += grid->Fx[indx+indy*grid->Nx]*weightX[xx+1]*weightY[yy+1];
					particle->Fy[i] += grid->Fy[indx+indy*grid->Nx]*weightX[xx+1]*weightY[yy+1];
				}
			}
		}

		
	}
}

