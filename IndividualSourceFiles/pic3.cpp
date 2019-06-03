#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "Grid.h"
#include "Function.h"

void Weight3d(struct grid3D *grid,struct particle3D *particle,int type){
	//type = 0/1/2  => NGP/CIC/TSC
	for(int i=0;i < particle->number;i++){
		int lx,ly,lz,sx,sy,sz;
		double shift = -grid->L/2;	//make (0,0) to be in the center of grid.
		lx = (particle->x[i]-shift)/grid->dx;
		ly = (particle->y[i]-shift)/grid->dy;
		lz = (particle->z[i]-shift)/grid->dz;

		if(type == 0){
			//NGP
			sx = particle->x[i]-shift-lx * grid->dx - 0.5*grid->dx + 1;
			sy = particle->y[i]-shift-ly * grid->dy - 0.5*grid->dy + 1;
			sz = particle->z[i]-shift-lz * grid->dz - 0.5*grid->dz + 1;
			grid->density[(lx+sx)+(ly+sy) * grid->Nx + (lz+sz) * grid->Nx * grid->Ny] += particle->mass[i];
		}else if(type == 1){
			//CIC
			for(int zz=0;zz<2;zz++){
				for(int j=0;j<4;j++){
					int p = j / 2;
					int q = j % 2;
					double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
					wFactor *= (1-fabs(particle->z[i]-shift-(lz+zz)*grid->dz)/grid->dz);
					grid->density[(lx+p)+(ly+q)*grid->Nx + (lz+zz) * grid->Nx * grid->Ny] += particle->mass[i] * wFactor;
				}
			}
			
		}else if(type == 2){
			//TSC
			//xxx
			//xox
			//xxx
			lx = (particle->x[i]-shift+0.5*grid->dx)/grid->dx;	//Find the nearest point in lattice index.
			ly = (particle->y[i]-shift+0.5*grid->dy)/grid->dy;
			lz = (particle->z[i]-shift+0.5*grid->dz)/grid->dz;
			double weightX[3];	//Weight factor in x direction for 3 affected points
			double weightY[3];
			double weightZ[3];
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
			for(int zz=-1;zz<2;zz++){
				double ddz = fabs(particle->z[i]-shift-(lz+zz)*grid->dz);
				if(ddz <= grid->dz/2){
					weightZ[zz+1] = 3.0/4 - pow(ddz / grid->dz,2);
				}else if(ddz <= grid->dz/2*3.0){
					weightZ[zz+1] = 0.5*pow(1.5-ddz / grid->dz,2);
				}else{
					weightZ[zz+1]=0.0;
				}
			}
			//Weight mass into density
			int indx,indy,indz;
			for(int xx=-1;xx<2;xx++){
				for(int yy=-1;yy<2;yy++){
					for(int zz=-1;zz<2;zz++){
						//Account for periodic boundary
						indx = ((lx+xx)+grid->Nx)%grid->Nx;
						indy = ((ly+yy)+grid->Ny)%grid->Ny;
						indz = ((lz+zz)+grid->Nz)%grid->Nz;
						int index = indx + indy*grid->Nx + indz*grid->Nx*grid->Ny;
						grid->density[index]+=weightX[xx+1]*weightY[yy+1]*weightZ[zz+1]*particle->mass[i];
					}
					
				}
			}
		}
		
	}
}
void WeightForce3d(struct grid3D *grid,struct particle3D *particle,int type){
	//type = 0/1/2  => NGP/CIC/TSC
	for(int i=0;i < particle->number;i++){
		int lx,ly,lz,sx,sy,sz;
		double shift = -grid->L/2;	//make (0,0) to be in the center of grid.
		lx = (particle->x[i]-shift)/grid->dx;
		ly = (particle->y[i]-shift)/grid->dy;
		lz = (particle->z[i]-shift)/grid->dz;
		particle->Fx[i]=0.0;
		particle->Fy[i]=0.0;
		particle->Fz[i]=0.0;

		if(type == 0){
			sx = particle->x[i]-shift-lx * grid->dx - 0.5*grid->dx + 1;
			sy = particle->y[i]-shift-ly * grid->dy - 0.5*grid->dy + 1;
			sz = particle->z[i]-shift-lz * grid->dz - 0.5*grid->dz + 1;
			int pos = (lx+sx)+(ly+sy) * grid->Nx + (lz+sz) * grid->Nx * grid->Ny;
			particle->Fx[i]=grid->Fx[pos];
			particle->Fy[i]=grid->Fy[pos];
			particle->Fz[i]=grid->Fz[pos];
		}else if(type == 1){
			for(int zz=0;zz<2;zz++){
				for(int j=0;j<4;j++){
					int p = j / 2;
					int q = j % 2;
					double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
					wFactor *= (1-fabs(particle->z[i]-shift-(lz+zz)*grid->dz)/grid->dz);
					int pos = (lx+p)+(ly+q)*grid->Nx + (lz+zz) * grid->Nx * grid->Ny;
					particle->Fx[i] += grid->Fx[pos] * wFactor;
					particle->Fy[i] += grid->Fy[pos] * wFactor;
					particle->Fz[i] += grid->Fz[pos] * wFactor;
				}
			}
			
		}else if(type == 2){
			//TSC
			//xxx
			//xox
			//xxx
			lx = (particle->x[i]-shift+0.5*grid->dx)/grid->dx;	//Find the nearest point in lattice index.
			ly = (particle->y[i]-shift+0.5*grid->dy)/grid->dy;
			lz = (particle->z[i]-shift+0.5*grid->dz)/grid->dz;
			double weightX[3];	//Weight factor in x direction for 3 affected points
			double weightY[3];
			double weightZ[3];
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
			for(int zz=-1;zz<2;zz++){
				double ddz = fabs(particle->z[i]-shift-(lz+zz)*grid->dz);
				if(ddz <= grid->dz/2){
					weightZ[zz+1] = 3.0/4 - pow(ddz / grid->dz,2);
				}else if(ddz <= grid->dz/2*3.0){
					weightZ[zz+1] = 0.5*pow(1.5-ddz / grid->dz,2);
				}else{
					weightZ[zz+1]=0.0;
				}
			}
			//Weight mass into density
			int indx,indy,indz;
			for(int xx=-1;xx<2;xx++){
				for(int yy=-1;yy<2;yy++){
					for(int zz=-1;zz<2;zz++){
						//Account for periodic boundary
						indx = ((lx+xx)+grid->Nx)%grid->Nx;
						indy = ((ly+yy)+grid->Ny)%grid->Ny;
						indz = ((lz+zz)+grid->Nz)%grid->Nz;
						int index = indx + indy*grid->Nx + indz*grid->Nx*grid->Ny;
						double weight = weightX[xx+1]*weightY[yy+1]*weightZ[zz+1];
						particle->Fx[i] += grid->Fx[index]*weight;
						particle->Fy[i] += grid->Fy[index]*weight;
						particle->Fz[i] += grid->Fz[index]*weight;
					}
					
				}
			}
		}

		
	}
}


