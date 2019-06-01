#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

struct particle3D{
	int number;
	double *mass;
	double *x;
	double *y;
	double *z;
	double *Fx;
	double *Fy;
	double *Fz;
	double *vx;
	double *vy;
	double *vz;
};
struct grid3D{
	double L;
	int Nx;
	int Ny;
	int Nz;
	int N;
	double dx;
	double dy;
	double dz;
	double *density;
	double *Fx;
	double *Fy;
	double *Fz;
};

void Weight(struct grid3D *grid,struct particle3D *particle,int type){
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
void WeightForce(struct grid3D *grid,struct particle3D *particle,int type){
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

int main(){
	int weightFunction = 1;  //0/1/2 : NGP/CIC/TSC

	//================Random number generator.
	//To use : d=gsl_rng_uniform(rng);
	gsl_rng *rng;
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,123456);//The seed is 123456.

	//================Grid 
	struct grid3D grid;
	grid.L = 10.0;		//Length of box (-L/2 ~ L/2)
	grid.Nx = 11;		//Number of grid in x direction. (Should be an odd number)
	grid.Ny = grid.Nx;
	grid.Nz = grid.Nx;
	grid.N = grid.Nx * grid.Ny * grid.Nz;
	grid.dx = grid.L / (grid.Nx-1);	//-1 is because boundary (make box closed)
	grid.dy = grid.L / (grid.Ny-1);
	grid.dz = grid.L / (grid.Nz-1);
	grid.density = (double*) malloc(grid.N*sizeof(double));	//Density on the grid
	grid.Fx = (double*) malloc(grid.N*sizeof(double));			//Force on grid comes from potential calculation
	grid.Fy = (double*) malloc(grid.N*sizeof(double));
	grid.Fz = (double*) malloc(grid.N*sizeof(double));

	

	//================Particles
	int NParticle=10;//Number of particles used in simulation
	double massParticle=1.0;
		
	struct particle3D myParticle;
	myParticle.number = NParticle;
	myParticle.mass = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.x = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.y = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.z = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.Fx = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.Fy = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.Fz = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.vx = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.vy = (double*)malloc(myParticle.number*sizeof(double));
	myParticle.vz = (double*)malloc(myParticle.number*sizeof(double));
	

	//Initialize density grid
	for(int i=0;i<grid.N;i++){
		grid.density[i]=0.0;		
	}

	//Initialize Initial Position of Particle.
	for (int i = 0; i < myParticle.number; ++i){
		myParticle.mass[i]=1.0;
		myParticle.x[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
		myParticle.y[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
		myParticle.z[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
		printf("At (%f,%f) \n",myParticle.x[i],myParticle.y[i]);
	}

	//Initialize Force field (for test use)
	for(int i=0;i<grid.N;i++){
		grid.Fx[i]=0.0;
		grid.Fy[i]=0.0;
		grid.Fz[i]=0.0;
		
	}
	grid.Fx[3+6*grid.Nx] = -1.0;
	grid.Fy[3+6*grid.Nx] = 2.0;

	Weight(&grid,&myParticle,weightFunction);

	for(int i=0;i < grid.Nx * grid.Ny;i++){
		if(i % grid.Nx == 0){
			printf("\n");
		}
		printf("%.2f\t",grid.density[i]);
	}
	printf("\n");

	//Remap the force to particle.
	WeightForce(&grid,&myParticle,weightFunction);

	
	//Print out the force on a particle.
	for(int i=0;i<myParticle.number;i++){
		printf("(Fx,Fy,Fz)=(%.2f,%.2f)\n",myParticle.Fx[i],myParticle.Fy[i],myParticle.Fz[i]);
	}



	return 0;
}

