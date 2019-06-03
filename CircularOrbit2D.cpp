//To compile should add "-lfftw3 -lgsl"
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
using namespace std;

#include "Grid.h"
#include "Function.h"

// struct particle2D{
// 	int number;
// 	double *mass;
// 	double *x;
// 	double *y;
// 	double *Fx;
// 	double *Fy;
// 	double *vx;
// 	double *vy;
// };
// struct grid2D{
// 	double L;
// 	int Nx;
// 	int Ny;
// 	int N;
// 	double dx;
// 	double dy;
// 	double *density;
// 	double *phi;
// 	double *Fx;
// 	double *Fy;
// };
void Weight(struct grid2D *grid, struct particle2D *particle,int type){
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
void WeightForce(struct grid2D *grid,struct particle2D *particle,int type){
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

// void _fftw_2d_(int const dim, struct grid2D *grid);

// void _2nd_order_diff_(struct grid2D *grid, int const ii, int const jj );

int main( int argc, char *argv[] ){
	//================Simulation Constants
	int weightFunction = 2;  	//0/1/2 : NGP/CIC/TSC
	int dim = 2;				
	double L = 10.0;				//Length of box (from -L/2 ~ L/2)
	int Nx = 11;				//Number of grid in x direction. (should be odd number)
	int NParticle=2;//Number of particles used in simulation
	double massParticle=1.0;

	//================Random number generator.
		//To use : d=gsl_rng_uniform(rng);
		gsl_rng *rng;
		rng = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(rng,123456);//The seed is 123456.

	//================Initialize Grid Parameter=========== 
		struct grid2D grid;
		grid.L = L;		
		grid.Nx = Nx;		
		grid.Ny = grid.Nx;
		grid.N = grid.Nx * grid.Ny;
		grid.dx = grid.L / (grid.Nx-1);	//-1 is because boundary (make box closed)
		grid.dy = grid.L / (grid.Ny-1);
		grid.density = (double*) malloc(grid.N*sizeof(double));	//Density on the grid
		grid.phi = (double*) malloc(grid.N*sizeof(double));		//Potential on the grid
		grid.Fx = (double*) malloc(grid.N*sizeof(double));			//Force on grid comes from potential calculation
		grid.Fy = (double*) malloc(grid.N*sizeof(double));

	//================Initialize Particles Parameter======
		struct particle2D myParticle;
		myParticle.number = NParticle;
		myParticle.mass = (double*)malloc(myParticle.number*sizeof(double));
		myParticle.x = (double*)malloc(myParticle.number*sizeof(double));
		myParticle.y = (double*)malloc(myParticle.number*sizeof(double));
		myParticle.Fx = (double*)malloc(myParticle.number*sizeof(double));
		myParticle.Fy = (double*)malloc(myParticle.number*sizeof(double));
		myParticle.vx = (double*)malloc(myParticle.number*sizeof(double));
		myParticle.vy = (double*)malloc(myParticle.number*sizeof(double));
	
	//================Initialize density field=============
		for(int i=0;i<grid.N;i++){
			grid.density[i]=0.0;		
		}

	//Initialize Initial Position of Particle.
		for (int i = 0; i < myParticle.number; ++i){
			myParticle.mass[i]=1.0;
			myParticle.x[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
			myParticle.y[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
			printf("At (%f,%f) \n",myParticle.x[i],myParticle.y[i]);
		}

	//Initialize Force field value
		for(int i=0;i<grid.N;i++){
			grid.Fx[i]=0.0;
			grid.Fy[i]=0.0;		
		}

	//Deposit Particles to grid
	Weight(&grid,&myParticle,weightFunction);
	
	//Print out Density field
		for(int i=0;i<grid.N;i++){
			if(i % grid.Nx == 0){
				printf("\n");
			}
			printf("%.2f\t",grid.density[i]);
		}
		printf("\n");

	//Use Fourier Transform to calculate potential and force.
	poisson_solver_fft_force_2d(dim, &grid);

	//Remap the force to particle.
	WeightForce(&grid,&myParticle,weightFunction);

	
	//Print out the force on a particle.
		for(int i=0;i<myParticle.number;i++){
			printf("(Fx,Fy)=(%.2f,%.2f)\n",myParticle.Fx[i],myParticle.Fy[i]);
		}

	return 0;
}

// void _fftw_2d_(int const dim, struct grid2D *grid) {

//   //double G_const = 6.67408e-8; // #g^-1 s^-2 cm^3 
// 	double G_const = 1.0; // #g^-1 s^-2 cm^3
// 	int Nx = grid->Nx;
// 	int Ny = grid->Ny; 
// 	int total_n = Nx * Ny;
// 	double dNx = (double) (Nx), dNy = (double) (Ny); // default Nx=Ny
// 	int ii, jj;

// 	fftw_complex *sigma_a, *fftsigma_a;
// 	fftw_complex *phika, *phia; 
// 	fftw_plan p, q;

// 	sigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
// 	fftsigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
// 	phika = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
// 	phia = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);

// 	for (ii=0; ii < Nx; ii+=1) {
// 		for (jj=0; jj < Ny; jj+=1) {
// 		  sigma_a[ ii + jj * Nx][0] = grid->density[ ii + jj*Nx]; 
// 		  sigma_a[ ii + jj * Nx][1] = 0.;
// 		}
// 	}


// 	/////////// fft ///////////
// 	p = fftw_plan_dft_2d(Nx, Ny, sigma_a, fftsigma_a, FFTW_FORWARD, FFTW_ESTIMATE);
// 	fftw_execute(p);
// 	fftw_destroy_plan(p);
// 	double kxx, kyy;
// 	for (ii=0; ii < Nx; ii+=1){
// 		for (jj=0; jj < Ny; jj+=1){
// 			kxx = pow((double)(ii+1)*2.*M_PI/dNx, 2.);
// 			kyy = pow((double)(jj+1)*2.*M_PI/dNy, 2.);
// 			phika[ii + jj*Nx][0] = -4. * M_PI * G_const * fftsigma_a[ii + jj * Nx][0] / ((kxx+kyy));
// 			phika[ii + jj*Nx][1] = -4. * M_PI * G_const * fftsigma_a[ii + jj * Nx][1] / ((kxx+kyy));
// 		}
// 	}
// 	/////////// inverse fft ///////////
// 	q = fftw_plan_dft_2d(Nx, Ny, phika, phia, FFTW_BACKWARD, FFTW_ESTIMATE);
// 	fftw_execute(q);
// 	fftw_destroy_plan(q);
// 	/////////// normalization ///////////
// 	for (ii=0; ii < Nx; ii+=1){
// 		for (jj=0; jj < Ny; jj+=1){
// 		  phia[ii + jj * Nx][0] /= grid->N;
// 		  phia[ii + jj * Nx][1] /= grid->N;
// 		}
// 	}
	
// 	////////// magnitude of potential //////////
// 	for (ii=0; ii < Nx; ii+=1){
// 		for (jj=0; jj < Ny; jj+=1){ 
// 		  grid->phi[ii + jj * Nx] = sqrt( pow(phia[ii + jj * Nx][0],2) + pow(phia[ii + jj * Nx][1],2) );
// 		}
// 	}  

// 	for (ii=0; ii < Nx; ii+=1){
// 		for (jj=0; jj < Ny; jj+=1){     
// 		  _2nd_order_diff_(grid, ii, jj);
// 		}
// 	}

// 	// for (ii=0; ii < Nx; ii+=1) {
// 	// 	for (jj=0; jj < Ny; jj+=1){
// 	// 	  printf("recover: %3d %3d %+.5e %+.5e I vs. %+.5e %+.5e I \n",
// 	// 	      ii, jj, sigma_a[ii + jj*Nx][0], sigma_a[ii + jj*Nx][1], 
// 	// 	      phia[ii + jj*Nx][0], phia[ii + jj*Nx][1]);    
// 	// 	}
// 	// }

// 	// for (ii=0; ii<Nx; ii+=1) {
// 	// 	for (jj=0; jj<Ny; jj+=1){
// 	// 	  printf("ii:%3d, jj:%3d, Fx:%.5e Fy:%.5e\n",
// 	// 	      ii, jj, Fx[ii + jj*Nx], Fy[ii + jj*Nx]);    
// 	// 	}
// 	// }

// 	fftw_cleanup();
// }

// void _2nd_order_diff_(struct grid2D *grid, int const ii, int const jj ) {

// 	double factor1 = -1./(1.+1.)/grid->dx;
// 	double factor2 = -1./(1.+1.)/grid->dy;
// 	int Nx = grid->Nx;
// 	int Ny = grid->Ny; 

// 	int l,r,t,d;	//Index of left , right , up and down.
// 	l = (Nx+ii+1)%Nx;
// 	r = (Nx+ii-1)%Nx;
// 	t = (Ny+jj+1)%Ny;
// 	d = (Ny+jj-1)%Ny;

// 	grid->Fx[ ii + jj*Nx] = factor1*(grid->phi[l + jj * Nx ] - grid->phi[r + jj * Nx] );
// 	grid->Fy[ ii + jj*Nx] = factor2*(grid->phi[ii + t * Nx ] - grid->phi[ii + d * Nx] );  

// }