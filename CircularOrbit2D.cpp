//To compile should add "-lfftw3 -lgsl"
#include <stdio.h>
#include <math.h>
// #include <gsl/gsl_rng.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>
#include <fftw3.h>
#include <complex.h>
using namespace std;

struct particle2D{
	int number;
	double *mass;
	double *x;
	double *y;
	double *Fx;
	double *Fy;
	double *vx;
	double *vy;
};
struct grid2D{
	double L;
	int Nx;
	int Ny;
	int N;
	double dx;
	double dy;
	double *density;
	double *phi;
	double *Fx;
	double *Fy;
};
void Weight(struct grid2D *grid,struct particle2D *particle,int type){

	//Initialize Density field
	for(int i=0 ; i < grid->N ; i++){
		grid->density[i]=0.0;
	}

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
			grid->density[(lx+sx)*grid->Nx+(ly+sy) ] += particle->mass[i];
		}else if(type == 1){
			//CIC
			//xx
			//ox
			for(int j=0;j<4;j++){
				int p = j / 2;
				int q = j % 2;
				double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
				grid->density[(lx+p)*grid->Nx+(ly+q)] += particle->mass[i] * wFactor;
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
					grid->density[indx*grid->Nx+indy]+=weightX[xx+1]*weightY[yy+1]*particle->mass[i];
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
			particle->Fx[i]=grid->Fx[(lx+sx)*grid->Nx+(ly+sy)]*particle->mass[i];
			particle->Fy[i]=grid->Fy[(lx+sx)*grid->Nx+(ly+sy)]*particle->mass[i];
		}else if(type == 1){
			for(int j=0;j<4;j++){
				int p = j / 2;
				int q = j % 2;
				double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
				particle->Fx[i] += grid->Fx[(lx+p)*grid->Nx+(ly+q)] * wFactor*particle->mass[i];
				particle->Fy[i] += grid->Fy[(lx+p)*grid->Nx+(ly+q)] * wFactor*particle->mass[i];
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
					particle->Fx[i] += grid->Fx[indx*grid->Nx+indy]*weightX[xx+1]*weightY[yy+1]*particle->mass[i];
					particle->Fy[i] += grid->Fy[indx*grid->Nx+indy]*weightX[xx+1]*weightY[yy+1]*particle->mass[i];
				}
			}
		}

		
	}
}

void _fftw_2d_(int const dim, struct grid2D *grid);

void _2nd_order_diff_(struct grid2D *grid, int const ii, int const jj );

void isolatedPotential(struct grid2D *grid);

void kick(struct particle2D *particle , double dt){
	double ax,ay;
	for(int i=0 ; i<particle->number ; i++){
		//Cauculate the acceleration of each particle.
		ax = particle->Fx[i] / particle->mass[i];
		ay = particle->Fy[i] / particle->mass[i];

		//Calculate velocity of each particle.
		particle->vx[i] += ax * dt;
		particle->vy[i] += ay * dt;
	}
}
void drift(struct particle2D *particle , double dt){
	for(int i=0 ; i<particle->number ; i++){
		particle->x[i] += particle->vx[i] * dt;
		particle->y[i] += particle->vy[i] * dt;
	}
}

int main( int argc, char *argv[] ){
	//================Simulation Constants
	int weightFunction = 1;  	//0/1/2 : NGP/CIC/TSC
	int orbitIntegrattion = 0;	//0/1/2 : KDK/DKD/RK4
	int dim = 2;				
	double L = 10.0;				//Length of box (from -L/2 ~ L/2)
	int Nx = 101;				//Number of grid in x direction. (should be odd number)
	int NParticle=2;//Number of particles used in simulation
	double massParticle=1.0;
	double dt = 1.0e-2;
	double G = 1.0;

	//Output to a file
	FILE *output;
	output = fopen("result.txt","w");
	//================Random number generator.
		//To use : d=gsl_rng_uniform(rng);
		// gsl_rng *rng;
		// rng = gsl_rng_alloc(gsl_rng_mt19937);
		// gsl_rng_set(rng,123456);//The seed is 123456.

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

	//Initialize mass of particles
		myParticle.mass[0]=2.0;
		myParticle.mass[1]=1.0;
	//Initialize Initial Position of Particle.
		// for (int i = 0; i < myParticle.number; ++i){
		// 	myParticle.x[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
		// 	myParticle.y[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
		// 	printf("At (%f,%f) \n",myParticle.x[i],myParticle.y[i]);
		// }
		myParticle.x[0] = 1.0;
		myParticle.y[0] = 0.0;
		myParticle.x[1] = -2.0;
		myParticle.y[1] = 0.0;
	//Initialize Force field value
		for(int i=0;i<grid.N;i++){
			grid.Fx[i]=0.0;
			grid.Fy[i]=0.0;		
		}

	Weight(&grid,&myParticle,weightFunction);
		
	isolatedPotential(&grid);				
		
	WeightForce(&grid,&myParticle,weightFunction);

	//Initialize Initial velocity
		myParticle.vx[0] = 0.0;
		myParticle.vy[0] = -sqrt(fabs(myParticle.Fx[0]*1.0)/myParticle.mass[0]);
		myParticle.vx[1] = 0.0;
		myParticle.vy[1] = sqrt(fabs(myParticle.Fx[1]*2.0)/myParticle.mass[1]);
		printf("%f\t%f\n",myParticle.Fx[0],myParticle.Fx[1]);

			

	//Time evolution loop
	double t = 0.0;
	for(int st=0;st<1000;st++){
		//Deposit Particles to grid
		Weight(&grid,&myParticle,weightFunction);

		//Use Fourier Transform to calculate potential and force.
		//_fftw_2d_(dim, &grid);
		isolatedPotential(&grid);

		//Remap the force to particle.
		WeightForce(&grid,&myParticle,weightFunction);

		//Move particle
		if(orbitIntegrattion == 0){
			//KDK scheme
			kick(&myParticle,dt/2);

			drift(&myParticle,dt);
			
			
			Weight(&grid,&myParticle,weightFunction);
			isolatedPotential(&grid);
			WeightForce(&grid,&myParticle,weightFunction);

			kick(&myParticle,dt/2);
		}else if(orbitIntegrattion == 1){
			//DKD scheme
			drift(&myParticle,dt/2);

			Weight(&grid,&myParticle,weightFunction);
			isolatedPotential(&grid);
			WeightForce(&grid,&myParticle,weightFunction);

			kick(&myParticle,dt);
			drift(&myParticle,dt/2);
		}else if(orbitIntegrattion == 2){
			//RK4
			//To-do
		}

		//print out the position of particle 1
		if(st % 1 == 0){
			fprintf(output,"%f\t%f\t",myParticle.x[0],myParticle.y[0]);
			fprintf(output,"%f\t%f\n",myParticle.x[1],myParticle.y[1]);
		}
		
	}

	
	
	// //Print out Density field
	// 	for(int i=0;i<grid.N;i++){
	// 		if(i % grid.Nx == 0){
	// 			printf("\n");
	// 		}
	// 		printf("%.2f\t",grid.density[i]);
	// 	}
	// 	printf("\n");

	// //Print out the force on a particle.
	// 	for(int i=0;i<myParticle.number;i++){
	// 		printf("(Fx,Fy)=(%.2f,%.2f)\n",myParticle.Fx[i],myParticle.Fy[i]);
	// 	}
	fclose(output);
	return 0;
}

void _fftw_2d_(int const dim, struct grid2D *grid) {

  //double G_const = 6.67408e-8; // #g^-1 s^-2 cm^3 
	double G_const = 1.0; // #g^-1 s^-2 cm^3
	int Nx = grid->Nx;
	int Ny = grid->Ny; 
	int total_n = Nx * Ny;
	double dNx = (double) (Nx), dNy = (double) (Ny); // default Nx=Ny
	int ii, jj;

	fftw_complex *sigma_a, *fftsigma_a;
	fftw_complex *phika, *phia; 
	fftw_plan p, q;

	sigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
	fftsigma_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
	phika = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);
	phia = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid->N);

	for (ii=0; ii < Nx; ii+=1) {
		for (jj=0; jj < Ny; jj+=1) {
		  sigma_a[ ii * Nx + jj][0] = grid->density[ ii * Nx + jj]; 
		  sigma_a[ ii * Nx + jj][1] = 0.;
		}
	}


	/////////// fft ///////////
	p = fftw_plan_dft_2d(Nx, Ny, sigma_a, fftsigma_a, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	double kxx, kyy;
	for (ii=0; ii < Nx; ii+=1){
		for (jj=0; jj < Ny; jj+=1){
			if(ii != 0 || jj != 0){
				kxx = pow((double)(ii)*2.*M_PI/grid->L, 2.);
				kyy = pow((double)(jj)*2.*M_PI/grid->L, 2.);
				phika[ii * Nx + jj][0] = -4. * M_PI * G_const * fftsigma_a[ii * Nx + jj][0] / ((kxx+kyy));
				phika[ii * Nx + jj][1] = -4. * M_PI * G_const * fftsigma_a[ii * Nx + jj][1] / ((kxx+kyy));
			}
			
		}
	}
	
	/////////// inverse fft ///////////
	q = fftw_plan_dft_2d(Nx, Ny, phika, phia, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(q);
	fftw_destroy_plan(q);
	/////////// normalization ///////////
	for (ii=0; ii < Nx; ii+=1){
		for (jj=0; jj < Ny; jj+=1){
		  phia[ii * Nx + jj][0] /= grid->N;
		  phia[ii * Nx + jj][1] /= grid->N;
		}
	}
	
	////////// magnitude of potential //////////
	for (ii=0; ii < Nx; ii+=1){
		for (jj=0; jj < Ny; jj+=1){ 
		  grid->phi[ii * Nx + jj] = -sqrt( pow(phia[ii * Nx + jj][0],2) + pow(phia[ii * Nx + jj][1],2) );
		}
	}  

	for (ii=0; ii < Nx; ii+=1){
		for (jj=0; jj < Ny; jj+=1){     
		  _2nd_order_diff_(grid, ii, jj);
		}
	}

	fftw_cleanup();
}

void _2nd_order_diff_(struct grid2D *grid, int const ii, int const jj ) {

	double factor1 = -1./(1.+1.)/grid->dx;
	double factor2 = -1./(1.+1.)/grid->dy;
	int Nx = grid->Nx;
	int Ny = grid->Ny; 

	int l,r,u,d;	//Index of left , right , up and down.
	l = (Nx+ii-1)%Nx;
	r = (Nx+ii+1)%Nx;
	u = (Ny+jj+1)%Ny;
	d = (Ny+jj-1)%Ny;

	grid->Fx[ ii * Nx + jj] = factor1*(grid->phi[r * Nx + jj] - grid->phi[l * Nx + jj] );
	grid->Fy[ ii * Nx + jj] = factor2*(grid->phi[ii * Nx + u] - grid->phi[ii * Nx + d] );  

}

void isolatedPotential(struct grid2D *grid){
	
	double *densityPad,*greenFunction;
	int N = grid->N;
	int Nx = grid->Nx;
	int Ny = grid->Ny;

	int NNx = 2*Nx;
	int NNy = 2*Ny;

	densityPad = (double*)malloc(4*N*sizeof(double));
	greenFunction = (double*)malloc(4*N*sizeof(double));

	fftw_complex *dp, *fftdp;	//For Padding density
	fftw_complex *gf, *fftgf; 	//For green function
	fftw_complex *phi,*ifftphi;	//For potential
	fftw_plan p1 , p2, q;

	dp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N);
	fftdp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N);
	gf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N);
	fftgf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N);
	phi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N);
	ifftphi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4*N);



	

	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			if(i < Nx && j < Ny){
				//Copy initial density
				densityPad[i*NNx+j] = grid->density[i*Nx+j];
			}else{
				//Padding 0s
				densityPad[i*NNx+j] = 0.0;
			}
		}
	}
	//Initialize green function array
	for(int i=0;i<2*Nx;i++){
		for(int j=0;j<2*Ny;j++){
			greenFunction[i*Nx+j] = 0.0;			
		}
	}
	//Construct green function
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			if(i != 0 || j != 0){
				greenFunction[i*NNx+j] = -1.0/sqrt(pow(i,2)+pow(j,2));
			}
		}
	}

	for(int i=Nx+1;i<NNx;i++){
		for(int j=0;j<Ny;j++){
			greenFunction[i*NNx+j] = greenFunction[(NNx-i)*NNx+j];
		}
	}
	for(int i=0;i<Nx;i++){
		for(int j=Ny+1;j<NNy;j++){
			greenFunction[i*NNx+j] = greenFunction[i*NNx+(NNy-j)];
		}
	}

	for(int i=Nx+1;i<NNx;i++){
		for(int j=Ny+1;j<NNy;j++){
			greenFunction[i*NNx+j] = greenFunction[(NNx-i)*NNx+(NNy-j)];
		}
	}

	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			dp[i*NNx+j][0] = densityPad[i*NNx+j];
			dp[i*NNx+j][1] = 0.0;
			gf[i*NNx+j][0] = greenFunction[i*NNx+j];
			gf[i*NNx+j][1] = 0.0;
		}
	}

	//fft
	p1 = fftw_plan_dft_2d(NNx, NNy, dp, fftdp, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p1);
	fftw_destroy_plan(p1);

	p2 = fftw_plan_dft_2d(NNx, NNy, gf, fftgf, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p2);
	fftw_destroy_plan(p2);

	for(int i=0;i<2*Nx;i++){
		for(int j=0;j<2*Ny;j++){
			//Multiply 2 imaginary numbers
			phi[i*NNx+j][0] = fftdp[i*NNx+j][0] * fftgf[i*NNx+j][0] - fftdp[i*NNx+j][1] * fftgf[i*NNx+j][1];
			phi[i*NNx+j][1] = fftdp[i*NNx+j][0] * fftgf[i*NNx+j][1] + fftdp[i*NNx+j][1] * fftgf[i*NNx+j][0];
		}
	}

	//ifft
	q = fftw_plan_dft_2d(NNx, NNy, phi, ifftphi, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(q);
	fftw_destroy_plan(q);

	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			grid->phi[i*Nx+j] = -grid->dx * grid->dy *sqrt(pow(ifftphi[i*NNx+j][0],2)+pow(ifftphi[i*NNx+j][1],2));
		}		
	}

	for (int i=0; i < Nx; i+=1){
		for (int j=0; j < Ny; j+=1){     
		  _2nd_order_diff_(grid, i, j);
		}
	}

	free(densityPad);
	free(greenFunction);
	free(dp);
	free(fftdp);
	free(gf);
	free(fftgf);
	free(phi);
	free(ifftphi);

}