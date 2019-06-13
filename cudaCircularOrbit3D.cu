//To compile should add "-lfftw3 -lgsl"
#include <stdio.h>
#include <math.h>
// #include <gsl/gsl_rng.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>
//#include <fftw3.h>
#include <complex>
#include <cufft.h>
using namespace std;

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
	double *phi;
	double *density;
	double *Fx;
	double *Fy;
	double *Fz;
};
struct rk43D{
    int step1;
    int step2;
    double *ax;
    double *ay;
	double *az;
    double *vx;
    double *vy;
	double *vz;
};
__global__
void phiToForce(double* phi,double* Fx,double* Fy,double* Fz,int Nx,double L,double nConst){
	int N = Nx*Nx*Nx;
	int const Ny = Nx;
  	int const Nz = Nx;
	double dx = L / Nx;
	double factor = -1./(2.0*dx) * nConst;

	int index = blockDim.x * blockIdx.x + threadIdx.x;
    
    while(index < N){
    	int ii = index / (Ny*Nz);
    	int jj = (index / Nz) % Ny;
    	int kk = index % Nz;
    	Fx[index] = factor*( phi[ ( (ii+1)%Nx )*Ny*Nz + jj*Nz + kk ] 
                              - phi[ ( (Nx+ii-1)%Nx )*Ny*Nz + jj*Nz + kk ] );

    	Fy[index] = factor*( phi[ ii*Ny*Nz + ( (jj+1)%Ny )*Nz + kk ] 
                              - phi[ ii*Ny*Nz + ( (Ny+jj-1)%Ny )*Nz + kk ] );  

  		Fz[index] = factor*( phi[ ii*Ny*Nz  + jj*Nz + ((kk+1)%Nz)] 
                              - phi[ ii*Ny*Nz  + jj*Nz + ((Nz+kk-1)%Nz)]);
    	index += blockDim.x*gridDim.x;
    }

}
__global__
void overk2(double *out,int Nx,int Ny,int Nzh){
	//Notice that *out actullay is a compelx <double> pointer array.

	int N = Nx*Ny*Nzh;
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int fi,fj;
	double kxx,kyy,kzz;

    while(index < N){
    	int ii = index / (Ny*Nzh);
    	int jj = (index / Nzh) % Ny;
    	int kk = index % Nzh;

    	if (2*ii < Nx) {fi = ii;}
		else           {fi = Nx-ii;}
		if (2*jj < Ny) {fj = jj;}
		else           {fj = Ny-jj;}

		kxx = 1.0*fi*fi;
		kyy = 1.0*fj*fj;
		kzz = 1.0*kk*kk;

		if(index != 0){
			//ii != 0 || jj != 0 || kk!=0
			out[2*index] = out[2*index] / (kxx+kyy+kzz);		//real part
			out[2*index+1] = out[2*index+1] / (kxx+kyy+kzz);	//imaginary part
		}   	

    	index += blockDim.x*gridDim.x;
    }
}

void Weight(struct grid3D *grid,struct particle3D *particle,int type);
void WeightForce(struct grid3D *grid,struct particle3D *particle,int type);
void poisson_solver_fft_force_3d(int const dim, struct grid3D *grid);
void _2nd_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk );
void _2nd_order_diff_3d_cuda(struct grid3D *grid);
void calculateGreenFFT(struct grid3D *grid, complex<double>* fftgf);
void isolatedPotential(struct grid3D *grid, complex<double>* fftgf);


void kick(struct particle3D *particle , double dt);
void drift(struct particle3D *particle , double dt);
void init_rk4(struct particle3D *particle, struct particle3D *buff, struct rk43D *rk4);
void rk4_mid(struct particle3D *particle, struct particle3D *buff, struct rk43D *rk4, double dt, int weighting);
void rk4_end(struct particle3D *particle, struct rk43D *rk4, double dt);
void periodic_boundary(double position, double length);
void boundary_check(int boundary, struct particle3D *particle, double L);

//Functions to locate memory and free memory of different struct.
void locateMemoryParticle(struct particle3D *particle,int N);
void freeMemoryParticle(struct particle3D *particle);
void locateMemoryRk4(struct rk43D *rk4,int N);
void freeMemoryRk4(struct rk43D *rk4);
void locateMemoryGrid(struct grid3D *grid);
void freeMemoryGrid(struct grid3D *grid);


int main( int argc, char *argv[] ){
	//================Simulation Constants
	int weightFunction = 1;  	//0/1/2 : NGP/CIC/TSC
	int orbitIntegration = 1;	//0/1/2 : KDK/DKD/RK4
	int poissonSolver = 0;		//0/1   : fft/isolated
	int boundary = 0;           //0/1/2 : periodic/isolated/no boundary
	int dim = 3;				
	double L = 10.0;				//Length of box (from -L/2 ~ L/2)
	int Nx = 256;				//Number of grid in x direction. (should be odd number)
	int NParticle=2;//Number of particles used in simulation
	//double massParticle=1.0;
	double dt = 1.0e-2;
	double G = 1.0;
	double T,r1_0,r2_0;
	cudaEvent_t start, stop;	//For cuda timing
	float totalTime;
	
	//================Structs
	struct grid3D grid;
	struct particle3D myParticle;
	struct particle3D buffParticle;		//If not RK4 mode , it will not be malloc and free
	struct rk43D myrk4;					//If not RK4 mode , it will not be malloc and free
	complex<double>* fftgf;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
 	cudaEventRecord(start,0);

 	

	//Output to a file
	FILE *output;
	output = fopen("result.txt","w");
	//================Random number generator.
		//To use : d=gsl_rng_uniform(rng);
		// gsl_rng *rng;
		// rng = gsl_rng_alloc(gsl_rng_mt19937);
		// gsl_rng_set(rng,123456);//The seed is 123456.

	//================Initialize Grid Parameter=========== 
		
		grid.L = L;		
		grid.Nx = Nx;		
		grid.Ny = grid.Nx;
		grid.Nz = grid.Nx;
		grid.N = grid.Nx * grid.Ny * grid.Nz;
		grid.dx = grid.L / (grid.Nx-1);	//-1 is because boundary (make box closed)
		grid.dy = grid.L / (grid.Ny-1);
		grid.dz = grid.L / (grid.Nz-1);
		locateMemoryGrid(&grid);

	//================Initialize Particles ======
		
		myParticle.number = NParticle;
		locateMemoryParticle(&myParticle,NParticle);

		if(orbitIntegration == 2){
	//================Initialize Particles for RK4 (buffer)======
			buffParticle.number = NParticle;
			locateMemoryParticle(&buffParticle,NParticle);
	//================Initialize Runge-Kutta Coefficient======
			myrk4.step1 = 0;
			myrk4.step2 = 1;
			locateMemoryRk4(&myrk4,NParticle);
		}

		if(poissonSolver == 1){
			fftgf = (complex<double>*) malloc(sizeof(complex<double>) * 8*grid.N);
			calculateGreenFFT(&grid,fftgf);
		}
      
	
	//Initialize mass of particles
		myParticle.mass[0]=2.0;
		myParticle.mass[1]=1.0;

		if(orbitIntegration == 2){
			buffParticle.mass[0] = myParticle.mass[0];
			buffParticle.mass[1] = myParticle.mass[1];
		}
	//Initialize Initial Position of Particle.
		// for (int i = 0; i < myParticle.number; ++i){
		// 	myParticle.x[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
		// 	myParticle.y[i]=gsl_rng_uniform(rng) * grid.L - grid.L/2;
		// 	printf("At (%f,%f) \n",myParticle.x[i],myParticle.y[i]);
		// }
		myParticle.x[0] = 1.0;
		myParticle.y[0] = 0.0;
		myParticle.z[0] = 0.0;
		myParticle.x[1] = -2.0;
		myParticle.y[1] = 0.0;
		myParticle.z[1] = 0.0;
		r1_0 = sqrt(pow(myParticle.x[0],2)+pow(myParticle.y[0],2)+pow(myParticle.z[0],2));
		r2_0 = sqrt(pow(myParticle.x[1],2)+pow(myParticle.y[1],2)+pow(myParticle.z[1],2));

		Weight(&grid,&myParticle,weightFunction);
		
		if(poissonSolver == 0){poisson_solver_fft_force_3d(dim,&grid);}
		else if ( poissonSolver == 1 ){isolatedPotential(&grid,fftgf);}
		
		WeightForce(&grid,&myParticle,weightFunction);

	// //Initialize Initial velocity
		myParticle.vx[0] = 0.0;
		myParticle.vy[0] = -sqrt(fabs(myParticle.Fx[0]*myParticle.x[0])/myParticle.mass[0]);
		myParticle.vz[0] = 0.0;
		myParticle.vx[1] = 0.0;
		myParticle.vy[1] = sqrt(fabs(myParticle.Fx[1]*myParticle.x[1])/myParticle.mass[1]);
		myParticle.vz[1] = 0.0;

		double F_0;
		F_0 = G * myParticle.mass[0] * myParticle.mass[1]/pow(r1_0+r2_0,2);
		T =sqrt(myParticle.mass[0]*4*pow(M_PI,2)*r1_0/F_0);

		//Check whether force is same magnitude but inverse direction.
		printf("%f\t%f\t%f\n",myParticle.Fx[0],myParticle.Fy[0],myParticle.Fz[0]);
		printf("%f\t%f\t%f\n",myParticle.Fx[1],myParticle.Fy[1],myParticle.Fz[1]);
		
	//Time evolution loop
	double t = 0.0;
	for(int st=0;st < 1000;st++){
	 	//Deposit Particles to grid
	 	Weight(&grid,&myParticle,weightFunction);

	 	//Use Fourier Transform to calculate potential and force.
		if ( poissonSolver == 0 )	poisson_solver_fft_force_3d(dim, &grid);
		else if ( poissonSolver == 1 ){isolatedPotential(&grid,fftgf);}

 		//Remap the force to particle.
		WeightForce(&grid,&myParticle,weightFunction);

      

 		//Move particle
		if(orbitIntegration == 0){
 		//KDK scheme
			kick(&myParticle,dt/2);
 			drift(&myParticle,dt);
			boundary_check(boundary, &myParticle, L);
 			Weight(&grid,&myParticle,weightFunction);
			if ( poissonSolver == 0 )       poisson_solver_fft_force_3d(dim, &grid);
			else if ( poissonSolver == 1 ){isolatedPotential(&grid,fftgf);}
			WeightForce(&grid,&myParticle,weightFunction);

 			kick(&myParticle,dt/2);
		}
		else if(orbitIntegration == 1){
 			//DKD scheme
 			drift(&myParticle,dt/2);
			boundary_check(boundary, &myParticle, L);
 			Weight(&grid,&myParticle,weightFunction);
			if ( poissonSolver == 0 )       poisson_solver_fft_force_3d(dim, &grid);
			else if ( poissonSolver == 1 )  isolatedPotential(&grid,fftgf);
			WeightForce(&grid,&myParticle,weightFunction);

 			kick(&myParticle,dt);
 			drift(&myParticle,dt/2);
 		}
		else if(orbitIntegration == 2){
 			//RK4	
			init_rk4(&myParticle,&buffParticle,&myrk4);
			
			rk4_mid(&myParticle,&buffParticle,&myrk4,dt/2,1);        //k1
			boundary_check(boundary, &buffParticle, L);
			Weight(&grid,&buffParticle,weightFunction);
			if ( poissonSolver == 0 )       poisson_solver_fft_force_3d(dim, &grid);
			else if ( poissonSolver == 1 )  isolatedPotential(&grid,fftgf);
			WeightForce(&grid,&buffParticle,weightFunction);
			
			rk4_mid(&myParticle,&buffParticle,&myrk4,dt/2,2);        //k2
			boundary_check(boundary, &buffParticle, L);
			Weight(&grid,&buffParticle,weightFunction);
			if ( poissonSolver == 0 )       poisson_solver_fft_force_3d(dim, &grid);
			else if ( poissonSolver == 1 )  isolatedPotential(&grid,fftgf);
			WeightForce(&grid,&buffParticle,weightFunction);
			
			rk4_mid(&myParticle,&buffParticle,&myrk4,dt,2);          //k3
			boundary_check(boundary, &buffParticle, L);
			Weight(&grid,&buffParticle,weightFunction);
			if ( poissonSolver == 0 )       poisson_solver_fft_force_3d(dim, &grid);
			else if ( poissonSolver == 1 )  isolatedPotential(&grid,fftgf);
			WeightForce(&grid,&buffParticle,weightFunction);
			
			rk4_mid(&myParticle,&buffParticle,&myrk4,dt,1);          //k4
			
			rk4_end(&myParticle,&myrk4,dt);

 		}
		//Boundary Condition
		#ifdef CHECK_THIS_LATER
		if (boundary == 0){
			for (int i=0; i<NParticle; i++){
			       	if ( abs(myParticle.x[i]) > L/2){
					periodic_boundary(myParticle.x[i],L);
				}
				if ( abs(myParticle.y[i]) > L/2){
					periodic_boundary(myParticle.y[i],L);
				}
				if ( abs(myParticle.z[i]) > L/2){
					periodic_boundary(myParticle.z[i],L);
				}
			}	
		}
		else if (boundary == 1){
			for (int i=0; i<NParticle; i++){
				if ( abs(myParticle.x[i]) > L/2 || abs(myParticle.y[i]) > L/2 || abs(myParticle.z[i]) > L/2){
					myParticle.mass[i] = 0;
					cout << "A particle reaches the boundary." << endl;
				}
			}
		}
 		#endif
		boundary_check(boundary, &myParticle, L);

 		//print out the position of particle 1
		if(st % 20 == 0){
			printf("Step:%d\n", st);
			// double momentum_x = 0;
			// double momentum_y = 0;
			// double momentum_z = 0;

   //      		for (int i=0; i<NParticle; i++){
			// 	momentum_x += myParticle.mass[i] * myParticle.vx[i];
			// 	momentum_y += myParticle.mass[i] * myParticle.vy[i];
			// 	momentum_z += myParticle.mass[i] * myParticle.vz[i];
			// }
			// cout << "(px , py,  pz) = (" << momentum_x << ", " << momentum_y << ", " << momentum_z << ")" << endl;
			fprintf(output,"%f\t%f\t",myParticle.x[0],myParticle.y[0]);
			fprintf(output,"%f\t%f\n",myParticle.x[1],myParticle.y[1]);
	
		}
		t+=dt;
 	}
 	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&totalTime, start, stop);
	printf("Total time : %f\n",totalTime/1000);
		
 	fclose(output);
	freeMemoryGrid(&grid);
	freeMemoryParticle(&myParticle);
	if(orbitIntegration == 2){
		freeMemoryParticle(&buffParticle);
	}
	if(poissonSolver == 1){
		free(fftgf);
	}
	return 0;
}

void poisson_solver_fft_force_3d(int const dim, struct grid3D *grid){
	
	int const Nx = grid->Nx;
	int const Ny = grid->Ny; 
	int const Nz = grid->Nz;
	int const Nzh = (Nz/2+1);
	cufftHandle p1,p2;
	
	double *in;
	complex<double> *out;
	in = (double*) malloc( sizeof(double) * Nx*Ny*Nz );
	out = (complex<double>*) malloc( sizeof(complex<double>) * Nx*Ny*Nzh);
		

	

	/////////// fft ///////////
	cufftDoubleReal *dataIn;
	cufftDoubleComplex *dataOut;
    cudaMalloc((void**)&dataIn, sizeof(cufftDoubleReal)*Nx*Ny*Nz);
    cudaMalloc((void**)&dataOut, sizeof(cufftDoubleComplex)*Nx*Ny*Nzh);
    cudaMemcpy(dataIn, grid->density, sizeof(cufftDoubleReal)*Nx*Ny*Nz, cudaMemcpyHostToDevice);
	
	//cufft
	if (cufftPlan3d(&p1,Nx,Ny,Nz, CUFFT_D2Z) != CUFFT_SUCCESS) {
        printf("CUFFT error: Plan D2Z creation failed.\n");
		exit(1);
    }
    if (cufftExecD2Z(p1, dataIn, dataOut) != CUFFT_SUCCESS) {
		printf("CUFFT error: ExecD2Z forward failed.\n");
		exit(1);
    }
    cudaMemcpy(out, dataOut,sizeof(cufftDoubleComplex)*Nx*Ny*Nzh, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(p1);
    

	double *d_out;
	cudaMalloc((void**)&d_out, sizeof(double)*2*Nx*Ny*Nzh);
	cudaMemcpy(d_out,out,sizeof(double)*2*Nx*Ny*Nzh,cudaMemcpyHostToDevice);

	overk2 <<<128,128>>> (d_out,Nx,Ny,Nzh);

	cudaMemcpy(out,d_out,sizeof(double)*2*Nx*Ny*Nzh,cudaMemcpyDeviceToHost);

	cudaFree(d_out);
	


	/////////// inverse fft ///////////
    cudaMemcpy(dataOut, out,sizeof(cufftDoubleComplex)*Nx*Ny*Nzh, cudaMemcpyHostToDevice);
	
	//cufft
	if (cufftPlan3d(&p2,Nx,Ny,Nz, CUFFT_Z2D) != CUFFT_SUCCESS) {
        printf("CUFFT error: Plan creation failed.\n");
		exit(1);
    }
    if (cufftExecZ2D(p2, dataOut, dataIn) != CUFFT_SUCCESS) {
		printf("CUFFT error: ExecZ2D  failed.\n");
		exit(1);
    }
    cudaMemcpy(in, dataIn, sizeof(cufftDoubleReal)*Nx*Ny*Nz, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(p2);


	/////////// normalization ///////////

	double nConst = -1.0 / M_PI/grid->L;	//Normalization constant.

	double* d_in;
	double* d_Fx; 
	double* d_Fy;
	double* d_Fz;

	int size = grid->N*sizeof(double);

	cudaMalloc((void**)&d_in, size);
	cudaMalloc((void**)&d_Fx, size);
	cudaMalloc((void**)&d_Fy, size);
	cudaMalloc((void**)&d_Fz, size);

	cudaMemcpy(d_in,in,size,cudaMemcpyHostToDevice);
	phiToForce <<<128,128>>>(d_in,d_Fx,d_Fy,d_Fz,grid->Nx,grid->L,nConst);
	cudaDeviceSynchronize();
	cudaMemcpy(grid->Fx,d_Fx,size, cudaMemcpyDeviceToHost);
	cudaMemcpy(grid->Fy,d_Fy,size, cudaMemcpyDeviceToHost);
	cudaMemcpy(grid->Fz,d_Fz,size, cudaMemcpyDeviceToHost);

	cudaFree(d_in);
	cudaFree(d_Fx);
	cudaFree(d_Fy);
	cudaFree(d_Fz);
	cudaFree(dataIn);
	cudaFree(dataOut);
	

	free(in);
	free(out);
	
}


void _2nd_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk ) {

  double factor1 = -1./(2.*grid->dx);
  double factor2 = -1./(2.*grid->dy);
  double factor3 = -1./(2.*grid->dz);
  int const Nx = grid->Nx;
  int const Ny = grid->Ny;
  int const Nz = grid->Nz;
  int index = ii*Ny*Nz + jj*Nz + kk;

  grid->Fx[ index ] = factor1*( grid->phi[ ( (Nx+ii+1)%Nx )*Nx*Ny + jj*Nx + kk ] 
                              - grid->phi[ ( (Nx+ii-1)%Nx )*Nx*Ny + jj*Nx + kk ] );

  grid->Fy[ index ] = factor2*( grid->phi[ ii*Nx*Ny + ( (Ny+jj+1)%Ny )*Nx + kk ] 
                              - grid->phi[ ii*Nx*Ny + ( (Ny+jj-1)%Ny )*Nx + kk ] );  

  grid->Fz[ index ] = factor3*( grid->phi[ ii*Nx*Ny  + jj*Nx + ((Nz+kk+1)%Nz)] 
                              - grid->phi[ ii*Nx*Ny  + jj*Nx + ((Nz+kk-1)%Nz)]);

}
void _2nd_order_diff_3d_cuda(struct grid3D *grid) {
	double* d_phi;
	double* d_Fx; 
	double* d_Fy;
	double* d_Fz;

	int size = grid->N * sizeof(double);

	cudaMalloc((void**)&d_phi, size);
	cudaMalloc((void**)&d_Fx, size);
	cudaMalloc((void**)&d_Fy, size);
	cudaMalloc((void**)&d_Fz, size);

	cudaMemcpy(d_phi,grid->phi,size,cudaMemcpyHostToDevice);
	phiToForce <<<128,128>>>(d_phi,d_Fx,d_Fy,d_Fz,grid->Nx,grid->L,1.0);
	cudaDeviceSynchronize();
	cudaMemcpy(grid->Fx,d_Fx,size, cudaMemcpyDeviceToHost);
	cudaMemcpy(grid->Fy,d_Fy,size, cudaMemcpyDeviceToHost);
	cudaMemcpy(grid->Fz,d_Fz,size, cudaMemcpyDeviceToHost);

	cudaFree(d_phi);
	cudaFree(d_Fx);
	cudaFree(d_Fy);
	cudaFree(d_Fz);
}

void calculateGreenFFT(struct grid3D *grid, complex<double> *fftgf){
	double *greenFunction;
	int N = grid->N;
	int Nx = grid->Nx;
	int Ny = grid->Ny;
	int Nz = grid->Nz;

	int NNx = 2*Nx;
	int NNy = 2*Ny;
	int NNz = 2*Nz;

	greenFunction = (double*)malloc(8*N*sizeof(double));
	

	complex<double> *gf; 	//For green function
	gf = (complex<double>*) malloc(sizeof(complex<double>) * 8*N);
	

	//Initialize green function array
	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			for(int k=0;k<NNz;k++){
				greenFunction[i*NNx*NNy+j*NNx+k] = 0.0;	
			}		
		}
	}
	//Construct green function
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				if(i != 0 || j != 0 || k != 0){
					greenFunction[i*NNx*NNy+j*NNx+k] = -1.0/sqrt(pow(i,2)+pow(j,2)+pow(k,2));
				}
			}
		}
	}

	for(int i=Nx+1;i<NNx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				greenFunction[i*NNx*NNy+j*NNx+k] = greenFunction[(NNx-i)*NNx*NNy+j*NNx+k];
			}
		}
	}
	for(int i=0;i<Nx;i++){
		for(int j=Ny+1;j<NNy;j++){
			for(int k=0;k<Nz;k++){
				greenFunction[i*NNx*NNy+j*NNx+k] = greenFunction[i*NNx*NNy+(NNy-j)*NNx+k];
			}
		}
	}

	for(int i=Nx+1;i<NNx;i++){
		for(int j=Ny+1;j<NNy;j++){
			for(int k=0;k<Nz;k++){
				greenFunction[i*NNx*NNy+j*NNx+k] = greenFunction[(NNx-i)*NNx*NNy+(NNy-j)*NNx+k];
			}
			
		}
	}

	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			for(int k=Nz+1;k<NNz;k++){
				greenFunction[i*NNx*NNy+j*NNx+k] = greenFunction[i*NNx*NNy + j*NNx + (NNz-k)];
			}
		}
	}

	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			for(int k=0;k<NNz;k++){
				gf[i*NNx*NNy + j*NNx + k] = complex<double>(0.,0.);
				gf[i*NNx*NNy + j*NNx + k] += greenFunction[i*NNx*NNy + j*NNx + k];
			}
			
		}
	}

	cufftHandle plan;
	cufftDoubleComplex *data;
    cudaMalloc((void**)&data, sizeof(cufftDoubleComplex)*8*N);
    cudaMemcpy(data, gf, sizeof(double)*8*N*2, cudaMemcpyHostToDevice);

	if (cufftPlan3d(&plan,NNx,NNy,NNz, CUFFT_Z2Z) != CUFFT_SUCCESS) {
        printf("CUFFT error: Plan creation failed.\n");
		exit(1);
    }
    if (cufftExecZ2Z(plan, data, data, CUFFT_FORWARD) != CUFFT_SUCCESS) {
		printf("CUFFT error: ExecZ2Z forward failed.\n");
		exit(1);
    }
    cudaMemcpy(fftgf, data, sizeof(double)*8*N*2, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(plan);
    cudaFree(data);
	
	free(gf);
	free(greenFunction);
	
}

void isolatedPotential(struct grid3D *grid, complex<double>* fftgf){
	double *densityPad;
	int N = grid->N;
	int Nx = grid->Nx;
	int Ny = grid->Ny;
	int Nz = grid->Nz;

	int NNx = 2*Nx;
	int NNy = 2*Ny;
	int NNz = 2*Nz;

	densityPad = (double*)malloc(8*N*sizeof(double));
	cufftHandle p1,p2;				//fft plan for cuFFT
	
	complex<double> *dp, *fftdp;	//For Padding density
	
	complex<double> *phi,*ifftphi;	//For potential
	
	dp = (complex<double>*) malloc(sizeof(complex<double>) * 8*N);
	fftdp = (complex<double>*) malloc(sizeof(complex<double>) * 8*N);
	phi = (complex<double>*) malloc(sizeof(complex<double>) * 8*N);
	ifftphi = (complex<double>*) malloc(sizeof(complex<double>) * 8*N);

	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			for(int k=0;k<NNz;k++){
				if(i < Nx && j < Ny && k < Nz){
					//Copy initial density
					densityPad[i*NNx*NNy+j*NNx+k] = grid->density[i*Nx*Ny+j*Nx+k];
				}else{
					//Padding 0s
					densityPad[i*NNx*NNy+j*NNx+k] = 0.0;
				}
			}
		}
	}
	

	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			for(int k=0;k<NNz;k++){
				dp[i*NNx*NNy + j*NNx + k] = complex<double>(0.,0.);
				dp[i*NNx*NNy + j*NNx + k] += densityPad[i*NNx*NNy + j*NNx + k];
			}
			
		}
	}

    cufftDoubleComplex *data;
    cudaMalloc((void**)&data, sizeof(cufftDoubleComplex)*8*N);
    cudaMemcpy(data, dp, sizeof(double)*8*N*2, cudaMemcpyHostToDevice);
	
	//cufft
	if (cufftPlan3d(&p1,NNx,NNy,NNz, CUFFT_Z2Z) != CUFFT_SUCCESS) {
        printf("CUFFT error: Plan creation failed.\n");
		exit(1);
    }
    if (cufftExecZ2Z(p1, data, data, CUFFT_FORWARD) != CUFFT_SUCCESS) {
		printf("CUFFT error: ExecZ2Z forward failed.\n");
		exit(1);
    }
    cudaMemcpy(fftdp, data, sizeof(double)*8*N*2, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(p1);
    cudaFree(data);
	

	for(int i=0;i<NNx;i++){
		for(int j=0;j<NNy;j++){
			for(int k=0;k<NNz;k++){
				//Multiply 2 imaginary numbers
				int index = i*NNx*NNy + j*NNx + k;
				// phi[index][0] = fftdp[index][0] * fftgf[index][0] - fftdp[index][1] * fftgf[index][1];
				// phi[index][1] = fftdp[index][0] * fftgf[index][1] + fftdp[index][1] * fftgf[index][0];
				phi[index] = fftdp[index] * fftgf[index];
			}
		}
	}

	cufftDoubleComplex *data2;
    cudaMalloc((void**)&data2, sizeof(cufftDoubleComplex)*8*N);
    cudaMemcpy(data, phi, sizeof(double)*8*N*2, cudaMemcpyHostToDevice);
	
	////ifft cufft
	if (cufftPlan3d(&p2,NNx,NNy,NNz, CUFFT_Z2Z) != CUFFT_SUCCESS) {
        printf("CUFFT error: Plan creation failed.\n");
		exit(1);
    }
    if (cufftExecZ2Z(p2, data2, data2, CUFFT_INVERSE) != CUFFT_SUCCESS) {
		printf("CUFFT error: ExecZ2Z forward failed.\n");
		exit(1);
    }
    cudaMemcpy(ifftphi, data2, sizeof(double)*8*N*2, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(p2);
    cudaFree(data2);

	
	

	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				int index1 = i*Nx*Ny+j*Nx+k;//index for N size grid.
				int index2 = i*NNx*NNy + j*NNx + k;// index for 2N size grid.
				grid->phi[index1] = -1.0/grid->dx / (8*N) * abs(ifftphi[index2]);
			}
		}		
	}

	// for (int i=0; i < Nx; i++){
	// 	for (int j=0; j < Ny; j++){  
	// 		for(int k=0 ; k < Nz ; k++){
	// 			_2nd_order_diff_3d(grid, i, j,k);
	// 		}   
	// 	}
	// }

	//Use GPU to calculate force from potential.
	_2nd_order_diff_3d_cuda(grid);

	free(densityPad);
	free(dp);
	free(fftdp);
	free(phi);
	free(ifftphi);

}

void Weight(struct grid3D *grid,struct particle3D *particle,int type){

	//Initialize Density field
	for(int i=0 ; i < grid->N ; i++){
		grid->density[i]=0.0;
	}

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
			grid->density[(lx+sx)*grid->Nx*grid->Ny + (ly+sy)*grid->Nx + (lz+sz) ] += particle->mass[i];
		}else if(type == 1){
			//CIC
			for(int zz=0;zz<2;zz++){
				for(int j=0;j<4;j++){
					int p = j / 2;
					int q = j % 2;
					double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
					wFactor *= (1-fabs(particle->z[i]-shift-(lz+zz)*grid->dz)/grid->dz);
					grid->density[(lx+p)*grid->Nx*grid->Ny + (ly+q)*grid->Nx + (lz+zz) ] += particle->mass[i] * wFactor;
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
						int index = indx*grid->Nx*grid->Ny + indy*grid->Nx + indz;
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
			int pos = (lx+sx)*grid->Nx*grid->Ny + (ly+sy)*grid->Nx + (lz+sz) ;
			particle->Fx[i]=grid->Fx[pos]*particle->mass[i];
			particle->Fy[i]=grid->Fy[pos]*particle->mass[i];
			particle->Fz[i]=grid->Fz[pos]*particle->mass[i];
		}else if(type == 1){
			for(int zz=0;zz<2;zz++){
				for(int j=0;j<4;j++){
					int p = j / 2;
					int q = j % 2;
					double wFactor = (1-fabs(particle->x[i]-shift-(lx+p)*grid->dx)/grid->dx)*(1-fabs(particle->y[i]-shift-(ly+q)*grid->dy)/grid->dy);
					wFactor *= (1-fabs(particle->z[i]-shift-(lz+zz)*grid->dz)/grid->dz);
					int pos = (lx+p)*grid->Nx*grid->Ny + (ly+q)*grid->Nx + (lz+zz) ;
					particle->Fx[i] += grid->Fx[pos] * wFactor*particle->mass[i];
					particle->Fy[i] += grid->Fy[pos] * wFactor*particle->mass[i];
					particle->Fz[i] += grid->Fz[pos] * wFactor*particle->mass[i];
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
						int index = indx*grid->Nx*grid->Ny + indy*grid->Nx + indz;
						double weight = weightX[xx+1]*weightY[yy+1]*weightZ[zz+1];
						particle->Fx[i] += grid->Fx[index]*weight*particle->mass[i];
						particle->Fy[i] += grid->Fy[index]*weight*particle->mass[i];
						particle->Fz[i] += grid->Fz[index]*weight*particle->mass[i];
					}
					
				}
			}
		}

		
	}
}
void kick(struct particle3D *particle , double dt){
	double ax,ay,az;
	for(int i=0 ; i<particle->number ; i++){
		//Cauculate the acceleration of each particle.
		ax = particle->Fx[i] / particle->mass[i];
		ay = particle->Fy[i] / particle->mass[i];
		az = particle->Fz[i] / particle->mass[i];

		//Calculate velocity of each particle.
		particle->vx[i] += ax * dt;
		particle->vy[i] += ay * dt;
		particle->vz[i] += az * dt;
	}
}
void drift(struct particle3D *particle , double dt){
	for(int i=0 ; i<particle->number ; i++){
		particle->x[i] += particle->vx[i] * dt;
		particle->y[i] += particle->vy[i] * dt;
		particle->z[i] += particle->vz[i] * dt;
	}
}
void init_rk4(struct particle3D *particle, struct particle3D *buff, struct rk43D *rk4){
	for(int i=0 ; i<particle->number ; i++){
		buff->vx[i] = particle->vx[i];
		buff->vy[i] = particle->vy[i];
		buff->vz[i] = particle->vz[i];
		buff->Fx[i] = particle->Fx[i];
		buff->Fy[i] = particle->Fy[i];
		buff->Fz[i] = particle->Fz[i];
		rk4->ax[i] = 0;
		rk4->ay[i] = 0;
		rk4->az[i] = 0;
		rk4->vx[i] = 0;
		rk4->vy[i] = 0;
		rk4->vz[i] = 0;
		rk4->step1 = 0;
		rk4->step2 = 1;
	}
}
void rk4_mid(struct particle3D *particle, struct particle3D *buff, struct rk43D *rk4, double dt, int weighting){
	for(int i=0 ; i<particle->number ; i++){
		buff->Fx[i] = rk4->step1 * buff->Fx[i] + rk4->step2 * particle->Fx[i];
		buff->Fy[i] = rk4->step1 * buff->Fy[i] + rk4->step2 * particle->Fy[i];
		buff->Fz[i] = rk4->step1 * buff->Fz[i] + rk4->step2 * particle->Fz[i];

		buff->x[i] = particle->x[i] + buff->vx[i] * dt;
		buff->y[i] = particle->y[i] + buff->vy[i] * dt;
		buff->z[i] = particle->z[i] + buff->vz[i] * dt;

		rk4->ax[i] += weighting/6.0 * buff->Fx[i]/particle->mass[i];
		rk4->ay[i] += weighting/6.0 * buff->Fy[i]/particle->mass[i];
		rk4->az[i] += weighting/6.0 * buff->Fz[i]/particle->mass[i];
		rk4->vx[i] += weighting/6.0 * buff->vx[i];
		rk4->vy[i] += weighting/6.0 * buff->vy[i];
		rk4->vz[i] += weighting/6.0 * buff->vz[i];

		buff->vx[i] = particle->vx[i] + buff->Fx[i]/particle->mass[i] * dt;
		buff->vy[i] = particle->vy[i] + buff->Fy[i]/particle->mass[i] * dt;
		buff->vz[i] = particle->vz[i] + buff->Fz[i]/particle->mass[i] * dt;

		rk4->step1 = 1;
		rk4->step2 = 0;
	}
}
void rk4_end(struct particle3D *particle, struct rk43D *rk4, double dt){
	for(int i=0 ; i<particle->number ; i++){
		particle->x[i]  += rk4->vx[i] * dt;
		particle->y[i]  += rk4->vy[i] * dt;
		particle->z[i]  += rk4->vz[i] * dt;
		particle->vx[i] += rk4->ax[i] * dt;
		particle->vy[i] += rk4->ay[i] * dt;
		particle->vz[i] += rk4->az[i] * dt;
	}
}
void periodic_boundary(double position, double length){
	int sign = position/abs(position);
	position = sign * remainder(abs(position + sign*length/2), length) - sign*length/2;
	cout << "A particle reaches the boundary." << endl;
}
void boundary_check(int boundary, struct particle3D *particle, double L){
	if (boundary == 0){
		for (int i=0; i<particle->number; i++){
			if ( abs(particle->x[i]) > L/2){
				periodic_boundary(particle->x[i],L);
			}
			if ( abs(particle->y[i]) > L/2){
				periodic_boundary(particle->y[i],L);
			}
			if ( abs(particle->z[i]) > L/2){
				periodic_boundary(particle->z[i],L);
			}
		}
	}
	else if (boundary == 1){
		for (int i=0; i<particle->number; i++){
			if ( abs(particle->x[i]) > L/2 || abs(particle->y[i]) > L/2 || abs(particle->z[i]) > L/2){
				particle->mass[i] = 0;
				cout << "A particle reaches the boundary." << endl;
			}
		}
	}
}
//Functions to locate memory and free memory of different struct.
void locateMemoryParticle(struct particle3D *particle,int N){
	particle->mass = (double*)malloc(N*sizeof(double));
    particle->x = (double*)malloc(N*sizeof(double));
    particle->y = (double*)malloc(N*sizeof(double));
    particle->z = (double*)malloc(N*sizeof(double));
	particle->Fx = (double*)malloc(N*sizeof(double));
    particle->Fy = (double*)malloc(N*sizeof(double));
    particle->Fz = (double*)malloc(N*sizeof(double));
	particle->vx = (double*)malloc(N*sizeof(double));
    particle->vy = (double*)malloc(N*sizeof(double));
	particle->vz = (double*)malloc(N*sizeof(double));
}
void freeMemoryParticle(struct particle3D *particle){
	free(particle->mass);
	free(particle->x);
	free(particle->y);
	free(particle->z);
	free(particle->Fx);
	free(particle->Fy);
	free(particle->Fz);
	free(particle->vx);
	free(particle->vy);
	free(particle->vz);
}
void locateMemoryRk4(struct rk43D *rk4,int N){
	rk4->ax = (double*)malloc(N*sizeof(double));
	rk4->ay = (double*)malloc(N*sizeof(double));
	rk4->az = (double*)malloc(N*sizeof(double));
	rk4->vx = (double*)malloc(N*sizeof(double));
	rk4->vy = (double*)malloc(N*sizeof(double));
	rk4->vz = (double*)malloc(N*sizeof(double));	
}
void freeMemoryRk4(struct rk43D *rk4){
	free(rk4->ax);
	free(rk4->ay);
	free(rk4->az);
	free(rk4->vx);
	free(rk4->vy);
	free(rk4->vz);
}
void locateMemoryGrid(struct grid3D *grid){
	grid->density = (double*) malloc(grid->N*sizeof(double));	//Density on the grid
	grid->phi = (double*) malloc(grid->N*sizeof(double));		//Potential on the grid
	grid->Fx = (double*) malloc(grid->N*sizeof(double));			//Force on grid comes from potential calculation
	grid->Fy = (double*) malloc(grid->N*sizeof(double));
	grid->Fz = (double*) malloc(grid->N*sizeof(double));
}
void freeMemoryGrid(struct grid3D *grid){
	free(grid->density);
	free(grid->phi);
	free(grid->Fx);
	free(grid->Fy);
	free(grid->Fz);
}
