#define _USE_MATH_DEFINES 
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <ctime>
#include <sys/time.h>
using namespace std;
// scheme selection
const int scheme	= 3;            // 1 KDK 2 DKD 3 RK4
// thread number
const int NThread 	= 3;

//constants
const double 	G	= 1.0;          // gravitational constant
const double 	dt	= 1.0e-2;       // time interval
const int	N	= 60;		// particle number   

//KDK Kick
void kdk_k(int p0, int p1, const double mass[], double pxyz[][3-0], double vxyz[][3-0]){
	double r12	= pow(pow(pxyz[p0][0]-pxyz[p1][0],2) + pow(pxyz[p0][1]-pxyz[p1][1],2) + pow(pxyz[p0][2]-pxyz[p1][2],2), 0.5);
	double force	= G*mass[p1]/(r12*r12);
	double ax	= force*( pxyz[p1][0] - pxyz[p0][0] )/r12;
	double ay	= force*( pxyz[p1][1] - pxyz[p0][1] )/r12;
	double az	= force*( pxyz[p1][2] - pxyz[p0][2] )/r12;
	
	vxyz[p0][0] 	= vxyz[p0][0] + ax*dt/2;
	vxyz[p0][1]	= vxyz[p0][1] + ay*dt/2;
	vxyz[p0][2]	= vxyz[p0][2] + az*dt/2;
}

//KDK Drift
void kdk_d(int p0, double pxyz[][3-0], double vxyz[][3-0]){
	pxyz[p0][0] += vxyz[p0][0]*dt;
	pxyz[p0][1] += vxyz[p0][1]*dt;
	pxyz[p0][2] += vxyz[p0][2]*dt;
}

//DKD Drift
void dkd_d(int p0, double pxyz[][3-0], double vxyz[][3-0]){
	pxyz[p0][0] += vxyz[p0][0]*dt/2;
	pxyz[p0][1] += vxyz[p0][1]*dt/2;
	pxyz[p0][2] += vxyz[p0][2]*dt/2;
}

//DKD Kick
void dkd_k(int p0, int p1, const double mass[], double pxyz[][3-0], double vxyz[][3-0]){
        double r12      = pow(pow(pxyz[p0][0]-pxyz[p1][0],2) + pow(pxyz[p0][1]-pxyz[p1][1],2) + pow(pxyz[p0][2]-pxyz[p1][2],2), 0.5);
        double force	= G*mass[p1]/(r12*r12);
        double ax       = force*( pxyz[p1][0] - pxyz[p0][0] )/r12;
        double ay       = force*( pxyz[p1][1] - pxyz[p0][1] )/r12;
        double az       = force*( pxyz[p1][2] - pxyz[p0][2] )/r12;

        vxyz[p0][0]     = vxyz[p0][0] + ax*dt;
        vxyz[p0][1]     = vxyz[p0][1] + ay*dt;
        vxyz[p0][2]     = vxyz[p0][2] + az*dt;
}

//RK4
//k1
void rk4_k1(int p0, int p1, const double mass[], double pxyz[][3-0], double vxyz[][3-0], double k1_a[][3-0], double k1_v[][3-0], double buf_p2[][3-0]){
	k1_v[p0][0]	= vxyz[p0][0];
	k1_v[p0][1]	= vxyz[p0][1];
	k1_v[p0][2]	= vxyz[p0][2];
	
	double kk1_r12		= pow(pow(pxyz[p0][0]-pxyz[p1][0],2) + pow(pxyz[p0][1]-pxyz[p1][1],2) + pow(pxyz[p0][2]-pxyz[p1][2],2), 0.5);
	double kk1_force	= G*mass[p1]/pow(kk1_r12,2);
	
	k1_a[p0][0]	= kk1_force * ( pxyz[p1][0] - pxyz[p0][0] )/kk1_r12;
	k1_a[p0][1]     = kk1_force * ( pxyz[p1][1] - pxyz[p0][1] )/kk1_r12;
	k1_a[p0][2]     = kk1_force * ( pxyz[p1][2] - pxyz[p0][2] )/kk1_r12;

	buf_p2[p0][0]	= pxyz[p0][0] + k1_v[p0][0] * dt/2;		//particle coordinates for k2
	buf_p2[p0][1]   = pxyz[p0][1] + k1_v[p0][1] * dt/2;
	buf_p2[p0][2]   = pxyz[p0][2] + k1_v[p0][2] * dt/2;
}

void rk4_k2(int p0, int p1, const double mass[], double pxyz[][3-0], double vxyz[][3-0], double k2_a[][3-0], double k2_v[][3-0], double k1_a[][3-0], double buf_p2[][3-0], double buf_p3[][3-0]){
	k2_v[p0][0]	= vxyz[p0][0] + k1_a[p0][0] * dt/2;
	k2_v[p0][1]     = vxyz[p0][1] + k1_a[p0][1] * dt/2;
	k2_v[p0][2]     = vxyz[p0][2] + k1_a[p0][2] * dt/2;

	double kk2_r12		= pow(pow(buf_p2[p0][0]-buf_p2[p1][0],2) + pow(buf_p2[p0][1]-buf_p2[p1][1],2) + pow(buf_p2[p0][2]-buf_p2[p1][2],2), 0.5);
	double kk2_force	= G*mass[p1]/pow(kk2_r12,2);

	k2_a[p0][0]     = kk2_force * ( buf_p2[p1][0] - buf_p2[p0][0] )/kk2_r12;
	k2_a[p0][1]     = kk2_force * ( buf_p2[p1][1] - buf_p2[p0][1] )/kk2_r12;
	k2_a[p0][2]     = kk2_force * ( buf_p2[p1][2] - buf_p2[p0][2] )/kk2_r12;
	
	buf_p3[p0][0]   = pxyz[p0][0] + k2_v[p0][0] * dt/2;             //particle coordinates for k2
	buf_p3[p0][1]   = pxyz[p0][1] + k2_v[p0][1] * dt/2;
	buf_p3[p0][2]   = pxyz[p0][2] + k2_v[p0][2] * dt/2;
}

void rk4_k3(int p0, int p1, const double mass[], double pxyz[][3-0], double vxyz[][3-0], double k3_a[][3-0], double k3_v[][3-0], double k2_a[][3-0], double buf_p3[][3-0], double buf_p4[][3-0]){
	k3_v[p0][0]     = vxyz[p0][0] + k2_a[p0][0] * dt/2;
	k3_v[p0][1]     = vxyz[p0][1] + k2_a[p0][1] * dt/2;
	k3_v[p0][2]     = vxyz[p0][2] + k2_a[p0][2] * dt/2;

	double kk3_r12         = pow(pow(buf_p3[p0][0]-buf_p3[p1][0],2) + pow(buf_p3[p0][1]-buf_p3[p1][1],2) + pow(buf_p3[p0][2]-buf_p3[p1][2],2), 0.5);
	double kk3_force       = G*mass[p1]/pow(kk3_r12,2);

	k3_a[p0][0]     = kk3_force * ( buf_p3[p1][0] - buf_p3[p0][0] )/kk3_r12;
	k3_a[p0][1]     = kk3_force * ( buf_p3[p1][1] - buf_p3[p0][1] )/kk3_r12;
	k3_a[p0][2]     = kk3_force * ( buf_p3[p1][2] - buf_p3[p0][2] )/kk3_r12;

	buf_p4[p0][0]   = pxyz[p0][0] + k3_v[p0][0] * dt;             //particle coordinates for k2
	buf_p4[p0][1]   = pxyz[p0][1] + k3_v[p0][1] * dt;
	buf_p4[p0][2]   = pxyz[p0][2] + k3_v[p0][2] * dt;
}

void rk4_k4(int p0, int p1, const double mass[], double pxyz[][3-0], double vxyz[][3-0], double k4_a[][3-0], double k4_v[][3-0], double k3_a[][3-0], double buf_p4[][3-0]){
        k4_v[p0][0]     = vxyz[p0][0] + k3_a[p0][0] * dt;
        k4_v[p0][1]     = vxyz[p0][1] + k3_a[p0][1] * dt;
        k4_v[p0][2]     = vxyz[p0][2] + k3_a[p0][2] * dt;

        double kk4_r12         = pow(pow(buf_p4[p0][0]-buf_p4[p1][0],2) + pow(buf_p4[p0][1]-buf_p4[p1][1],2) + pow(buf_p4[p0][2]-buf_p4[p1][2],2), 0.5);
        double kk4_force       = G*mass[p1]/pow(kk4_r12,2);

        k4_a[p0][0]     = kk4_force * ( buf_p4[p1][0] - buf_p4[p0][0] )/kk4_r12;
        k4_a[p0][1]     = kk4_force * ( buf_p4[p1][1] - buf_p4[p0][1] )/kk4_r12;
        k4_a[p0][2]     = kk4_force * ( buf_p4[p1][2] - buf_p4[p0][2] )/kk4_r12;		
}

void rk4_anchor(int p0, double pxyz[][3-0], double vxyz[][3-0], double k4_a[][3-0], double k3_a[][3-0], double k2_a[][3-0], double k1_a[][3-0], double k4_v[][3-0],  double k3_v[][3-0], double k2_v[][3-0], double k1_v[][3-0]){
	vxyz[p0][0] += (k1_a[p0][0] + 2*k2_a[p0][0] + 2*k3_a[p0][0] + k4_a[p0][0]) * dt/6;
	vxyz[p0][1] += (k1_a[p0][1] + 2*k2_a[p0][1] + 2*k3_a[p0][1] + k4_a[p0][1]) * dt/6;
	vxyz[p0][2] += (k1_a[p0][2] + 2*k2_a[p0][2] + 2*k3_a[p0][2] + k4_a[p0][2]) * dt/6;

	pxyz[p0][0] += (k1_v[p0][0] + 2*k2_v[p0][0] + 2*k3_v[p0][0] + k4_v[p0][0]) * dt/6;
	pxyz[p0][1] += (k1_v[p0][1] + 2*k2_v[p0][1] + 2*k3_v[p0][1] + k4_v[p0][1]) * dt/6;
	pxyz[p0][2] += (k1_v[p0][2] + 2*k2_v[p0][2] + 2*k3_v[p0][2] + k4_v[p0][2]) * dt/6;
}

int main (int argc, char *argv[])
{
	// measure performance
        struct timeval start, end;
        gettimeofday(&start, NULL);

	// set the number of threads
        //const int NThread = 2;
        omp_set_num_threads(NThread);
	// initial condition
	const double Mass[N]= {2.0, 1.0, 0.5};
	double Position[N][3-0], Velocity[N][3-0];
	double K1_a[N][3-0], K2_a[N][3-0], K3_a[N][3-0], K4_a[N][3-0];
	double K1_v[N][3-0], K2_v[N][3-0], K3_v[N][3-0], K4_v[N][3-0];	
	double Buf_p1[N][3-0], Buf_p2[N][3-0], Buf_p3[N][3-0], Buf_p4[N][3-0];
	
	for (int i=0; i<N;   i++)
	for (int j=0; j<3-0; j++){
		Position[i][j] = 0.0;
		Velocity[i][j] = 0.0;
	}

	Position[0][1] = -1.250;
	Position[1][1] = 1.0;
	Velocity[0][1] = -0.13;
	Velocity[1][0] = 1.0;
	
	double t = 0;	
	int i, j;
	for (int st=0; st<100000; st++){

		//KDK
		if (scheme == 1){
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					kdk_k(i, j, Mass, Position, Velocity);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					kdk_d(i, Position, Velocity);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					kdk_k(i, j, Mass, Position, Velocity);
				}
			}
		}

		//DKD
		if (scheme == 2){
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					dkd_d(i, Position, Velocity);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					dkd_k(i, j, Mass, Position, Velocity);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					dkd_d(i, Position, Velocity);
				}
			}
		}

		//RK4
		if (scheme == 3){
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					rk4_k1(i, j, Mass, Position, Velocity, K1_a, K1_v, Buf_p2);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					rk4_k2(i, j, Mass, Position, Velocity, K2_a, K2_v, K1_a, Buf_p2, Buf_p3);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					rk4_k3(i, j, Mass, Position, Velocity, K3_a, K3_v, K2_a, Buf_p3, Buf_p4);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++)
			for (int j=0; j<N; j++){
				if (i != j){
					rk4_k4(i, j, Mass, Position, Velocity, K4_a, K4_v, K3_a, Buf_p4);
				}
			}
			# pragma omp parallel for private(i,j)
			for (int i=0; i<N; i++){
				rk4_anchor(i, Position, Velocity, K4_a, K3_a, K2_a, K1_a, K4_v, K3_v, K2_v, K1_v);
			}
		}
		st += dt;
	}
	gettimeofday(&end, NULL);
	double time_taken;
	time_taken = (end.tv_sec - start.tv_sec) * 1e6;
	time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
	cout << "Process took " << fixed << time_taken << " seconds" << '\n';

	return 0;
}

