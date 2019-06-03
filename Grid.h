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
