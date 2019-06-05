
void Weight2d(struct grid2D *grid, struct particle2D *particle, int type);

void WeightForce2d(struct grid2D *grid,struct particle2D *particle,int type);

void Weight3d(struct grid3D *grid,struct particle3D *particle, int type);

void WeightForce3d(struct grid3D *grid,struct particle3D *particle,int type);

void poisson_solver_fft_force_2d(int const dim, struct grid2D *grid);

void _2nd_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

void _4th_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

void _6th_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

void poisson_solver_fft_force_3d(int const dim, struct grid3D *grid);

void _2nd_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk );

void _4th_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk );

void _6th_order_diff_3d(struct grid3D *grid, int const ii, int const jj, int const kk );
