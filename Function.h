
void Weight(struct grid2D *grid, struct particle2D *particle,int type);

void WeightForce(struct grid2D *grid,struct particle2D *particle,int type);

void poisson_solver_fft_force_2d(int const dim, struct grid2D *grid);

void _2nd_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

void _4th_order_diff_2d(struct grid2D *grid, int const ii, int const jj );

void _6th_order_diff_2d(struct grid2D *grid, int const ii, int const jj );
