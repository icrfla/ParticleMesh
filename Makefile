SIZEOBJ = ParticleMesh.o Weighting2d.o poisson_solver_fft_force_2d.o PotentialToForce2d.o

CC=g++-9

CFLAGS += -lfftw3
# CFLAGS += -fopenmp
CFLAGS += -lgsl

FFTW_PATH += /usr/local/Cellar/fftw/3.3.8_1/lib
FFTW_INCLUDE_PATH += /usr/local/Cellar/fftw/3.3.8_1/include

GSL_PATH += /usr/local/Cellar/gsl/2.5/lib
GSL_INCLUDE_PATH += /usr/local/Cellar/gsl/2.5/include/gsl

ParticleMesh: $(SIZEOBJ)
	$(CC) -o $@ $(SIZEOBJ) $(CFLAGS) -L$(FFTW_PATH) -I$(FFTW_INCLUDE_PATH) -L$(GSL_PATH) -I$(GSL_INCLUDE_PATH)

clean:
	rm -f $(SIZEOBJ) *~ ParticleMesh
