SIZEOBJ = ParticleMesh.o Weighting.cpp poisson_solver_fft_2d.o PotentialToForce2d.o

CC=g++

ParticleMesh: $(SIZEOBJ)
	$(CC) -o $@ $(SIZEOBJ) -lgsl -lfftw3

clean:
	rm -f $(SIZEOBJ) *~ ParticleMesh