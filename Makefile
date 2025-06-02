include make.inc

all:direct_solver iterative_solver MPI_direct_solver MPI_iterative_solver

LDLIBS += -L /usr/lib/x86_64-linux-gnu/libm.a -lm



direct_solver:direct_solver.c
	$(GCC) $< -o $@ $(CFLAGS)

iterative_solver:iterative_solver.c
	$(GCC) $< -o $@ $(CFLAGS)

MPI_direct_solver:MPI_direct_solver.c
	$(MPICC) $< -o $@ $(MPICCFLAGS) 

MPI_iterative_solver:MPI_iterative_solver.c
	$(MPICC) $< -o $@ $(MPICCFLAGS)

clean:
	rm -f direct_solver iterative_solver MPI_direct_solver MPI_iterative_solver


