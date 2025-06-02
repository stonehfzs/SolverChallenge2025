## direct_solver.c

Solve the system directly. (LU Factorization)

## iterative_solver.c

Solve the system iteratively.

## MPI_direct_solver.c / MPI_iterative_solver.c

Using MPI to accelerate the process.

## mmio.h

Read and write the matrix market standard format matrices.

## mmio_highlevel.h

Depend on mmio.h, directly change the mmio.h into CSR format.
`read_mtx_header` could be used to judge the format of the matrices.


## Remark

The actual and detailed process is defined in the file `mysolver.h`. So you should finish that.

Makefile will compile all needed files to make sure all structure are in the whole process.