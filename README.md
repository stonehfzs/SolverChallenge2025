# SolverChallenge25 sample

-------------------

## Introduction

There are four sample codes for SolveChallenge25: direct_solver.c , iterative_solver.c , MPI_direct_solver.c and MPI_iterative_solver.c .
The sample of solving a matrix using different methods and the specific output form of the submitted results.

## Structure of code

```
README                        instructions on installation
direct_solver.c               direct solver code sample
MPI_direct_solver.c           MPI direct solver code sample
iterative_solver.c            iterative solver code sample
MPI_iterative_solver.c        MPI iterative solver code sample
mysolver.h                    your solver code
mmio.h                        read matrix file
mmio_highlevel.h              read matrix interface
utlise.h                      check answer and memory function file
utlise_long.h                 check answer and memory function file (long double version)
make.inc                      compiler and compile path
Makefile                      Compile file
```

## Installation

You need to update compiler and compile path in file make.inc and run the command:
> **make**

## run command

After installation, four executable files will be generated as follows:

```
direct_solver                 direct solver sample
MPI_direct_solver             MPI direct solver sample
iterative_solver              iterative solver sample
MPI_iterative_solver          MPI iterative solver sample
```

The executable has five parameters, three of which are optional:

> **./executable_file matrix_name.mtx vector_name.rhs base-0/1[options] out_type[options] sys_type[options]**

For base-0/1 if no base is entered, the default is 1.

```
base-0                        read matrix base-0
base-1                        read matrix base-1
```

A total of three outputs are provided for the output type (default is 0):

```
type:                         direct solver        / iterative solver
0                             end to end time      / end to end time
1                             solver + solve time  / preprocess + solver time
2                             solve time           / solver time
```

The sys_type is type of algebraic system (default is 0):

```
sys_type:
0                             real    system
1                             complex system
```

## Output

If the matrix is solved by the direct method, the output is as follows:

```
matrix name :      matrix_name.mtx
vectorb name :     vector_name.rhs
read function :    base-X
type :             X
sys_type :         X
------------------------------------------
CHECK [option out_type time]        X.XXXXXX ms          check type-X time
Total memory usage:                 XXXX.XXXXX MB        check the memory cost
Check || b - Ax || 2             =  X.XXXXXXe-XX         check || b - Ax || 2 answer
Check || b - Ax || MAX           =  X.XXXXXXe-XX         check || b - Ax || MAX answer
Check || b - Ax || 2 / || b || 2 =  X.XXXXXXe-XX         check || b - Ax || 2 / || b || 2 answer
Check MAX { |b - Ax|_i / |b_i| } =  X.XXXXXXe-XX         check MAX { |b - Ax|_i / |b_i| } answer
LD-Check || b - Ax || 2             =  X.XXXXXXe-XX      check || b - Ax || 2 answer with long double 
LD-Check || b - Ax || MAX           =  X.XXXXXXe-XX      check || b - Ax || MAX answer with long double 
LD-Check || b - Ax || 2 / || b || 2 =  X.XXXXXXe-XX      check || b - Ax || 2 / || b || 2 answer with long double
<!-- LD-Check MAX { |b - Ax|_i / |b_i| } =  X.XXXXXXe-XX      check MAX { |b - Ax|_i / |b_i| } answer with long double -->
```

If the matrix is solved by the iterative method, the output is as follows:

```
matrix name :      matrix_name.mtx
vectorb name :     vector_name.rhs
read function :    base-X
type :             X
sys_type :         X
------------------------------------------
CHECK [option out_type time]        X.XXXXXX ms          check type-X time
Total memory usage:                 XXXX.XXXXX MB        check the memory cost
Check || b - Ax || 2             =  X.XXXXXXe-XX         check || b - Ax || 2 answer
Check || b - Ax || MAX           =  X.XXXXXXe-XX         check || b - Ax || MAX answer
Check || b - Ax || 2 / || b || 2 =  X.XXXXXXe-XX         check || b - Ax || 2 / || b || 2 answer
Check MAX { |b - Ax|_i / |b_i| } =  X.XXXXXXe-XX         check MAX { |b - Ax|_i / |b_i| } answer
LD-Check || b - Ax || 2             =  X.XXXXXXe-XX      check || b - Ax || 2 answer with long double 
LD-Check || b - Ax || MAX           =  X.XXXXXXe-XX      check || b - Ax || MAX answer with long double 
LD-Check || b - Ax || 2 / || b || 2 =  X.XXXXXXe-XX      check || b - Ax || 2 / || b || 2 answer with long double
<!-- LD-Check MAX { |b - Ax|_i / |b_i| } =  X.XXXXXXe-XX      check MAX { |b - Ax|_i / |b_i| } answer with long double -->
```

Samples use types double and long double to calculate the final result.

## Execution of Sample

### run the sample

for direct solver

```
./direct_solver  matrix_name.mtx vector_name.rhs base-0/1[options] out_type[options] sys_type[options]
```

for iterative solver

```
./iterative_solver  matrix_name.mtx vector_name.rhs base-0/1[options] out_type[options] sys_type[options]
```

### run the MPI sample

for MPI direct solver

```
mpirun -np MPI_number ./MPI_direct_solver matrix_name.mtx vector_name.rhs base-0/1[options] out_type[options] sys_type[options]
```

for MPI iterative solver

```
mpirun -np MPI_number ./MPI_iterative_solver matrix_name.mtx vector_name.rhs base-0/1[options] out_type[options] sys_type[options]
```
