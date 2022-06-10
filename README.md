# fenkar_program
OpenMP parallelized 3D Fenton-Karma finite difference code. Not well optimized, but functions.

To compile:
`gfortran -Ofast -fopenmp -fopt-info-vec -fopenacc fenkar_prog.f90`
To run:
 `export OMP_NUM_THREADS=48; ./a.out`
To plot:
`julia -O3 readData_makie.jl`
