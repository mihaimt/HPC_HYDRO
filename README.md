# HPC HYDRO
This is a work repository for HPC1b.

Contributors: Mihai Tomozeiu, Christian Reinhardt, Rafael Kung.


## how to compile & run

### on your local university desktop opensuse

```
cd HYDRO_C/Src
make lmpi

mpirun -n 4 ../Bin/hydro_lmpi -i ../Input/input_sedov_100x10.nml
```

The binary is in `HYDRO_C/Bin`. You can use `mpirun` to run the program, or use `HYDRO_C/Bin/run`, which runs the command above.


### on a cluster

We store the outputs from the simulations in ../Output/.

How to run a simulation:

There are in principe two possibilities:

1) Allocate one (or more) nodes with salloc
2) Submit a slurm batch job with sbatch

ZBOX:
1) Allocate with salloc --partition=zbox and run with srun -n cores hydro_lmpi
2) sbatch example-zbox.job (see example-zbox.job in this directory for details)

Important: Always compile with make lmpi and run hydro_lmpi!

Dora:
1) Allocate with salloc (no partition) and run with aprun -n cores hydro_mpi
2) sbatch example-dora.job (see example-dora.job in this directory for details and be sure to change ACCOUNT to uzh8)

Important: Always compile with make mpi and run hydro_mpi!




## Generating plots and movies with multi_output_to_movies.py

###proper call:

python multi_output_to_movies.py (path to output directory) (number of steps) (delay for movie frames)

###example:
```
python multi_output_to_movies.py ~/git/HPC_Hydro/HYDRO_C/Output/Input_sedov_100x10 10 100
```
.png and .mpg files will be created in the same directory as the analysed output files


