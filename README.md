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
