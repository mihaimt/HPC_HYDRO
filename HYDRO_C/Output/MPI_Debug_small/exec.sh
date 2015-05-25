#!/bin/bash
srun -n 2 /zbox/data/creinh/HPC_HYDRO/HYDRO_C/Bin/hydro_lmpi -i input_sedov_50x10.nml
python ../../../tools/multi_output_to_movies.py . 1000 10
