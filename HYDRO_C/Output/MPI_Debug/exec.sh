#!/bin/bash
srun -n 4 /zbox/user/creinh/HPC_HYDRO/HYDRO_C/Bin/hydro_lmpi -i input_sedov_1000x10.nml
python ../../../tools/multi_output_to_movies.py . 1000 10
