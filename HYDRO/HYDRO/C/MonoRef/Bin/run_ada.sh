#!/bin/ksh
set -x

# If no argument, use input_sedov_noio_10000x10000_4it.nml as input file
if [ $# != 1 ]
then
INPUT=Others/input_sedov_noio_10000x10000_4it.nml
else
INPUT=$1
fi
DIM=${INPUT%.nml};DIM=${DIM##*[a-z,A-Z]_}

# Execution hydro with $INPUT as input file
cd ../Output
rm output_*
MYDATE=$(date +"%Y_%m_%d:%H_%M_%S")
../Src/hydro -i ../Input/$INPUT | tee ../Bin/hydro_monoref_${DIM}_$MYDATE.res
