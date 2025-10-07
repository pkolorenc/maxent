ulimit -s unlimited

# set the intel compiler and MKL environment variables
source /opt/intel/oneapi/compiler/latest/env/vars.sh
source /opt/intel/oneapi/mkl/2022.1.0/env/vars.sh

ifort -V
icc -V
echo $MKLROOT

# set the number of threads available for multi-threaded MKL routines
export OMP_STACKSIZE=1G
export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4
