#!/bin/bash

pahtOfScript=$(cd "$(dirname "$0")"; pwd)
echo $pathOfScript

rm *.o *.mod *.exe

compiler=gfortran44
#optimization_flag="-g -frange-check"  # for debug
optimization_flag=-O2
openmp_flag=-fopenmp
netcdf_path=/SAS002/zerocustom/netcdf/REAS/gcc44_gfortran44/  # for REAS
netcdf_flag="-lnetcdf -lnetcdff"  # for REAS
lapack_path="/SAS002/zerocustom/LAPACK/V3.4.2/build_REAS/gfortran_v7/liblapack.a /SAS002/zerocustom/BLAS/OpenBLAS/V0.2.14/build_REAS/gfortran_v1/lib/libopenblas.a"


$compiler  -fdefault-real-8 -ffree-line-length-512 -c mod_basicUtility.f90 ${optimization_flag}
$compiler  ${openmp_flag} ${optimization_flag} -ffree-line-length-512 *.o main_WRFLETKF.f90 ${lapack_path} -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib/ ${netcdf_flag}

