#!/bin/bash

pahtOfScript=$(cd "$(dirname "$0")"; pwd)
echo $pathOfScript

rm *.o *.mod *.exe

compiler=ifort
#optimization_flag="-g -check bounds -limf -debug inline-debug-info"  # for debug
optimization_flag="-O2 -limf"
openmp_flag=-openmp
netcdf_path=/SAS002/zerocustom/netcdf/REAS/icc_ifort_15.0.1.133/
netcdf_flag="-lnetcdf -lnetcdff"  # for REAS
lapack_path="/SAS002/zerocustom/LAPACK/V3.4.2/build_REAS/ifort_v8/liblapack.a /SAS002/zerocustom/BLAS/OpenBLAS/V0.2.14/build_REAS/ifort_v1/lib/libopenblas.a"


$compiler  -r8 -c mod_basicUtility.f90 ${optimization_flag} ${openmp_flag}
$compiler  ${openmp_flag} ${optimization_flag} *.o main_WRFLETKF.f90 ${lapack_path} -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib/ ${netcdf_flag}

