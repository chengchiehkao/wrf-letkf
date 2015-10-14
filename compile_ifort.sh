#!/bin/bash

cd /SAS002/zerocustom/WRFLETKF/v0/

rm *.o *.mod *.exe

compiler=ifort
#optimization_flag="-g -check bounds -limf -debug inline-debug-info"
optimization_flag="-O2 -limf"
openmp_flag=-openmp
#netcdf_path=/SAS002/zerocustom/netcdf/REAS/gcc_ifort/
netcdf_path=/SAS002/zerocustom/netcdf/REAS/icc_ifort_15.0.1.133/
#netcdf_path=/SAS002/zerocustom/20140823/gcc_ifort_4.1.2_afterGCC44Installed/
netcdf_flag="-lnetcdf -lnetcdff"  # for REAS
lapack_path=/SAS002/zerocustom/LAPACK/V3.4.2/build_REAS/ifort_v1/liblapack.a

$compiler  -r8 -c mod_basicUtility.f90 ${optimization_flag} ${openmp_flag}
#$compiler  -r8 -c mod_basicUtility.f90 ${optimization_flag} -recursive #${openmp_flag}
#$compiler  -r8 -c mod_basicUtility.f90 ${optimization_flag}
#$compiler  ${openmp_flag} ${optimization_flag} *.o main_WRFLETKF.f90 -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib/ ${netcdf_flag}
$compiler  ${openmp_flag} ${optimization_flag} *.o main_WRFLETKF.f90 ${lapack_path} -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib/ ${netcdf_flag}

