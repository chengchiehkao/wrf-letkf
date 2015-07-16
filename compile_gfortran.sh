#!/bin/bash

cd /SAS002/zerocustom/WRFLETKF/v0/

rm *.o *.mod *.exe

compiler=gfortran44
#optimization_flag="-g -frange-check"
optimization_flag=-O2
openmp_flag=-fopenmp
#netcdf_path=/opt/pgi-libs/  # for REAS
#netcdf_path=/usr/local/netcdf/  # for 245
#netcdf_path=/work/zerocustom/netcdf/netcdf-fortran-4.2/  # for 201
#netcdf_path=/SAS002/zerocustom/netcdf/  # for 201
#netcdf_path=/SAS002/zerocustom/netcdf/201/gcc_gfortran/  # for 201
#netcdf_path=/SAS002/zerocustom/netcdf/REAS/gcc44_gfortran44/  # for REAS
#netcdf_path=/SAS002/zerocustom/20140823/gcc_gfortran_4.4.7/  # for REAS
netcdf_path=/SAS002/zerocustom/netcdf/REAS/buildEnv/netcdf-fortran-4.2/  # for REAS
netcdf_flag="-lnetcdf -lnetcdff"  # for REAS
#netcdf_flag="-lnetcdff"  # for 201 & 245
lapack_path=/SAS002/zerocustom/LAPACK/V3.4.2/build_REAS/gfortran_v5/liblapack.a

#export LD_LIBRARY_PATH=/SAS002/zerocustom/20150713/:/SAS002/zerocustom/20140823/gcc_gfortran_4.4.7/lib/:$LD_LIBRARY_PATH

$compiler  -fdefault-real-8 -ffree-line-length-512 -c mod_basicUtility.f90 ${optimization_flag}
$compiler  ${openmp_flag} ${optimization_flag} -ffree-line-length-512 *.o main_WRFLETKF.f90 ${lapack_path} -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib/ ${netcdf_flag}
#$compiler  ${openmp_flag} ${optimization_flag} -ffree-line-length-512 *.o /${netcdf_path}/lib/libnetcdff.a /${netcdf_path}/lib/libnetcdf.a main_WRFLETKF.f90 ${lapack_path} -o wrfletkf.exe -I${netcdf_path}/include/

