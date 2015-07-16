#!/bin/bash

cd /SAS002/zerocustom/WRFLETKF/v0/

rm *.o *.mod *.exe

compiler=pgf90
#optimization_flag="-g -Mchkstk -Mchkfpstk -C -Mdwarf2 -Mpgicoff"
optimization_flag="-fastsse -Mrecursive -Mreentrant"
#optimization_flag="-fastsse -Mipa=fast,inline"
openmp_flag="-mp=bind"
netcdf_path=/opt/pgi-libs/  # for REAS
#netcdf_path=/SAS002/zerocustom/netcdf/REAS/pgcc_pgfortran/
#netcdf_path=/SAS002/zerocustom/netcdf/201/gcc_pgf90/
#netcdf_path=/usr/local/netcdf/  # for 245
#netcdf_path=/usr/local/netcdf3-pgi/  # for 201
netcdf_flag="-lnetcdf -lnetcdff"  # for REAS
#netcdf_flag="-lnetcdf"  # for 201 & 245
lapack_path=/SAS002/zerocustom/LAPACK/V3.4.2/build_REAS/pgf_v9/liblapack.a

$compiler  -r8 -c mod_basicUtility.f90 ${optimization_flag} ${openmp_flag}
#$compiler  -r8 -c mod_basicUtility.f90 ${optimization_flag} -Mrecursive #${openmp_flag}
#$compiler   ${openmp_flag} ${optimization_flag} *.o main_WRFLETKF.f90 liblapack.a -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib ${netcdf_flag}
$compiler   ${openmp_flag} ${optimization_flag} *.o main_WRFLETKF.f90 ${lapack_path} -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib ${netcdf_flag}
#$compiler   ${openmp_flag} ${optimization_flag} *.o main_WRFLETKF.f90 -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib ${netcdf_flag} -llapack -lblas

