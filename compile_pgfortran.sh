#!/bin/bash

pahtOfScript=$(cd "$(dirname "$0")"; pwd)
echo $pathOfScript

rm *.o *.mod *.exe

compiler=pgf90
#optimization_flag="-g -Mchkstk -Mchkfpstk -C -Mdwarf2 -Mpgicoff"  # for debug
optimization_flag="-fastsse -Mrecursive -Mreentrant"
openmp_flag="-mp=bind"
netcdf_path=/opt/pgi-libs/  # for REAS
netcdf_flag="-lnetcdf -lnetcdff"  # for REAS
lapack_path=/SAS002/zerocustom/LAPACK/V3.4.2/build_REAS/pgf_v9/liblapack.a


$compiler  -r8 -c mod_basicUtility.f90 ${optimization_flag} ${openmp_flag}
$compiler   ${openmp_flag} ${optimization_flag} *.o main_WRFLETKF.f90 ${lapack_path} -o wrfletkf.exe -I${netcdf_path}/include/ -L${netcdf_path}/lib ${netcdf_flag}

