

FC = gfortran44 -ffree-line-length-512
preprocess_flag = -cpp -E
optimization_flag = -O2
#optimization_flag = -g -frange-check
debug_flag = -g -frange-check
openmp_flag = -fopenmp
netcdf_path = /SAS002/zerocustom/netcdf/REAS/gcc44_gfortran44/
netcdf_flag = -I$(netcdf_path)/include/ -L$(netcdf_path)/lib/  -lnetcdf -lnetcdff  # for REAS
lapack_path = /SAS002/zerocustom/LAPACK/V3.4.2/build_REAS/gfortran_v7/liblapack.a /SAS002/zerocustom/BLAS/OpenBLAS/V0.2.14/build_REAS/gfortran_v1/lib/libopenblas.a


wrfletkf.exe: mod_derivedtype.o mod_basicutility.o mod_ioutility.o mod_systemutility.o mod_math.o mod_assimilationutility.o main_wrfletkf.o
	$(FC) $(optimization_flag) $(openmp_flag) $(netcdf_flag) *.o $(lapack_path) \
		-o wrfletkf.exe

main_wrfletkf.o: main_WRFLETKF.f90
	$(FC) $(optimization_flag) $(openmp_flag) -c main_WRFLETKF.f90 \
		-o main_wrfletkf.o

mod_assimilationutility.o: mod_derivedtype.o mod_basicutility.o mod_systemutility.o mod_math.o mod_assimilationUtility.f90 assimilationUtility/*.f90
	$(FC) $(optimization_flag) $(openmp_flag) -c mod_assimilationUtility.f90 \
		-o mod_assimilationutility.o

mod_math.o: mod_math.f90 math/*.f90
	$(FC) $(optimization_flag) $(openmp_flag) -c mod_math.f90 \
		-o mod_math.o

mod_systemutility.o: mod_derivedtype.o mod_basicutility.o mod_systemUtility.f90 systemUtility/*.f90 systemUtility/*.F90
	(cd systemUtility/; ls *.F90 | cut -d . -f 1 | xargs -t -I@ sh -c "$(FC) $(preprocess_flag) @.F90 > @.f90" ); \
	$(FC) $(optimization_flag) $(openmp_flag) -c mod_systemUtility.f90  \
		-o mod_systemutility.o

mod_ioutility.o: mod_derivedtype.o mod_IOUtility.f90 IOUtility/*.f90
	$(FC) $(optimization_flag) $(openmp_flag) $(netcdf_flag) -c mod_IOUtility.f90 \
		-o mod_ioutility.o

mod_basicutility.o: mod_basicUtility.f90 basicUtility/*.f90
	$(FC) $(optimization_flag) $(openmp_flag) -fdefault-real-8 -c mod_basicUtility.f90 \
		-o mod_basicutility.o

mod_derivedtype.o: mod_derivedType.f90
	$(FC) $(optimization_flag) $(openmp_flag) -c mod_derivedType.f90 \
		-o mod_derivedtype.o

clean:
	(cd systemUtility/; rm *.s; ls *.F90 | cut -d . -f 1 | xargs -t -I@ sh -c "rm @.f90" ); \
	rm *.o *.mod *.exe

