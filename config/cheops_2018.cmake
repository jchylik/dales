# CHEOPS GCC - v. 2018
# use with: 
# module load gnu/4.8.2 netcdf/4.1.3-gcc-4.8.2 intelmpi/2018
# module load hdf5/1.8.11 szlib
# module load cmake

set(CMAKE_Fortran_COMPILER "gfortran" ) # "ifort") #"gfortran")
set(Fortran_COMPILER_WRAPPER mpif90 ) # mpiifort)

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

set(NETCDF_INCLUDE_DIR "/opt/rrzk/lib/netcdf/4.1.3-gcc-4.8.2/include") # set(NETCDF_INCLUDE_DIR "/opt/rrzk/lib/netcdf/4.1.3/include")
set(NETCDF_LIB_1       "/opt/rrzk/lib/netcdf/4.1.3-gcc-4.8.2/lib/libnetcdff.so") # set(NETCDF_LIB_1       "/opt/rrzk/lib/netcdf/4.1.3/lib/libnetcdff.so")
set(NETCDF_LIB_2       "/opt/rrzk/lib/netcdf/4.1.3-gcc-4.8.2/lib/libnetcdf.so") # set(NETCDF_LIB_2       "/opt/rrzk/lib/netcdf/4.1.3/lib/libnetcdf.so")
# set(NETCDF_LIB_3       "/opt/rrzk/lib/netcdf/netcdf-4.1.3/f90/netcdf.mod")
set(HDF5_LIB_1         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5_hl.so")  # set(HDF5_LIB_1         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5_hl.so")
set(HDF5_LIB_2         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5.so") # set(HDF5_LIB_2         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5.so")
set(SZIP_LIB           "/opt/rrzk/lib/szip/szip-2.1/lib/libsz.so")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2}  ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
# set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${NETCDF_LIB_3} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)  
#



