cmake_minimum_required(VERSION 3.28)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})

#set(ADOLC $ENV{ADOLC_DIR})

#set(LAPACK $ENV{LAPACK_DIR})

#set(HDF5 $ENV{HDF5_DIR})

set(SRC $ENV{ParMeltMantle})

# PKG_CONFIG_PATH is set mannually in bashrc

# Remove the lines below if you do not wish to have PETSc determine the compilers
execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
SET(CMAKE_C_COMPILER ${C_COMPILER})
execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
if (CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
endif (CXX_COMPILER)
execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
if (FORTRAN_COMPILER)
  SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
  enable_language(Fortran)
endif (FORTRAN_COMPILER)

# Set Optimization flags
set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra")
#set(CMAKE_CXX_FLAGS "-g")

project(test)

# Find packages using find_package 
find_package(PkgConfig REQUIRED)
#find_package(MPI REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_search_module(LAPACKE REQUIRED IMPORTED_TARGET lapacke)
#pkg_search_module(HDF5 REQUIRED IMPORTED_TARGET hdf5)
#pkg_search_module(CGNS REQUIRED IMPORTED_TARGET cgns)

# Build library
# Create a library from the files in general folder
add_library(general_lib
        ${SRC}/src/general/integral.cpp 
        ${SRC}/src/general/mesh.c 
        ${SRC}/src/general/util.cpp
)

target_include_directories(general_lib
    PUBLIC ${SRC}/src/general
	 PUBLIC ${PETSC_DIR}/include
)

target_link_libraries(general_lib
    PkgConfig::LAPACKE
    PkgConfig::PETSC
)

# Create library from io folder 
add_library(io_lib
        ${SRC}/src/io/input.cpp
        ${SRC}/src/io/output.c
#		  ${SRC}/io/cgns_io.c
)

target_include_directories(io_lib
		  PUBLIC ${SRC}/src/io
		  PUBLIC ${PETSC_DIR}/include
)

target_link_libraries(io_lib
    general_lib
    PkgConfig::PETSC
	 #    PkgConfig::HDF5
#	 PkgConfig::CGNS
)

#Create a library from phase folder
add_library(phase_lib
		  ${SRC}/src/phase/eutectic_rescaled.cpp
)

target_include_directories(phase_lib
		  PUBLIC ${SRC}/src/phase
		  PUBLIC ${SRC}/src/include
		  )

#Create a library from files in mfem file
add_library(mfem_lib
        ${SRC}/src/mfem/basis.cpp
        ${SRC}/src/mfem/Hdivmixed.cpp
        ${SRC}/src/mfem/brmixed.cpp
		  ${SRC}/src/mfem/clean/serial_util.cpp
		  ${SRC}/src/mfem/clean/serial_bndry.cpp
		  ${SRC}/src/mfem/clean/serial_locmat.cpp
		  ${SRC}/src/mfem/clean/serial_assemble.cpp
)

target_include_directories(mfem_lib
		  PUBLIC ${SRC}/src/mfem
		  PUBLIC ${SRC}/src/mfem/clean
)

target_link_libraries(mfem_lib
    phase_lib
    general_lib
)

# ===================================================================
# Create a library from the files in the mlweno directory
add_library(mlweno_lib 
		  ${SRC}/src/mlweno_array/tensorstencilpoly.cpp
		  ${SRC}/src/mlweno_array/reconstruction.cpp
		  ${SRC}/src/mlweno_array/error.cpp
)

# Specify any additional include directories if needed
target_include_directories(mlweno_lib 
		  PUBLIC ${SRC}/src/mlweno_array
		  PUBLIC ${SRC}/src/general
		  #PUBLIC ${PETSC_DIR}/include 
)

# Link the PETSc libraries to the mlweno_lib library
target_link_libraries(mlweno_lib 
        PkgConfig::PETSC
)

# ===================================================================
add_library(transport_lib
		  ${SRC}/src/transport/flux.cpp
		  ${SRC}/src/transport/transport.cpp
		  ${SRC}/src/transport/print.cpp
)

# Specify any additional include directories if needed
target_include_directories(transport_lib 
		  PUBLIC ${SRC}/src/transport
		  PUBLIC ${SRC}/src/mlweno_array
		  PUBLIC ${SRC}/src/general
    PUBLIC ${PETSC_DIR}/include 
)

target_link_libraries(transport_lib
        io_lib
        general_lib
        mlweno_lib
        PkgConfig::PETSC
)

# ===================================================================
add_library(couple_lib
		  ${SRC}/driver/couple.cpp
		  ${SRC}/driver/print.cpp
		  porosity.cpp
		  flow.cpp
		  trans.cpp
		  )

target_include_directories(couple_lib
		  PUBLIC ${SRC}/src/phase
		  PUBLIC ${SRC}/src/general
		  PUBLIC ${SRC}/src/mfem
		  PUBLIC ${SRC}/src/mfem/clean
		  PUBLIC ${SRC}/src/mlweno
		  PUBLIC ${SRC}/src/transport
		  PUBLIC ${SRC}/driver
)

target_link_libraries(couple_lib
		  io_lib
		  general_lib
		  phase_lib
		  mfem_lib
		  transport_lib
		  PkgConfig::PETSC
)

add_executable(driver
		  ${SRC}/driver/driver.cpp
)

target_include_directories(driver
		  PUBLIC ${SRC}/driver
		  )

target_link_libraries(driver
		  io_lib
		  phase_lib
		  mfem_lib
		  couple_lib
		  mlweno_lib
		  transport_lib
		  PkgConfig::LAPACKE
		  PkgConfig::PETSC
		  phase_lib
)
