# (C) Copyright 2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# OpenMP
# ------

if( HAVE_OMP )

  if( CMAKE_Fortran_COMPILER_ID MATCHES "Clang" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -homp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp=omp" )
    set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=omp" )

  else()

    message( STATUS "Fortran compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")

  endif()

else()

  if( CMAKE_Fortran_COMPILER_ID MATCHES "Clang" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-openmp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hnoomp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-openmp")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp-stubs")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  ")

  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp=noomp" )
    set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=noomp" )

  else()

    message( STATUS "Fortran compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")

  endif()

endif()
