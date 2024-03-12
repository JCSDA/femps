################################################################################
# FEMPS
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

find_package( ecbuild 3.6 REQUIRED )

project( femps VERSION 1.2.0 LANGUAGES Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECTORY_LABELS "femps" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

set( FEMPS_LINKER_LANGUAGE Fortran )

################################################################################
# Standard compiler flags
################################################################################

include( femps_extra_macros )

################################################################################
# Dependencies
################################################################################

# # MPI
# ecbuild_find_mpi( COMPONENTS Fortran REQUIRED )
# ecbuild_include_mpi()
# link_libraries(${MPI_Fortran_LIBRARIES})

# NetCDF
find_package( jedicmake QUIET )
find_package( NetCDF REQUIRED COMPONENTS Fortran )

################################################################################
# Export package info
################################################################################

list( APPEND FEMPS_TPLS NetCDF)

set( FEMPS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
                        ${CMAKE_CURRENT_SOURCE_DIR}/src
                        ${CMAKE_Fortran_MODULE_DIRECTORY})

set( FEMPS_LIBRARIES femps )

get_directory_property( FEMPS_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${FEMPS_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND FEMPS_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND FEMPS_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND FEMPS_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

include_directories( ${FEMPS_INCLUDE_DIRS}
                     ${FEMPS_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )
add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME femps )

# print the summary of the configuration
ecbuild_print_summary()
