list( APPEND femps_src_files
femps_const_mod.f90
femps_fv3_mod.f90
femps_grid_mod.f90
femps_kinds_mod.f90
femps_operators_mod.f90
femps_solve_mod.f90
femps_testgrid_mod.f90
femps_utils_mod.f90
)

#Append includes and libraries with particular model
list (APPEND FEMPS_LIB_DEP ${NETCDF_LIBRARIES})

#Append src with particular model
list (APPEND FEMPS_SRC_DEP ${femps_src_files} )

#Add the library
ecbuild_add_library( TARGET          femps
                     SOURCES         ${femps_src_files}
                     PUBLIC_LIBS     ${NetCDF_LIBRARIES}
                     PUBLIC_INCLUDES ${NetCDF_INCLUDE_DIRS}
                     LINKER_LANGUAGE ${FEMPS_LINKER_LANGUAGE}
                   )

target_include_directories( femps INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
                                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Fortran module output directory for build and install interfaces
set( MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION} )
set_target_properties( ${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR} )
install( DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR} )
target_include_directories( femps INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}> )
