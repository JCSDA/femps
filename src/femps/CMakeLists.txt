list( APPEND femps_src_files
femps_const_mod.f90
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
                     LIBS            ${NETCDF_LIBRARIES}
                     LINKER_LANGUAGE ${FEMPS_LINKER_LANGUAGE}
                   )
