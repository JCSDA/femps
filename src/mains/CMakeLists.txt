ecbuild_add_executable( TARGET  femps_test.x
                        SOURCES femps_test.f90
                        LIBS	${NETCDF_LIBRARIES} femps
                      )
