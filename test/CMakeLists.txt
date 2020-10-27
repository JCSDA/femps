file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

list( APPEND fv3_grid_data
Data/fv3grid_c0006.nc4
Data/fv3grid_c0012.nc4
Data/fv3grid_c0024.nc4
Data/fv3grid_c0048.nc4
Data/fv3grid_c0096.nc4
)

LINK_FILES( "${fv3_grid_data}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

ecbuild_add_test( TARGET  test_femps_csgrid
                  SOURCES mains/test_femps_csgrid.f90
                  LIBS    femps )
