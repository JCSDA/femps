program fempoisson_driver

use femps_grid_mod
use femps_operators_mod
use femps_testgrid_mod
use femps_solve_mod, only: preliminary, testpoisson
use femps_inout_mod, only: readgridoprs

!use mpi
!use netcdf

implicit none
type(fempsgrid) :: grid, gridread
type(fempsoprs) :: oprs, oprsread

character(len=2) :: gridtype
integer :: igrid

! Choose grid type
! ----------------
gridtype = 'cs'      !Cubesphere (cs) or icosahedral hexagons (ix)


! Create a grid
! -------------
call grid%setup(gridtype)

if (gridtype == 'cs') then
  call cstestgrid(grid,1,3)
endif

! Build the operators
! -------------------
call oprs%setup(grid)
call oprs%build(grid)


! Perform all the setup
! ---------------------
call preliminary(grid,oprs)
print *,'done preliminary'


! Test the poisson solver
! -----------------------
call testpoisson(grid,oprs)
print *,'done testpoisson'


end program fempoisson_driver
