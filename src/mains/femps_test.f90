program fempoisson_driver

use femps_types_mod, only: fempsgrid, fempspbops
use femps_cstestgrid_mod, only: cdtestgrid
use femps_mod, only: class_femps => femps

!use mpi
!use netcdf

implicit none
type(fempsgrid)  :: grid
type(fempspbops) :: pbobs
character(len=2) :: gridtype
logical :: readgrid, readpbops

! Setup
! -----
gridtype = 'cs'      !Cubesphere (cs) or icosahedral hexagons (ix)
readgrid = .false.   !Read the grid from file
readpbops = .false.  !Read the pre-build operators from file


! Allocate the grid variables
! ---------------------------
if (gridtype == 'cs') then
  call set_fempsgrid_cs(grid)
elseif (gridtype == 'ix') then
  call set_fempsgrid_ix(grid)
else
  print*, "gridtype must be either cs or ix"
  return
endif

call allocate_grid(grid)


! Allocate the fempsgrid type
! ---------------------------
if (.not. readgrid) then
  if (gridtype == 'cs') then
    call cstestgrid(grid)
  elseif (gridtype == 'ix') then
    call ixtestgrid(grid)
  endif
  call writegrid(grid)
else
  call readgrid(grid)
endif

!

! Perform all the setup
! ---------------------
call femps%preliminary()
print *,'done preliminary'

! Test the poisson solver
! -----------------------
call femps%testpoisson()
print *,'done testpoisson'

! Clean up
! --------
call femps%delete()

end program fempoisson_driver
