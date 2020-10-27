! (C) Copyright 2019 UCAR and 2011-2018 John Thuburn, University of Exeter, UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

program fempoisson_driver

use femps_kinds_mod
use femps_grid_mod
use femps_operators_mod
use femps_testgrid_mod
use femps_solve_mod
use femps_fv3_mod

implicit none
type(fempsgrid) :: grid_fem, grid_fv3
type(fempsoprs) :: oprs_fem, oprs_fv3

type testdata
  real(kind=kind_real), allocatable, dimension(:) :: zeta, psi0, psi
end type testdata

real(kind=kind_real) :: rms_fem, rms_fv3, rms_ref, rms_rel, tol

type(testdata) :: test_fem
type(testdata) :: test_fv3

! Create the FEMPS test grid
! --------------------------
call grid_fem%setup('cs',ngrids=6,cube=96,niter=10)
call cstestgrid(grid_fem)


! Create the FV3 test grid
! ------------------------
call grid_fv3%setup('cs',ngrids=5,cube=96,niter=50)
call fv3grid_to_ugrid(grid_fv3,'Data/')


! Build the connectivity and extra geom
! -------------------------------------
call grid_fem%build_cs(1,1)
call grid_fv3%build_cs(1,1)


! Write grids
! -----------
call grid_fem%writegrid('grid_fem.nc4')
call grid_fv3%writegrid('grid_fv3.nc4')


! Perform all the setup
! ---------------------
call preliminary(grid_fem,oprs_fem)
call preliminary(grid_fv3,oprs_fv3)

call oprs_fem%pdelete() !Partial deletes
call oprs_fv3%pdelete() !Partial deletes

call oprs_fem%writeoperators(grid_fem,'operators_fem.nc4')
call oprs_fv3%writeoperators(grid_fv3,'operators_fv3.nc4')


! Set up test data
! ----------------
call test_setup(test_fem,grid_fem)
call test_setup(test_fv3,grid_fv3)


! Laplacian
! ---------
call laplace(grid_fem,oprs_fem,grid_fem%ngrids,test_fem%psi0,test_fem%zeta)
call laplace(grid_fv3,oprs_fv3,grid_fv3%ngrids,test_fv3%psi0,test_fv3%zeta)

! Inverse Laplacian
! -----------------
call inverselaplace(grid_fem,oprs_fem,grid_fem%ngrids,test_fem%zeta,test_fem%psi,.true.)
call inverselaplace(grid_fv3,oprs_fv3,grid_fv3%ngrids,test_fv3%zeta,test_fv3%psi,.true.)


! Print output
! ------------
call test_endloop(test_fem,grid_fem,rms_fem)
call test_endloop(test_fv3,grid_fv3,rms_fv3)


! Test pass/fail - femps
! --------------
tol = 1.0e-11_kind_real

rms_ref = 5.632313769573670E-004_kind_real
rms_rel = abs((rms_fem - rms_ref)/rms_ref)


print*, ' '
print*, ' '
print*, 'Result for femps internal grid'
print*, 'Expected RMS', rms_ref
print*, 'Actual RMS', rms_fem
print*, 'Relative difference', rms_rel

if (abs(rms_rel) > tol) then
  print*, ' '
  print*, 'Test failed with requirement that relative difference to expected RMS <=',tol
  print*, ' '
  stop 1
else
  print*, ' '
  print*, 'Test passed'
  print*, ' '
endif

! Test pass/fail - fv3
! --------------
rms_ref = 6.481994273258480E-005_kind_real
rms_rel = abs((rms_fv3 - rms_ref)/rms_ref)

print*, 'Result for fv3 grid'
print*, 'Expected RMS', rms_ref
print*, 'Actual RMS', rms_fv3
print*, 'Relative difference', rms_rel

if (abs(rms_rel) > tol) then
  print*, ' '
  print*, 'Test failed with requirement that relative difference to expected RMS <=',tol
  print*, ' '
  stop 1
else
  print*, ' '
  print*, 'Test passed'
  print*, ' '
endif

! Clean up
! --------
call oprs_fem%delete()
call grid_fem%delete()
call oprs_fv3%delete()
call grid_fv3%delete()

call test_delete(test_fem)
call test_delete(test_fv3)

stop 0

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine test_setup(test,grid)

implicit none

type(testdata),  intent(inout) :: test
type(fempsgrid), intent(in)    :: grid

integer :: n
real(kind=kind_real) :: psibar

allocate(test%zeta(grid%nface(grid%ngrids)))
allocate(test%psi0(grid%nface(grid%ngrids)))
allocate(test%psi (grid%nface(grid%ngrids)))

do n = 1, grid%nface(grid%ngrids)
 test%psi0(n) = cos(grid%flat(n,grid%ngrids))*sin(grid%flong(n,grid%ngrids)) ! Large-scale part
enddo
test%psi0(10) = 10.0_kind_real*test%psi0(10) ! Plus small-scale part

! Remove global mean (to ensure unambiguous result)
psibar = sum(test%psi0*grid%farea(:,grid%ngrids))/sum(grid%farea(:,grid%ngrids))
test%psi0 = test%psi0 - psibar

! Convert to area integrals
test%psi0 = test%psi0*grid%farea(:,grid%ngrids)

end subroutine test_setup

! --------------------------------------------------------------------------------------------------

subroutine test_delete(test)

implicit none

type(testdata), intent(inout) :: test

deallocate(test%zeta, test%psi0, test%psi)

end subroutine test_delete

! --------------------------------------------------------------------------------------------------

subroutine test_endloop(test,grid,rms)

implicit none

type(testdata),       intent(inout) :: test
type(fempsgrid),      intent(in)    :: grid
real(kind=kind_real), intent(out)   :: rms

integer :: n, nprt = 5

print *,' '
do n = 1,nprt
 print *,'Original field psi0     = ', test%psi0(n)
 print *,'Soln of Poisson eqn psi = ', test%psi(n)
enddo
test%psi0 = test%psi0/grid%farea(:,grid%ngrids)
test%psi  = test%psi/grid%farea(:,grid%ngrids)
print *,' '

rms = sqrt(sum((test%psi0-test%psi)*(test%psi0-test%psi))/grid%nface(grid%ngrids))
print *,'RMS err in global problem = ', rms

end subroutine test_endloop

! --------------------------------------------------------------------------------------------------

end program fempoisson_driver
