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
  integer :: npass
  real(kind=kind_real), allocatable, dimension(:) :: zeta, psi0, psi
end type testdata

real(kind=kind_real) :: rms_fem, rms_fv3, rms_ref, rms_rel

type(testdata) :: test_fem
type(testdata) :: test_fv3

! Create the FEMPS test grid
! --------------------------
call grid_fem%setup('cs',ngrids=6,cube=96)
call cstestgrid(grid_fem)


! Create the FV3 test grid
! ------------------------
call grid_fv3%setup('cs',ngrids=5,cube=96)
call fv3grid_to_ugrid(grid_fv3,'/gpfsm/dnb31/drholdaw/JediDev/fv3-bundle/build-intel-17.0.7.259-release-default/fv3-jedi/test')


! Build the connectivity and extra geom
! -------------------------------------
call grid_fem%build_cs(1,1)
call grid_fv3%build_cs(1,1)


! Write grids
! -----------
call grid_fem%writegrid('/gpfsm/dnb31/drholdaw/JediScratch/femps/griddata/grid_fem.nc4')
call grid_fv3%writegrid('/gpfsm/dnb31/drholdaw/JediScratch/femps/griddata/grid_fv3.nc4')


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
call inverselaplace(grid_fem,oprs_fem,grid_fem%ngrids,test_fem%npass,test_fem%zeta,test_fem%psi)
call inverselaplace(grid_fv3,oprs_fv3,grid_fv3%ngrids,test_fv3%npass,test_fv3%zeta,test_fv3%psi)


! Print output
! ------------
call test_endloop(test_fem,grid_fem,rms_fem)
call test_endloop(test_fv3,grid_fv3,rms_fv3)


! Test pass/fail - femps
! --------------
rms_ref = 5.632313769573670E-004_kind_real
rms_rel = (rms_fem - rms_ref)

if (rms_rel > 1.0e-15_kind_real) then
  print*, ' '
  print*, 'Final RMS', rms_fem
  print*, 'Reference final RMS', rms_ref
  print*, 'Relative difference', rms_rel
  print*, 'Failed with requirement relative difference <= 1.0e-16'
  print*, ' '
  stop 1
else
  print*, 'Test passed'
endif

! Test pass/fail - fv3
! --------------
rms_ref = 9.367543386191256E-002
rms_rel = (rms_fv3 - rms_ref)

if (rms_rel > 1.0e-15_kind_real) then
  print*, ' '
  print*, 'Final RMS', rms_fv3
  print*, 'Reference final RMS', rms_ref
  print*, 'Relative difference', rms_rel
  print*, 'Failed with requirement relative difference <= 1.0e-16'
  print*, ' '
  stop 1
else
  print*, 'Test passed'
endif


! Clean up
! --------
call oprs_fem%delete()
call grid_fem%delete()
call oprs_fv3%delete()
call grid_fv3%delete()

call test_delete(test_fem)
call test_delete(test_fv3)

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine test_setup(test,grid)

implicit none

type(testdata),  intent(inout) :: test
type(fempsgrid), intent(in)    :: grid

integer :: n
real(kind=kind_real) :: psibar

test%npass = 10 ! Number of iterations

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
