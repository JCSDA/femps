program fempoisson_driver

use femps_kinds_mod
use femps_grid_mod
use femps_operators_mod
use femps_testgrid_mod
use femps_solve_mod

implicit none
type(fempsgrid) :: grid
type(fempsoprs) :: oprs

integer :: n, npass, nprt
real(kind=kind_real), allocatable, dimension(:) :: zeta, psi0, psi
real(kind=kind_real) :: psibar, rms, rms_ref, rms_rel


! Create a grid
! -------------
call grid%setup('cs',cube=96)
call cstestgrid(grid,3,3)

call grid%writegrid('grid.nc4')

! Perform all the setup
! ---------------------
call preliminary(grid,oprs)

call oprs%writeoperators(grid,'operators.nc4')

! Problem setup
! -------------
npass = 10 ! Number of iterations
nprt = 5   ! Number of values to print for testing


! Set up test data
! ----------------
allocate(zeta(grid%nface(grid%ngrids)))
allocate(psi0(grid%nface(grid%ngrids)))
allocate(psi (grid%nface(grid%ngrids)))

do n = 1, grid%nface(grid%ngrids)
  psi0(n) = cos(grid%flat(n,grid%ngrids))*sin(grid%flong(n,grid%ngrids)) ! Large-scale part
enddo
psi0(10) = 10.0_kind_real*psi0(10) ! Plus small-scale part

! Remove global mean (to ensure unambiguous result)
psibar = sum(psi0*grid%farea(:,grid%ngrids))/sum(grid%farea(:,grid%ngrids))
psi0 = psi0 - psibar

! Convert to area integrals
psi0 = psi0*grid%farea(:,grid%ngrids)


! Laplacian
! ---------
call laplace(grid,oprs,grid%ngrids,psi0,zeta)


! Inverse Laplacian
! -----------------
call inverselaplace(grid,oprs,grid%ngrids,npass,zeta,psi)


! Print output
! ------------
print *,' '
do n = 1,nprt
  print *,'Original field psi0     = ', psi0(n)
  print *,'Soln of Poisson eqn psi = ', psi(n)
enddo
psi0 = psi0/grid%farea(:,grid%ngrids)
psi  = psi/grid%farea(:,grid%ngrids)

print *,' '
rms = sqrt(sum((psi0-psi)*(psi0-psi))/grid%nface(grid%ngrids))
print *,'RMS err in global problem = ', rms


! Clean up
! --------
call oprs%delete()
call grid%delete()


! Test pass/fail
! --------------
rms_ref = 5.632313769573670E-004_kind_real
rms_rel = (rms - rms_ref)

if (rms_rel > 1.0e-16_kind_real) then
  print*, ' '
  print*, 'Final RMS', rms
  print*, 'Reference final RMS', rms_ref
  print*, 'Relative difference', rms_rel
  print*, 'Failed with requirement relative difference <= 1.0e-16'
  print*, ' '
  stop 1
else
  print*, 'Test passed'
endif

deallocate(zeta, psi0, psi)

end program fempoisson_driver
