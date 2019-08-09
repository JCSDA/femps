program fempoisson_driver

use femps_kinds_mod
use femps_grid_mod
use femps_operators_mod
use femps_testgrid_mod
use femps_solve_mod

implicit none
type(fempsgrid) :: grid
type(fempsoprs) :: oprs

integer :: npass = 10 ! Number of passes
integer :: nf, ne, nv, if1, ipass, nprt
real(kind=kind_real), allocatable :: psi0(:), zeta(:), psi(:), ff1(:), ff2(:), ff3(:), ff4(:), &
                                     temp1(:), temp2(:)
real(kind=kind_real) :: long, lat, psibar
real(kind=kind_real) :: rms, rms_ref, rms_rel


! Create a grid
! -------------
call grid%setup('cs')
call cstestgrid(grid,1,3)


! Perform all the setup
! ---------------------
call preliminary(grid,oprs)
print *,'done preliminary'


! Run a test comparing the Laplace versus inverse Laplace
! -------------------------------------------------------

nf = grid%nface(grid%ngrids)
ne = grid%nedge(grid%ngrids)
nv = grid%nvert(grid%ngrids)

print *,' '
print *,'--------------------------'
print *,' '
print *,'Testing mgsolve '
print *,' '

! Number of values to print for testing
nprt = 5

allocate(psi0(nf),zeta(nf),psi(nf),ff1(nf),ff2(nf),ff3(nf),ff4(nf),temp1(ne),temp2(ne))

! Set up test data
! Large-scale part
do if1 = 1, nf
  long = grid%flong(if1,grid%ngrids)
  lat = grid%flat(if1,grid%ngrids)
  ! psi0(if1) = SIN(lat)
  psi0(if1) = COS(lat)*SIN(long)
enddo
! Plus small-scale part
psi0(10) = 10.0_kind_real*psi0(10)
! Remove global mean (to ensure unambiguous result)
psibar = SUM(psi0*grid%farea(:,grid%ngrids))/SUM(grid%farea(:,grid%ngrids))
psi0 = psi0 - psibar
! Convert to area integrals
ff1 = psi0*grid%farea(:,grid%ngrids)
print *,'Original field ff1 =     ',ff1(1:nprt)
print *,' '

! Calculate laplacian
call laplace(grid,oprs,ff1,zeta,grid%ngrids,nf,ne)


! Now invert laplacian to check we get back to where we started

! Initialize result to zero
! Note psi will be stream function times grid cell area
psi = 0.0_kind_real
temp2 = 0.0_kind_real

! Iterate several passes
do ipass = 1, npass

  print *,'Pass ',ipass

  ! Compute residual based on latest estimate
  call massL(oprs,psi,ff2,grid%ngrids,nf)

  call Ddual1(grid,ff2,temp1,grid%ngrids,nf,ne)

  ! Improve the estimate temp2 that we obtained in the previous pass
  call massMinv(grid,oprs,temp1,temp2,grid%ngrids,ne,4)

  call Dprimal2(grid,oprs,temp2,ff3,grid%ngrids,ne,nf)

  ! Residual
  ff4 = zeta - ff3

  ! Now solve the approximate Poisson equation
  ! D2 xminv D1bar L psi' = residual
  call mgsolve(grid,oprs,ff3,ff4,grid%ngrids)

  ! And increment best estimate
  ! *** We could think about adding beta*ff3 here and tuning beta for optimal convergence ***
  psi = psi + ff3

  ! Remove global mean (to ensure unambiguous result)
  psibar = SUM(psi)/SUM(grid%farea(:,grid%ngrids))
  psi = psi - psibar*grid%farea(:,grid%ngrids)

  print *,'Original field ff1      = ',ff1(1:nprt)
  print *,'Soln of Poisson eqn psi = ', psi(1:nprt)
  ff4 = (ff1-psi)/grid%farea(:,grid%ngrids)
  rms = sqrt(sum(ff4*ff4)/nf)
  print *,'RMS err in global problem = ', rms
  print *,' '

enddo

! Clean up
! --------
call oprs%delete()
call grid%delete()


rms_ref = 5.632313769573670E-004_kind_real
rms_rel = (rms - rms_ref)

print *,'done testpoisson'

if (rms_rel > 1.0e-16_kind_real) then
  print*, ' '
  print*, 'Final RMS', rms
  print*, 'Reference final RMS', rms_ref
  print*, 'Relative difference', rms_rel
  print*, 'Failed with requirement relative difference <= 1.0e-16'
  print*, ' '
  stop 1
endif


end program fempoisson_driver
