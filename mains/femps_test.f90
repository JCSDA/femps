program fempoisson_driver

use femps_mod, only: class_femps => femps
use mpi
use netcdf

implicit none
type(class_femps) :: femps


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
