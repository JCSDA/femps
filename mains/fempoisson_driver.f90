program fempoisson_driver

use fempoisson_mod, only: fempoisson

implicit none
type(fempoisson) :: femp


! Perform all the setup
! ---------------------
call femp%preliminary()
print *,'done preliminary'

! Test the poisson solver
! -----------------------
call femp%testpoisson()
print *,'done testpoisson'

! Clean up
! --------
call femp%delete()

end program fempoisson_driver
