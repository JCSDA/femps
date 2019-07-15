module femps_constants_mod

use femps_kinds_mod

implicit none
private

! Pi
real(kind=kind_real), parameter, public :: pi = 3.14159265358979323_kind_real
real(kind=kind_real), parameter, public :: piby2 = 0.5_kind_real*pi

! Earth's radius
real(kind=kind_real), parameter, public :: rearth = 6371220.0_kind_real

end module femps_constants_mod
