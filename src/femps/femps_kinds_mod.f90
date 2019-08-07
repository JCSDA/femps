module femps_kinds_mod

 use, intrinsic :: iso_c_binding
 implicit none

 private
 public kind_real

 integer, parameter :: kind_real=c_double

end module femps_kinds_mod
