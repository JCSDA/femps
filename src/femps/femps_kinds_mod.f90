! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_kinds_mod

 use, intrinsic :: iso_c_binding
 implicit none

 private
 public kind_real

 integer, parameter :: kind_real=c_double

end module femps_kinds_mod
