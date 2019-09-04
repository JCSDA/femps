module femps_fv3_mod

use netcdf

use femps_kinds_mod, only: kind_real
use femps_grid_mod, only: fempsgrid
use femps_utils_mod

implicit none

private
public fv3grid_to_ugrid

interface fv3grid_to_ugrid
   module procedure fv3grid_to_ugrid_read
   module procedure fv3grid_to_ugrid_
end interface

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine fv3grid_to_ugrid_read(grid,path)

implicit none

type(fempsgrid),  intent(inout) :: grid
character(len=*), intent(in)    :: path

integer :: iv, fxdim, vxdim
character(len=2055) :: filename

real(kind=kind_real), allocatable, dimension(:,:,:) :: flong
real(kind=kind_real), allocatable, dimension(:,:,:) :: flat
real(kind=kind_real), allocatable, dimension(:,:,:) :: vlong
real(kind=kind_real), allocatable, dimension(:,:,:) :: vlat

do iv = 1,grid%ngrids

  print*, path

  write(filename,"(A,A10,I0.4,A4)") trim(path),'/fv3grid_c',grid%ncube(iv),'.nc4'

  call message('Reading fv3 grid from '//trim(filename), trace)

  call fv3grid_read(filename,fxdim,vxdim,flong,flat,vlong,vlat)

  call fv3grid_to_ugrid_(grid,iv,fxdim,vxdim,flat,flong,vlat,vlong)

  deallocate(flong,flat,vlong,vlat)

enddo

end subroutine fv3grid_to_ugrid_read

! --------------------------------------------------------------------------------------------------

subroutine fv3grid_read(filename,fxdim,vxdim,flong,flat,vlong,vlat)

implicit none

character(len=*),                                    intent(in)    :: filename
integer,                                             intent(out)   :: fxdim
integer,                                             intent(out)   :: vxdim
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: flong
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: flat
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: vlong
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: vlat

integer :: ncid, dimid, varid

! Read the fv3 grid from file
! ---------------------------

call nccheck ( nf90_open(trim(filename), NF90_NOWRITE, ncid), "nf90_open"//trim(filename) )

! Get size of cube
call nccheck ( nf90_inq_dimid(ncid, "fxdim",  dimid), "nf90_inq_dimid fxdim" )
call nccheck ( nf90_inquire_dimension(ncid, dimid, len = fxdim), "nf90_inquire_dimension fxdim" )

call nccheck ( nf90_inq_dimid(ncid, "vxdim",  dimid), "nf90_inq_dimid vxdim" )
call nccheck ( nf90_inquire_dimension(ncid, dimid, len = vxdim), "nf90_inquire_dimension vxdim" )

! Allocate arrays
allocate(flong(fxdim,fxdim,6))
allocate(flat (fxdim,fxdim,6))
allocate(vlong(vxdim,vxdim,6))
allocate(vlat (vxdim,vxdim,6))

! Get flong and flat
call nccheck ( nf90_inq_varid (ncid, 'flons', varid), "nf90_inq_varid flons" )
call nccheck ( nf90_get_var   (ncid, varid, flong), "nf90_get_var flons" )
call nccheck ( nf90_inq_varid (ncid, 'flats', varid), "nf90_inq_varid flats" )
call nccheck ( nf90_get_var   (ncid, varid, flat), "nf90_get_var flats" )

! Get vlong and vlat
call nccheck ( nf90_inq_varid (ncid, 'vlons', varid), "nf90_inq_varid vlons" )
call nccheck ( nf90_get_var   (ncid, varid, vlong), "nf90_get_var vlons" )
call nccheck ( nf90_inq_varid (ncid, 'vlats', varid), "nf90_inq_varid vlats" )
call nccheck ( nf90_get_var   (ncid, varid, vlat), "nf90_get_var vlats" )

call nccheck ( nf90_close   (ncid), "nf90_close" )

end subroutine fv3grid_read

! --------------------------------------------------------------------------------------------------

subroutine fv3grid_to_ugrid_(grid,iv,fxdim,vxdim,flat,flong,vlat,vlong)

implicit none

type(fempsgrid),      intent(inout) :: grid
integer,              intent(in)    :: iv
integer,              intent(in)    :: fxdim
integer,              intent(in)    :: vxdim
real(kind=kind_real), intent(in)    :: flong(fxdim,fxdim,6)
real(kind=kind_real), intent(in)    :: flat (fxdim,fxdim,6)
real(kind=kind_real), intent(in)    :: vlong(vxdim,vxdim,6)
real(kind=kind_real), intent(in)    :: vlat (vxdim,vxdim,6)

integer :: ji,jj,jt,n

! Subroutine for converting the fv3 grid into the
! unstructured grid used by femps

! Grid faces
n = 0
do jt = 1,6
  do jj = 1,fxdim
    do ji = 1,fxdim
      n = n+1
      grid%flong(n,iv) = flong(ji,jj,jt)
      grid%flat(n,iv)  = flat (ji,jj,jt)
    enddo
  enddo
enddo

! Grid vertices
n = 0
do jt = 1,6
  do jj = 1,fxdim
    do ji = 1,fxdim
      n = n+1
      grid%vlong(n,iv) = vlong(ji,jj,jt)
      grid%vlat(n,iv)  = vlat (ji,jj,jt)
    enddo
  enddo
enddo

! Two vertices are neglected in above sweep
grid%vlong(grid%nvert(iv)-1,iv) = vlong(1,fxdim+1,1)
grid%vlat (grid%nvert(iv)-1,iv) = vlat (1,fxdim+1,1)
grid%vlong(grid%nvert(iv)  ,iv) = vlong(fxdim+1,1,2)
grid%vlat (grid%nvert(iv)  ,iv) = vlat (fxdim+1,1,2)

end subroutine fv3grid_to_ugrid_

! --------------------------------------------------------------------------------------------------

end module femps_fv3_mod
