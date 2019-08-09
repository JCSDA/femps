module femps_grid_mod

use netcdf

use femps_kinds_mod

implicit none
private
public fempsgrid

! Type to hold all the grid information
! -------------------------------------
type fempsgrid

  character(len=2) :: gtype  ! cs (cubed-sphere), ih (icosahedral hexagons), di (diamonds)

  integer :: ngrids      ! Number of grids in multigrid hierarchy

  integer :: nfacex      ! Number of faces on each grid
  integer :: nedgex      ! Number of edges on each grid
  integer :: nvertx      ! Number of vertices on each edge

  integer, allocatable, dimension(:) :: nface ! faces on each grid
  integer, allocatable, dimension(:) :: nedge ! edges on each grid
  integer, allocatable, dimension(:) :: nvert ! vertices on each edge

  integer :: n0

  integer :: nefmx, nevmx                       ! maximum in neoff and neofv
  integer, allocatable, dimension(:,:) :: neoff ! number of edges and vertices of each face on each grid
  integer, allocatable, dimension(:,:) :: neofv ! number of edges of each vertex on each grid

  ! Dimensions for the following arrays
  ! -----------------------------------
  integer, private :: dimfnxtf, dimeoff, dimvoff, dimfnxte, dimvofe, dimfofv, dimeofv

  ! Connectivity
  ! ------------
  integer, allocatable, dimension(:,:,:) :: fnxtf   ! faces next to each face on each grid
  integer, allocatable, dimension(:,:,:) :: eoff    ! edges of each face on each grid
  integer, allocatable, dimension(:,:,:) :: voff    ! vertices of each face on each grid
  integer, allocatable, dimension(:,:,:) :: fnxte   ! faces either side of each edge on each grid
  integer, allocatable, dimension(:,:,:) :: vofe    ! vertices at the ends of each edge on each grid
  integer, allocatable, dimension(:,:,:) :: fofv    ! faces around each vertex on each grid
  integer, allocatable, dimension(:,:,:) :: eofv    ! edges incident on each vertex on each grid

  ! Coordinates and geometrical information
  ! ---------------------------------------
  real(kind=kind_real), allocatable, dimension(:,:) :: flong ! longitude of faces on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: flat  ! latitude of faces on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: vlong ! longitude of vertices on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: vlat  ! latitude of vertices on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: farea ! area of faces on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: ldist ! primal edge length, i.e. distance between neighbouring vertices
  real(kind=kind_real), allocatable, dimension(:,:) :: ddist ! dual edge length, i.e. distance between neighbouring face centres

  contains

   procedure, public :: setup
   procedure, public :: delete
   procedure, public :: writegrid

end type fempsgrid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setup(self,gridtype,cube)

implicit none
class(fempsgrid),  intent(inout) :: self
character(len=2),  intent(in)    :: gridtype
integer, optional, intent(in)    :: cube

integer :: nx

self%gtype = gridtype

if (self%gtype=='cs') then

  ! Set fempsgrid for a cubed sphere geometry
  ! -----------------------------------------
  self%ngrids = 6

  if (.not.present(cube)) then
    self%n0 = 3 !Default value (c96)
  else
    self%n0 = cube/2**5
  endif
  nx = self%n0*(2**(self%ngrids-1))

  self%nfacex = 6*nx*nx
  self%nedgex = 2*self%nfacex
  self%nvertx = self%nfacex + 2

  self%dimfnxtf = 4
  self%dimeoff  = 4
  self%dimvoff  = 4
  self%dimfnxte = 2
  self%dimvofe  = 2
  self%dimfofv  = 4
  self%dimeofv  = 4

elseif (self%gtype=='ih') then

  ! Set fempsgrid for a icosahedral hexagons geometry
  ! -------------------------------------------------
  self%ngrids = 7

  self%nfacex=5*2**(2*self%ngrids-1)+2
  self%nedgex=15*2**(2*self%ngrids-1)
  self%nvertx=5*2**(2*self%ngrids)

  self%dimfnxtf = 6
  self%dimeoff  = 6
  self%dimvoff  = 6
  self%dimfnxte = 2
  self%dimvofe  = 2
  self%dimfofv  = 3
  self%dimeofv  = 3

elseif (self%gtype=='di') then

  print*, "Diamong grid cell not implemented"
  stop

else

  print*, "Grid type should be cs (cubed-sphere), ih (icosahedral hexagons) or di (diamonds)"
  stop

endif

! Allocate grid variables
allocate(self%nface   (self%ngrids))
allocate(self%nedge   (self%ngrids))
allocate(self%nvert   (self%ngrids))
allocate(self%neoff   (self%nfacex,self%ngrids))
allocate(self%neofv   (self%nvertx,self%ngrids))
allocate(self%fnxtf   (self%nfacex,self%dimfnxtf,self%ngrids))
allocate(self%eoff    (self%nfacex,self%dimeoff ,self%ngrids))
allocate(self%voff    (self%nfacex,self%dimvoff ,self%ngrids))
allocate(self%fnxte   (self%nedgex,self%dimfnxte,self%ngrids))
allocate(self%vofe    (self%nedgex,self%dimvofe ,self%ngrids))
allocate(self%fofv    (self%nvertx,self%dimfofv ,self%ngrids))
allocate(self%eofv    (self%nvertx,self%dimeofv ,self%ngrids))
allocate(self%flong   (self%nfacex,self%ngrids))
allocate(self%flat    (self%nfacex,self%ngrids))
allocate(self%vlong   (self%nvertx,self%ngrids))
allocate(self%vlat    (self%nvertx,self%ngrids))
allocate(self%farea   (self%nfacex,self%ngrids))
allocate(self%ldist   (self%nedgex,self%ngrids))
allocate(self%ddist   (self%nedgex,self%ngrids))

end subroutine setup

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fempsgrid), intent(inout) :: self

! Deallocate grid variables
if(allocated(self%nface)) deallocate(self%nface)
if(allocated(self%nedge)) deallocate(self%nedge)
if(allocated(self%nvert)) deallocate(self%nvert)
if(allocated(self%neoff)) deallocate(self%neoff)
if(allocated(self%neofv)) deallocate(self%neofv)
if(allocated(self%fnxtf)) deallocate(self%fnxtf)
if(allocated(self%eoff )) deallocate(self%eoff )
if(allocated(self%voff )) deallocate(self%voff )
if(allocated(self%fnxte)) deallocate(self%fnxte)
if(allocated(self%vofe )) deallocate(self%vofe )
if(allocated(self%fofv )) deallocate(self%fofv )
if(allocated(self%eofv )) deallocate(self%eofv )
if(allocated(self%flong)) deallocate(self%flong)
if(allocated(self%flat )) deallocate(self%flat )
if(allocated(self%vlong)) deallocate(self%vlong)
if(allocated(self%vlat )) deallocate(self%vlat )
if(allocated(self%farea)) deallocate(self%farea)
if(allocated(self%ldist)) deallocate(self%ldist)
if(allocated(self%ddist)) deallocate(self%ddist)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine writegrid(self,filename)

implicit none
class(fempsgrid), intent(in) :: self
character(len=*),  intent(in) :: filename

integer :: ncid, vc, varid(1000)
integer :: ngrids_dimid, nfacex_dimid, nedgex_dimid, nvertx_dimid
integer :: dimfnxtf_dimid, dimeoff_dimid, dimvoff_dimid, dimfnxte_dimid, &
           dimvofe_dimid, dimfofv_dimid, dimeofv_dimid


! Create file
! -----------
call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_CLOBBER), ncid), "nf90_create" )


! Define dimensions
! -----------------
call nccheck( nf90_def_dim(ncid, "ngrids",   self%ngrids,     ngrids_dimid), "nf90_def_dim ngrids"  )
call nccheck( nf90_def_dim(ncid, "nfacex",   self%nfacex,     nfacex_dimid), "nf90_def_dim nfacex"  )
call nccheck( nf90_def_dim(ncid, "nvertx",   self%nvertx,     nvertx_dimid), "nf90_def_dim nvertx"  )
call nccheck( nf90_def_dim(ncid, "nedgex",   self%nedgex,     nedgex_dimid), "nf90_def_dim nedgex"  )
call nccheck( nf90_def_dim(ncid, "dimfnxtf", self%dimfnxtf, dimfnxtf_dimid), "nf90_def_dim dimfnxtf")
call nccheck( nf90_def_dim(ncid, "dimeoff",  self%dimeoff,   dimeoff_dimid), "nf90_def_dim dimeoff" )
call nccheck( nf90_def_dim(ncid, "dimvoff",  self%dimvoff,   dimvoff_dimid), "nf90_def_dim dimvoff" )
call nccheck( nf90_def_dim(ncid, "dimfnxte", self%dimfnxte, dimfnxte_dimid), "nf90_def_dim dimfnxte")
call nccheck( nf90_def_dim(ncid, "dimvofe",  self%dimvofe,   dimvofe_dimid), "nf90_def_dim dimvofe" )
call nccheck( nf90_def_dim(ncid, "dimfofv",  self%dimfofv,   dimfofv_dimid), "nf90_def_dim dimfofv" )
call nccheck( nf90_def_dim(ncid, "dimeofv",  self%dimeofv,   dimeofv_dimid), "nf90_def_dim dimeofv" )


! Define variables
! ----------------
vc = 1

call nccheck( nf90_def_var(ncid, "nface", NF90_FLOAT, (/ ngrids_dimid /), varid(vc)), "nf90_def_var nface" ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "nedge", NF90_FLOAT, (/ ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "nvert", NF90_FLOAT, (/ ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "neoff", NF90_FLOAT, (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "neofv", NF90_FLOAT, (/ nvertx_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "fnxtf", NF90_FLOAT, (/ nfacex_dimid, dimfnxtf_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "eoff" , NF90_FLOAT, (/ nfacex_dimid, dimeoff_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "voff" , NF90_FLOAT, (/ nfacex_dimid, dimvoff_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "fnxte", NF90_FLOAT, (/ nedgex_dimid, dimfnxte_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "vofe" , NF90_FLOAT, (/ nedgex_dimid, dimvofe_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "fofv" , NF90_FLOAT, (/ nvertx_dimid, dimfofv_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "eofv" , NF90_FLOAT, (/ nvertx_dimid, dimeofv_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "flong", NF90_FLOAT, (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "flat" , NF90_FLOAT, (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "vlong", NF90_FLOAT, (/ nvertx_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "vlat" , NF90_FLOAT, (/ nvertx_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "farea", NF90_FLOAT, (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "ldist", NF90_FLOAT, (/ nedgex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "ddist", NF90_FLOAT, (/ nedgex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1


! End define
! ----------
call nccheck( nf90_enddef(ncid), "nf90_enddef" )


! Write variables
! ---------------
vc = 1
call nccheck( nf90_put_var( ncid, varid(vc), self%nface ), "nf90_put_var nface" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%nedge ), "nf90_put_var nedge" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%nvert ), "nf90_put_var nvert" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%neoff ), "nf90_put_var neoff" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%neofv ), "nf90_put_var neofv" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%fnxtf ), "nf90_put_var fnxtf" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%eoff  ), "nf90_put_var eoff " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%voff  ), "nf90_put_var voff " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%fnxte ), "nf90_put_var fnxte" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%vofe  ), "nf90_put_var vofe " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%fofv  ), "nf90_put_var fofv " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%eofv  ), "nf90_put_var eofv " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%flong ), "nf90_put_var flong" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%flat  ), "nf90_put_var flat " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%vlong ), "nf90_put_var vlong" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%vlat  ), "nf90_put_var vlat " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%farea ), "nf90_put_var farea" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%ldist ), "nf90_put_var ldist" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%ddist ), "nf90_put_var ddist" ); vc = vc + 1

call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine writegrid

! --------------------------------------------------------------------------------------------------

subroutine nccheck(status,iam)

use netcdf

implicit none
integer,                    intent (in) :: status
character(len=*), optional, intent (in) :: iam

character(len=1024) :: error_descr

if(status /= nf90_noerr) then

  error_descr = "NetCDF error, aborting ... "

  if (present(iam)) then
    error_descr = trim(error_descr)//", "//trim(iam)
  endif

  error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))

  print*, "Aborting: ", trim(error_descr)

  call abort()

end if

end subroutine nccheck

! --------------------------------------------------------------------------------------------------

end module femps_grid_mod

! Conventions
!
! Latitude and longitude in radians.
! Area and lengths are for the unit sphere.
!
! 1. fnxtf and eoff are ordered anticlockwise and such that
! the i'th edge lies between the face in question and its i'th
! neighbour.
!
! 2. voff are ordered anticlockwise such that the k'th vertex
! is the common vertex of the k'th and (k+1)'th edge.
!
! 3. eofv are ordered anticlockwise.
!
! 4. fofv are ordered anticlockwise such that the k'th face lies
! between the k'th and (k+1)'th edge.
!
! 5. The positive normal direction n points from
! fnxte(:,1,:) -> fnxte(:,2,:)
! and the positive tangential direction t points from
! vofe(:,1,:) -> vofe(:,2,:)
! such that t = k x n (where k is vertical).
