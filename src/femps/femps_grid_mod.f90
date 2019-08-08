module femps_grid_mod

use femps_kinds_mod

implicit none
private
public fempsgrid

! Type to hold all the grid information
! -------------------------------------
type fempsgrid

  character(len=2) :: gtype  ! cs (cubed-sphere), ih (icosahedral hexagons), di (diamonds)

  integer :: ngrids      !Number of grids in multigrid hierarchy

  integer :: nfacex      !Number of faces on each grid
  integer :: nedgex      !Number of edges on each grid
  integer :: nvertx      !Number of vertices on each edge

  integer, allocatable, dimension(:) :: nface, nedge, nvert

  integer :: n0 = 3

  integer, allocatable, dimension(:,:) :: neoff !number of edges and vertices of each face on each grid
  integer, allocatable, dimension(:,:) :: neofv !number of edges of each vertex on each grid

  ! Dimensions for the following arrays
  ! -----------------------------------
  integer :: dimfnxtf, dimeoff, dimvoff, dimfnxte, dimvofe, dimfofv, dimeofv  !*

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

  real(kind=kind_real), allocatable, dimension(:) :: fareamin !*

  integer :: nefmx, nevmx  !*

  contains

   procedure :: setup
   procedure :: delete

end type fempsgrid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setup(self,gridtype)

implicit none
class(fempsgrid), intent(inout) :: self
character(len=2), intent(in)    :: gridtype

integer :: nx

self%gtype = gridtype

if (self%gtype=='cs') then

  ! Set fempsgrid for a cubed sphere geometry
  ! -----------------------------------------
  self%ngrids = 6

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
allocate(self%neoff   (self%nfacex,self%ngrids))
allocate(self%neofv   (self%nvertx,self%ngrids))
allocate(self%nface   (self%ngrids))
allocate(self%nedge   (self%ngrids))
allocate(self%nvert   (self%ngrids))
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
allocate(self%fareamin(self%ngrids))

end subroutine setup

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fempsgrid), intent(inout) :: self

! Deallocate grid variables
deallocate(self%neoff)
deallocate(self%neofv)
deallocate(self%nface)
deallocate(self%nedge)
deallocate(self%nvert)
deallocate(self%fnxtf)
deallocate(self%eoff)
deallocate(self%voff)
deallocate(self%fnxte)
deallocate(self%vofe)
deallocate(self%fofv)
deallocate(self%eofv)
deallocate(self%flong)
deallocate(self%flat)
deallocate(self%vlong)
deallocate(self%vlat)
deallocate(self%farea)
deallocate(self%ldist)
deallocate(self%ddist)
deallocate(self%fareamin)

end subroutine delete

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
