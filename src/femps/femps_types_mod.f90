module femps_types_mod

use femps_kinds_mod

implicit none
private
public fempsgrid, fempspbops

! Type to hold all the grid information
! -------------------------------------
type fempsgrid

  ! Number of grids in multigrid hierarchy
  ! --------------------------------------
  integer :: ngrids

  ! nface: number of faces on each grid
  ! nedge: number of edges on each grid
  ! nvert: number of vertices on each edge
  ! --------------------------------------
  integer :: nfacex, nedgex, nvertx
  integer, allocatable, dimension(:) :: nface, nedge, nvert

  ! neoff: number of edges and vertices of each face on each grid
  ! neofv: number of edges of each vertex on each grid
  ! --------------------------------------------------
  integer, allocatable, dimension(:,:) :: neoff, neofv

  ! demensions for the following arrays
  ! -----------------------------------
  integer :: dimfnxtf, dimeoff, dimvoff, dimfnxte, dimvofe, dimfofv, dimeofv

  ! CONNECTIVITY
  ! fnxtf: faces next to each face on each grid
  ! eoff : edges of each face on each grid
  ! voff : vertices of each face on each grid
  ! fnxte: faces either side of each edge on each grid
  ! vofe : vertices at the ends of each edge on each grid
  ! fofv : faces around each vertex on each grid
  ! eofv : edges incident on each vertex on each grid
  ! -------------------------------------------------
  integer, allocatable, dimension(:,:,:) :: fnxtf, eoff, voff, fnxte, &
                                            vofe, fofv, eofv

  ! COORDINATES AND GEOMETRICAL INFORMATION
  ! flong: longitude of faces on each grid
  ! flat : latitude of faces on each grid
  ! vlong: longitude of vertices on each grid
  ! vlat : latitude of vertices on each grid
  ! farea: area of faces on each grid
  ! ldist: primal edge length, i.e. distance between neighbouring vertices
  ! ddist: dual edge length, i.e. distance between neighbouring face centres
  ! ------------------------------------------------------------------------
  real(kind=kind_real), allocatable, dimension(:,:) :: flong, flat, vlong, vlat, &
                                                       farea, ldist, ddist

  real(kind=kind_real), allocatable, dimension(:) :: fareamin

  ! Grid file
  character*31 :: ygridfile = 'gridopermap_cube_0000013824.dat'

  integer :: nefmx, nevmx

  contains
   procedure :: set_grid_cs
   procedure :: set_grid_ih
   procedure :: allocate_grid
   procedure :: deallocate_grid

end type fempsgrid

! Type to hold the pre-built operators
! ------------------------------------
type fempspbops

  ! Dimensions
  ! ----------
  integer :: nlsmx, nmsmx, njsmx, nhsmx, nrsmx, nrxsmx, nwsmx, ntsmx, nxmisx, ninjmx

  ! varea: area of dual faces on each grid
  ! --------------------------------------
  real(kind=kind_real), allocatable :: varea(:,:)

  ! eoffin(f,j,:)   Indicates whether the normal at the j'th edge is
  !                 inward or outward relative to face f.
  ! eofvin(v,j,:)   Indicates whether the tangent at the j'th edge is
  !                 inward or outward relative to vertex v.
  ! -------------------------------------------------------
  integer, allocatable :: eoffin(:,:,:), eofvin(:,:,:)

  ! HODGE STAR, MASS MATRIX, AND RELATED OPERATORS
  ! nlsten   : number of faces in stencil for L mass matrix
  ! lsten    : stencil for L mass matrix
  ! lmass    : coefficients for L mass matrix
  ! nmsten   : number of faces in stencil for M mass matrix
  ! msten    : stencil for M mass matrix
  ! mmass    : coefficients for M mass matrix
  ! njsten   : number of vertices in stencil for J operator
  ! jsten    : stencil for J operator
  ! jstar    : coefficients for J operator
  ! nhsten   : number of edges in stencil for H operator
  ! hsten    : stencil for H operator
  ! hstar    : coefficients for H operator
  ! nrsten   : number of vertices in stencil for R operator (= self%neoff)
  ! rsten    : stencil for R operator (= voff)
  ! rcoeff   : coefficients for R operator
  ! nrxsten  : number of faces in stencil for R transpose operator (= self%neofv)
  ! rxsten   : stencil for R transpose operator (= fofv)
  ! rxcoeff  : coefficients for R transpose operator
  ! nwsten   : number of edges in stencil for W operator
  ! wsten    : stencil for W operator
  ! wcoeff   : coefficients for W operator
  ! ntsten   : number of edges in stencel for T operator
  ! tsten    : stencil for T operator
  ! tcoeff   : coefficients for T operator
  ! jlump    : coefficients of lumped J matrix
  ! mlump    : coefficients of lumped M matrix
  ! hlump    : coefficients of lumped H matrix
  ! nxminvten: number of edges in stencil for approximate inverse of M
  ! xminvsten: stencil for approximate inverse of M
  ! xminv    : coefficients for approximate inverse of M
  ! velcoeff : coefficients to reconstruct velocity vector in cells
  ! ---------------------------------------------------------------
  integer, allocatable, dimension(:,:) :: nlsten, nmsten, njsten, &
                                          nhsten, nrsten, nrxsten, &
                                          nwsten, ntsten, nxminvsten

  integer, allocatable, dimension(:,:,:) :: lsten, msten, jsten, &
                                            hsten, rsten, rxsten, &
                                            wsten, tsten, xminvsten

  real(kind=kind_real), allocatable, dimension(:,:)     :: jlump, mlump, hlump

  real(kind=kind_real), allocatable, dimension(:,:,:)   :: lmass, mmass, jstar, &
                                                           hstar, rcoeff, rxcoeff, &
                                                           wcoeff, xminv, velcoeff

  real(kind=kind_real), allocatable, dimension(:,:,:,:) :: tcoeff


  ! RESTRICTION AND PROLONGATION OPERATORS FOR MULTIGRID
  ! ninj   : number of faces in stencil for restriction operator
  ! injsten: stencil for restriction operator
  ! injwgt : weights for restriction operator
  ! -----------------------------------------
  integer, allocatable, dimension(:,:)   :: ninj
  integer, allocatable, dimension(:,:,:) :: injsten
  real(kind=kind_real), allocatable, dimension(:,:,:) :: injwgt

  real(kind=kind_real), allocatable :: lapdiag(:,:), underrel(:)

  contains
   procedure :: allocate_pbops
   procedure :: deallocate_pbops

end type fempspbops

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine set_grid_cs(self)

implicit none
class(fempsgrid), intent(inout) :: self

integer :: n0, nx

! Set fempsgrid for a cubed sphere geometry
! -----------------------------------------
self%ngrids = 6

n0 = 3
nx = n0*(2**(self%ngrids-1))

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

end subroutine set_grid_cs

! --------------------------------------------------------------------------------------------------

subroutine set_grid_ih(self)

implicit none
class(fempsgrid), intent(inout) :: self

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

end subroutine set_grid_ih

! --------------------------------------------------------------------------------------------------

subroutine allocate_grid(self)

implicit none
class(fempsgrid), intent(inout) :: self

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

end subroutine allocate_grid

! --------------------------------------------------------------------------------------------------

subroutine deallocate_grid(self)

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

end subroutine deallocate_grid

! --------------------------------------------------------------------------------------------------

subroutine allocate_pbops(self,grid)

implicit none
class(fempspbops), intent(inout) :: self
type(fempsgrid),   intent(in)    :: grid

allocate(self%lapdiag(grid%nfacex,grid%ngrids))
allocate(self%underrel(grid%ngrids))

end subroutine allocate_pbops

! --------------------------------------------------------------------------------------------------

subroutine deallocate_pbops(self)

implicit none
class(fempspbops), intent(inout) :: self

if (allocated(self%lapdiag   )) deallocate(self%lapdiag   )
if (allocated(self%underrel  )) deallocate(self%underrel  )

end subroutine deallocate_pbops

! --------------------------------------------------------------------------------------------------

end module femps_types_mod

! Conventions
!
! Latitude and longitude in radians.
! Area and lengths are for the unit sphere.


! Conventions
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
