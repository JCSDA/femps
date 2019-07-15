module femps_types_mod

use femps_kinds_mod

implicit none
public

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

end type fempsgrid

! Type to hold the pre-build operators
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

end type fempspbops

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine set_fempsgrid_cs(grid)

implicit none
type(fempsgrid), intent(inout) :: grid

! Set fempsgrid for a cubed sphere geometry
! -----------------------------------------
grid%ngrids = 6

grid%nfacex = 6*nx*nx
grid%nedgex = 2*nfacex
grid%nvertx = nfacex + 2

grid%dimfnxtf = 4
grid%dimeoff  = 4
grid%dimvoff  = 4
grid%dimfnxte = 2
grid%dimvofe  = 2
grid%dimfofv  = 4
grid%dimeofv  = 4

end subroutine set_fempsgrid_cs

! --------------------------------------------------------------------------------------------------

subroutine set_fempsgrid_ih(grid)

implicit none
type(fempsgrid), intent(inout) :: grid

! Set fempsgrid for a icosahedral hexagons geometry
! -------------------------------------------------
grid%ngrids = 7

grid%nfacex=5*2**(2*ngrids-1)+2
grid%nedgex=15*2**(2*ngrids-1)
grid%nvertx=5*2**(2*ngrids)

grid%dimfnxtf = 6
grid%dimeoff  = 6
grid%dimvoff  = 6
grid%dimfnxte = 2
grid%dimvofe  = 2
grid%dimfofv  = 3
grid%dimeofv  = 3

end subroutine set_fempsgrid_ih

! --------------------------------------------------------------------------------------------------

subroutine allocate_fempsgrid(grid)

implicit none
type(fempsgrid), intent(inout) :: grid

! Allocate grid variables
allocate(grid%neoff(grid%nfacex,grid%ngrids))
allocate(grid%neofv(grid%nvertx,grid%ngrids))
allocate(grid%nface(grid%ngrids))
allocate(grid%nedge(grid%ngrids))
allocate(grid%nvert(grid%ngrids))
allocate(grid%fnxtf(nfacex,grid%dimfnxtf,grid%ngrids))
allocate(grid%eoff (nfacex,grid%dimeoff ,grid%ngrids))
allocate(grid%voff (nfacex,grid%dimvoff ,grid%ngrids))
allocate(grid%fnxte(nedgex,grid%dimfnxte,grid%ngrids))
allocate(grid%vofe (nedgex,grid%dimvofe ,grid%ngrids))
allocate(grid%fofv (nvertx,grid%dimfofv ,grid%ngrids))
allocate(grid%eofv (nvertx,grid%dimeofv ,grid%ngrids))
allocate(grid%flong(nfacex,ngrids))
allocate(grid%flat (nfacex,ngrids))
allocate(grid%vlong(nvertx,ngrids))
allocate(grid%vlat (nvertx,ngrids))
allocate(grid%farea(nfacex,ngrids))
allocate(grid%ldist(nedgex,ngrids))
allocate(grid%ddist(nedgex,ngrids))
allocate(grid%fareamin(grid%ngrids)

end subroutine allocate_fempsgrid

! --------------------------------------------------------------------------------------------------

subroutine deallocate_fempsgrid(grid)

implicit none
type(fempsgrid), intent(inout) :: grid

! Deallocate grid variables
deallocate(grid%neoff)
deallocate(grid%neofv)
deallocate(grid%nface)
deallocate(grid%nedge)
deallocate(grid%nvert)
deallocate(grid%fnxtf)
deallocate(grid%eoff)
deallocate(grid%voff)
deallocate(grid%fnxte)
deallocate(grid%vofe)
deallocate(grid%fofv)
deallocate(grid%eofv)
deallocate(grid%flong)
deallocate(grid%flat)
deallocate(grid%vlong)
deallocate(grid%vlat)
deallocate(grid%farea)
deallocate(grid%ldist)
deallocate(grid%ddist)
deallocate(grid%fareamin)

end subroutine deallocate_fempsgrid

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
