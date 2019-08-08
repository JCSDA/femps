module femps_inout_mod

use femps_grid_mod
use femps_operators_mod
use femps_const_mod

implicit none

private
public readgridoprs

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine readgridoprs(grid,oprs)

! To allocate array space for the grid information in module grid
! and to read the information from file

implicit none
type(fempsgrid), intent(inout) :: grid
type(fempsoprs), intent(inout) :: oprs

integer :: if0, ie0, iv0, igrid, ix, ixx, if1, if2, iv1, iv2, ie1

character(len=31) :: ygridfile
integer, parameter :: changrid = 25
integer, parameter :: channml = 20

! Read namelists
! --------------
namelist /rundata/ ygridfile

open(channml,file='poissonnml.in',delim='apostrophe')
read(channml,rundata)
close(channml)

! Open file for reading
open(changrid,file=ygridfile,FORM='UNFORMATTED')

! First read grid%ngrids
read(changrid) grid%ngrids

! Allocate nface, nedge, nvert
allocate(grid%nface(grid%ngrids), grid%nedge(grid%ngrids), grid%nvert(grid%ngrids))

! Read numbers of faces, edges and vertices on each grid
READ(changrid) grid%nface
READ(changrid) grid%nedge
READ(changrid) grid%nvert

! Find maximum values in order to allocated subsequent arrays
grid%nfacex = MAXVAL(grid%nface)
grid%nedgex = MAXVAL(grid%nedge)
grid%nvertx = MAXVAL(grid%nvert)

! Allocate neoff, neofv
allocate(grid%neoff(grid%nfacex,grid%ngrids), grid%neofv(grid%nvertx,grid%ngrids))

! Read the numbers of edges of each face and vertex on each grid
grid%neoff = 0
grid%neofv = 0
READ(changrid) ((grid%neoff(if0,igrid),          &
                    if0 = 1, grid%nface(igrid)), &
                    igrid = 1, grid%ngrids)
READ(changrid) ((grid%neofv(iv0,igrid),          &
                    iv0 = 1, grid%nvert(igrid)), &
                    igrid = 1, grid%ngrids)

! Find maximum values in order to allocate subsequent arrays
grid%nefmx = MAXVAL(grid%neoff)
grid%nevmx = MAXVAL(grid%neofv)


! Allocate connectivity arrays arrays
allocate(grid%fnxtf(grid%nfacex,grid%nefmx,grid%ngrids), grid%eoff(grid%nfacex,grid%nefmx,grid%ngrids), &
         grid%voff(grid%nfacex,grid%nefmx,grid%ngrids),  grid%fnxte(grid%nedgex,2,grid%ngrids),    &
         grid%vofe(grid%nedgex,2,grid%ngrids),      grid%fofv(grid%nvertx,grid%nevmx,grid%ngrids), &
         grid%eofv(grid%nvertx,grid%nevmx,grid%ngrids))

! Read the connectivity arrays
READ(changrid) (((grid%fnxtf(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, grid%nefmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%eoff(if0,ix,igrid),           &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, grid%nefmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%voff(if0,ix,igrid),           &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, grid%nefmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%fnxte(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%vofe(ie0,ix,igrid),           &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%fofv(iv0,ix,igrid),           &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, grid%nevmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%eofv(iv0,ix,igrid),           &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, grid%nevmx),            &
                     igrid = 1, grid%ngrids)


! Allocate the geometrical information arrays
allocate(grid%flong(grid%nfacex,grid%ngrids), grid%flat(grid%nfacex,grid%ngrids),  &
         grid%vlong(grid%nvertx,grid%ngrids), grid%vlat(grid%nvertx,grid%ngrids),  &
         grid%farea(grid%nfacex,grid%ngrids), oprs%varea(grid%nvertx,grid%ngrids), &
         grid%ldist(grid%nedgex,grid%ngrids), grid%ddist(grid%nedgex,grid%ngrids), &
         grid%fareamin(grid%ngrids))

! Read the geometrical information arrays
READ(changrid) ((grid%flong(if0,igrid),              &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%flat(if0,igrid),               &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%vlong(iv0,igrid),              &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%vlat(iv0,igrid),               &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%farea(if0,igrid),              &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%varea(iv0,igrid),              &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%ldist(ie0,igrid),              &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%ddist(ie0,igrid),              &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)

! Dimensionalize
grid%farea = grid%farea*rearth*rearth
oprs%varea = oprs%varea*rearth*rearth
grid%ldist = grid%ldist*rearth
grid%ddist = grid%ddist*rearth

! Determine smallest face area on each grid
do igrid = 1, grid%ngrids
  grid%fareamin(igrid) = MINVAL(grid%farea(1:grid%nface(igrid),igrid))
enddo


! Allocate arrays for size of operator stencils
allocate(oprs%nlsten(grid%nfacex,grid%ngrids), oprs%nmsten(grid%nedgex,grid%ngrids), &
         oprs%njsten(grid%nvertx,grid%ngrids), oprs%nhsten(grid%nedgex,grid%ngrids), &
         oprs%nrsten(grid%nfacex,grid%ngrids), oprs%nrxsten(grid%nvertx,grid%ngrids), &
         oprs%nwsten(grid%nedgex,grid%ngrids), oprs%ntsten(grid%nfacex,grid%ngrids))

! Read the sizes of the operator stencils on each grid
oprs%nlsten = 0
oprs%nmsten = 0
oprs%njsten = 0
oprs%nhsten = 0
oprs%nrsten = 0
oprs%nrxsten = 0
oprs%nwsten = 0
oprs%ntsten = 0
READ(changrid) ((oprs%nlsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nmsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%njsten(iv0,igrid),             &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nhsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nrsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nrxsten(iv0,igrid),            &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nwsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%ntsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)

! Find maximum values in order to allocate subsequent arrays
oprs%nlsmx = MAXVAL(oprs%nlsten)
oprs%nmsmx = MAXVAL(oprs%nmsten)
oprs%njsmx = MAXVAL(oprs%njsten)
oprs%nhsmx = MAXVAL(oprs%nhsten)
oprs%nrsmx = MAXVAL(oprs%nrsten)
oprs%nrxsmx = MAXVAL(oprs%nrxsten)
oprs%nwsmx = MAXVAL(oprs%nwsten)
oprs%ntsmx = MAXVAL(oprs%ntsten)

PRINT *,'Maximum stencil sizes:'
PRINT *,'massL ...  ',oprs%nlsmx
PRINT *,'massM ...  ',oprs%nmsmx
PRINT *,'HodgeJ ... ',oprs%njsmx
PRINT *,'HodgeH ... ',oprs%nhsmx
PRINT *,'operR ...  ',oprs%nrxsmx
PRINT *,'operW ...  ',oprs%nwsmx
PRINT *,'operT ...  ',oprs%ntsmx
PRINT *,' '


! Allocate arrays for operator stencils and coefficients
allocate(oprs%lsten(grid%nfacex,oprs%nlsmx,grid%ngrids), oprs%msten(grid%nedgex,oprs%nmsmx,grid%ngrids), &
         oprs%jsten(grid%nvertx,oprs%njsmx,grid%ngrids), oprs%hsten(grid%nedgex,oprs%nhsmx,grid%ngrids), &
         oprs%rsten(grid%nfacex,oprs%nrsmx,grid%ngrids), oprs%rxsten(grid%nvertx,oprs%nrxsmx,grid%ngrids), &
         oprs%wsten(grid%nedgex,oprs%nwsmx,grid%ngrids), oprs%tsten(grid%nfacex,oprs%ntsmx,grid%ngrids))
allocate(oprs%lmass(grid%nfacex,oprs%nlsmx,grid%ngrids), oprs%mmass(grid%nedgex,oprs%nmsmx,grid%ngrids), &
         oprs%jstar(grid%nvertx,oprs%njsmx,grid%ngrids), oprs%hstar(grid%nedgex,oprs%nhsmx,grid%ngrids), &
         oprs%rcoeff(grid%nfacex,oprs%nrsmx,grid%ngrids), oprs%rxcoeff(grid%nvertx,oprs%nrxsmx,grid%ngrids), &
         oprs%wcoeff(grid%nedgex,oprs%nwsmx,grid%ngrids), oprs%tcoeff(grid%nfacex,oprs%ntsmx,oprs%ntsmx,grid%ngrids), &
         oprs%jlump(grid%nvertx,grid%ngrids), oprs%mlump(grid%nedgex,grid%ngrids), oprs%hlump(grid%nedgex,grid%ngrids))

! Read the operator stencils and coefficients
READ(changrid) (((oprs%lsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nlsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%msten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nmsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%jsten(iv0,ix,igrid),          &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%njsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%hsten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nhsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nrsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rxsten(iv0,ix,igrid),         &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%nrxsmx),           &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%wsten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nwsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%tsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%ntsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%lmass(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nlsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%mmass(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nmsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%jstar(iv0,ix,igrid),          &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%njsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%hstar(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nhsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rcoeff(if0,ix,igrid),         &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nrsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rxcoeff(iv0,ix,igrid),        &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%nrxsmx),           &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%wcoeff(ie0,ix,igrid),         &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nwsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) ((((oprs%tcoeff(if0,ix,ixx,igrid),    &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%ntsmx),            &
                     ixx = 1, oprs%ntsmx),           &
                     igrid = 1, grid%ngrids)

! Dimensionalize
oprs%lmass = oprs%lmass/(rearth*rearth)

! Construct the tables eoffin and eofvin
allocate(oprs%eoffin(grid%nfacex,grid%nefmx,grid%ngrids), oprs%eofvin(grid%nvertx,grid%nevmx,grid%ngrids))
do igrid = 1, grid%ngrids

  do if1 = 1, grid%nface(igrid)
    do ix = 1, grid%neoff(if1,igrid)
      ie1 = grid%eoff(if1,ix,igrid)
      if2 = grid%fnxte(ie1,1,igrid)
      if (if1 == if2) then
        ! This edge points out of face if1
        oprs%eoffin(if1,ix,igrid) = -1.0d0
      else
        ! This edge points into face if1
        oprs%eoffin(if1,ix,igrid) = 1.0d0
      endif
    enddo
  enddo

  do iv1 = 1, grid%nvert(igrid)
    do ix = 1, grid%neofv(iv1,igrid)
      ie1 = grid%eofv(iv1,ix,igrid)
      iv2 = grid%vofe(ie1,1,igrid)
      if (iv1 == iv2) then
        ! This edge points away from vertex iv1
        oprs%eofvin(iv1,ix,igrid) = -1.0d0
      else
        ! This edge points towards vertex iv1
        oprs%eofvin(iv1,ix,igrid) = 1.0d0
      endif
    enddo
  enddo

enddo


! Allocate array for size of restriction stencil
allocate(oprs%ninj(grid%nfacex,grid%ngrids-1))

! Read the size of the restriction stencil on each grid
oprs%ninj = 0
READ(changrid) ((oprs%ninj(if0,igrid),              &
                    if0 = 1, grid%nface(igrid)),    &
                    igrid = 1, grid%ngrids-1)

! Find maximum value in order to allocate subsequent arrays
oprs%ninjmx = MAXVAL(oprs%ninj)

! Allocate arrays for restriction stencils and weights
allocate(oprs%injsten(grid%nfacex,oprs%ninjmx,grid%ngrids-1))
allocate(oprs%injwgt(grid%nfacex,oprs%ninjmx,grid%ngrids-1))

! Read the restriction stencil and weights
READ(changrid) (((oprs%injsten(if0,ix,igrid),       &
                    if0 = 1, grid%nface(igrid)),    &
                    ix = 1, oprs%ninjmx),           &
                    igrid = 1, grid%ngrids-1)
READ(changrid) (((oprs%injwgt(if0,ix,igrid),        &
                    if0 = 1, grid%nface(igrid)),    &
                    ix = 1, oprs%ninjmx),           &
                    igrid = 1, grid%ngrids-1)

allocate(oprs%lapdiag(grid%nfacex,grid%ngrids))
allocate(oprs%underrel(grid%ngrids))

end subroutine readgridoprs

! --------------------------------------------------------------------------------------------------

subroutine readoprs(grid,oprs)

! To allocate array space for the grid information in module grid
! and to read the information from file

implicit none
type(fempsgrid), intent(inout) :: grid
type(fempsoprs), intent(inout) :: oprs

integer :: if0, ie0, iv0, igrid, ix, ixx, if1, if2, iv1, iv2, ie1

character(len=31) :: ygridfile
integer, parameter :: changrid = 25
integer, parameter :: channml = 20

! Read namelists
! --------------
namelist /rundata/ ygridfile

open(channml,file='poissonnml.in',delim='apostrophe')
read(channml,rundata)
close(channml)

! Open file for reading
open(changrid,file=ygridfile,FORM='UNFORMATTED')

! First read grid%ngrids
read(changrid) grid%ngrids

! Allocate nface, nedge, nvert
allocate(grid%nface(grid%ngrids), grid%nedge(grid%ngrids), grid%nvert(grid%ngrids))

! Read numbers of faces, edges and vertices on each grid
READ(changrid) grid%nface
READ(changrid) grid%nedge
READ(changrid) grid%nvert

! Find maximum values in order to allocated subsequent arrays
grid%nfacex = MAXVAL(grid%nface)
grid%nedgex = MAXVAL(grid%nedge)
grid%nvertx = MAXVAL(grid%nvert)

! Allocate neoff, neofv
allocate(grid%neoff(grid%nfacex,grid%ngrids), grid%neofv(grid%nvertx,grid%ngrids))

! Read the numbers of edges of each face and vertex on each grid
grid%neoff = 0
grid%neofv = 0
READ(changrid) ((grid%neoff(if0,igrid),          &
                    if0 = 1, grid%nface(igrid)), &
                    igrid = 1, grid%ngrids)
READ(changrid) ((grid%neofv(iv0,igrid),          &
                    iv0 = 1, grid%nvert(igrid)), &
                    igrid = 1, grid%ngrids)

! Find maximum values in order to allocate subsequent arrays
grid%nefmx = MAXVAL(grid%neoff)
grid%nevmx = MAXVAL(grid%neofv)


! Allocate connectivity arrays arrays
allocate(grid%fnxtf(grid%nfacex,grid%nefmx,grid%ngrids))
allocate(grid%eoff (grid%nfacex,grid%nefmx,grid%ngrids))
allocate(grid%voff (grid%nfacex,grid%nefmx,grid%ngrids))
allocate(grid%fnxte(grid%nedgex,2,grid%ngrids))
allocate(grid%vofe (grid%nedgex,2,grid%ngrids))
allocate(grid%fofv (grid%nvertx,grid%nevmx,grid%ngrids))
allocate(grid%eofv (grid%nvertx,grid%nevmx,grid%ngrids))

! Read the connectivity arrays
READ(changrid) (((grid%fnxtf(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, grid%nefmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%eoff(if0,ix,igrid),           &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, grid%nefmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%voff(if0,ix,igrid),           &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, grid%nefmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%fnxte(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%vofe(ie0,ix,igrid),           &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%fofv(iv0,ix,igrid),           &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, grid%nevmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((grid%eofv(iv0,ix,igrid),           &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, grid%nevmx),            &
                     igrid = 1, grid%ngrids)


! Allocate the geometrical information arrays
allocate(grid%flong(grid%nfacex,grid%ngrids), grid%flat(grid%nfacex,grid%ngrids),  &
         grid%vlong(grid%nvertx,grid%ngrids), grid%vlat(grid%nvertx,grid%ngrids),  &
         grid%farea(grid%nfacex,grid%ngrids), oprs%varea(grid%nvertx,grid%ngrids), &
         grid%ldist(grid%nedgex,grid%ngrids), grid%ddist(grid%nedgex,grid%ngrids), &
         grid%fareamin(grid%ngrids))

! Read the geometrical information arrays
READ(changrid) ((grid%flong(if0,igrid),              &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%flat(if0,igrid),               &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%vlong(iv0,igrid),              &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%vlat(iv0,igrid),               &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%farea(if0,igrid),              &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%varea(iv0,igrid),              &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%ldist(ie0,igrid),              &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((grid%ddist(ie0,igrid),              &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)

! Dimensionalize
grid%farea = grid%farea*rearth*rearth
oprs%varea = oprs%varea*rearth*rearth
grid%ldist = grid%ldist*rearth
grid%ddist = grid%ddist*rearth

! Determine smallest face area on each grid
do igrid = 1, grid%ngrids
  grid%fareamin(igrid) = MINVAL(grid%farea(1:grid%nface(igrid),igrid))
enddo


! Allocate arrays for size of operator stencils
allocate(oprs%nlsten(grid%nfacex,grid%ngrids), oprs%nmsten(grid%nedgex,grid%ngrids), &
         oprs%njsten(grid%nvertx,grid%ngrids), oprs%nhsten(grid%nedgex,grid%ngrids), &
         oprs%nrsten(grid%nfacex,grid%ngrids), oprs%nrxsten(grid%nvertx,grid%ngrids), &
         oprs%nwsten(grid%nedgex,grid%ngrids), oprs%ntsten(grid%nfacex,grid%ngrids))

! Read the sizes of the operator stencils on each grid
oprs%nlsten = 0
oprs%nmsten = 0
oprs%njsten = 0
oprs%nhsten = 0
oprs%nrsten = 0
oprs%nrxsten = 0
oprs%nwsten = 0
oprs%ntsten = 0
READ(changrid) ((oprs%nlsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nmsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%njsten(iv0,igrid),             &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nhsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nrsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nrxsten(iv0,igrid),            &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%nwsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((oprs%ntsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)

! Find maximum values in order to allocate subsequent arrays
oprs%nlsmx = MAXVAL(oprs%nlsten)
oprs%nmsmx = MAXVAL(oprs%nmsten)
oprs%njsmx = MAXVAL(oprs%njsten)
oprs%nhsmx = MAXVAL(oprs%nhsten)
oprs%nrsmx = MAXVAL(oprs%nrsten)
oprs%nrxsmx = MAXVAL(oprs%nrxsten)
oprs%nwsmx = MAXVAL(oprs%nwsten)
oprs%ntsmx = MAXVAL(oprs%ntsten)

PRINT *,'Maximum stencil sizes:'
PRINT *,'massL ...  ',oprs%nlsmx
PRINT *,'massM ...  ',oprs%nmsmx
PRINT *,'HodgeJ ... ',oprs%njsmx
PRINT *,'HodgeH ... ',oprs%nhsmx
PRINT *,'operR ...  ',oprs%nrxsmx
PRINT *,'operW ...  ',oprs%nwsmx
PRINT *,'operT ...  ',oprs%ntsmx
PRINT *,' '


! Allocate arrays for operator stencils and coefficients
allocate(oprs%lsten(grid%nfacex,oprs%nlsmx,grid%ngrids), oprs%msten(grid%nedgex,oprs%nmsmx,grid%ngrids), &
         oprs%jsten(grid%nvertx,oprs%njsmx,grid%ngrids), oprs%hsten(grid%nedgex,oprs%nhsmx,grid%ngrids), &
         oprs%rsten(grid%nfacex,oprs%nrsmx,grid%ngrids), oprs%rxsten(grid%nvertx,oprs%nrxsmx,grid%ngrids), &
         oprs%wsten(grid%nedgex,oprs%nwsmx,grid%ngrids), oprs%tsten(grid%nfacex,oprs%ntsmx,grid%ngrids))
allocate(oprs%lmass(grid%nfacex,oprs%nlsmx,grid%ngrids), oprs%mmass(grid%nedgex,oprs%nmsmx,grid%ngrids), &
         oprs%jstar(grid%nvertx,oprs%njsmx,grid%ngrids), oprs%hstar(grid%nedgex,oprs%nhsmx,grid%ngrids), &
         oprs%rcoeff(grid%nfacex,oprs%nrsmx,grid%ngrids), oprs%rxcoeff(grid%nvertx,oprs%nrxsmx,grid%ngrids), &
         oprs%wcoeff(grid%nedgex,oprs%nwsmx,grid%ngrids), oprs%tcoeff(grid%nfacex,oprs%ntsmx,oprs%ntsmx,grid%ngrids), &
         oprs%jlump(grid%nvertx,grid%ngrids), oprs%mlump(grid%nedgex,grid%ngrids), oprs%hlump(grid%nedgex,grid%ngrids))

! Read the operator stencils and coefficients
READ(changrid) (((oprs%lsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nlsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%msten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nmsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%jsten(iv0,ix,igrid),          &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%njsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%hsten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nhsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nrsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rxsten(iv0,ix,igrid),         &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%nrxsmx),           &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%wsten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nwsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%tsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%ntsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%lmass(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nlsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%mmass(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nmsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%jstar(iv0,ix,igrid),          &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%njsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%hstar(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nhsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rcoeff(if0,ix,igrid),         &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%nrsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%rxcoeff(iv0,ix,igrid),        &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, oprs%nrxsmx),           &
                     igrid = 1, grid%ngrids)
READ(changrid) (((oprs%wcoeff(ie0,ix,igrid),         &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, oprs%nwsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) ((((oprs%tcoeff(if0,ix,ixx,igrid),    &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, oprs%ntsmx),            &
                     ixx = 1, oprs%ntsmx),           &
                     igrid = 1, grid%ngrids)

! Dimensionalize
oprs%lmass = oprs%lmass/(rearth*rearth)

! Construct the tables eoffin and eofvin
allocate(oprs%eoffin(grid%nfacex,grid%nefmx,grid%ngrids), oprs%eofvin(grid%nvertx,grid%nevmx,grid%ngrids))
do igrid = 1, grid%ngrids

  do if1 = 1, grid%nface(igrid)
    do ix = 1, grid%neoff(if1,igrid)
      ie1 = grid%eoff(if1,ix,igrid)
      if2 = grid%fnxte(ie1,1,igrid)
      if (if1 == if2) then
        ! This edge points out of face if1
        oprs%eoffin(if1,ix,igrid) = -1.0d0
      else
        ! This edge points into face if1
        oprs%eoffin(if1,ix,igrid) = 1.0d0
      endif
    enddo
  enddo

  do iv1 = 1, grid%nvert(igrid)
    do ix = 1, grid%neofv(iv1,igrid)
      ie1 = grid%eofv(iv1,ix,igrid)
      iv2 = grid%vofe(ie1,1,igrid)
      if (iv1 == iv2) then
        ! This edge points away from vertex iv1
        oprs%eofvin(iv1,ix,igrid) = -1.0d0
      else
        ! This edge points towards vertex iv1
        oprs%eofvin(iv1,ix,igrid) = 1.0d0
      endif
    enddo
  enddo

enddo


! Allocate array for size of restriction stencil
allocate(oprs%ninj(grid%nfacex,grid%ngrids-1))

! Read the size of the restriction stencil on each grid
oprs%ninj = 0
READ(changrid) ((oprs%ninj(if0,igrid),              &
                    if0 = 1, grid%nface(igrid)),    &
                    igrid = 1, grid%ngrids-1)

! Find maximum value in order to allocate subsequent arrays
oprs%ninjmx = MAXVAL(oprs%ninj)

! Allocate arrays for restriction stencils and weights
allocate(oprs%injsten(grid%nfacex,oprs%ninjmx,grid%ngrids-1))
allocate(oprs%injwgt(grid%nfacex,oprs%ninjmx,grid%ngrids-1))

! Read the restriction stencil and weights
READ(changrid) (((oprs%injsten(if0,ix,igrid),       &
                    if0 = 1, grid%nface(igrid)),    &
                    ix = 1, oprs%ninjmx),           &
                    igrid = 1, grid%ngrids-1)
READ(changrid) (((oprs%injwgt(if0,ix,igrid),        &
                    if0 = 1, grid%nface(igrid)),    &
                    ix = 1, oprs%ninjmx),           &
                    igrid = 1, grid%ngrids-1)

allocate(oprs%lapdiag(grid%nfacex,grid%ngrids))
allocate(oprs%underrel(grid%ngrids))

end subroutine readoprs

! --------------------------------------------------------------------------------------------------

end module femps_inout_mod
