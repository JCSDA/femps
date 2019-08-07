module femps_inout_mod

use femps_kinds_mod
use femps_types_mod

implicit none
private
public

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine write_grid(grid)

type(fempsgrid), intent(in) :: grid

! Write out coordinates of edges for plotting
OPEN(44,FILE='primalgrid.dat',FORM='FORMATTED')
DO j = 1, nedgex
  iv = vofe(j,1,ngrids)
  long = vlong(iv,ngrids)
  lat = vlat(iv,ngrids)
  CALL ll2xyz(long,lat,x1,y1,z1)
  iv = vofe(j,2,ngrids)
  long = vlong(iv,ngrids)
  lat = vlat(iv,ngrids)
  CALL ll2xyz(long,lat,x2,y2,z2)
  WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
ENDDO
CLOSE(44)
OPEN(44,FILE='dualgrid.dat',FORM='FORMATTED')
DO j = 1, nedgex
  if1 = fnxte(j,1,ngrids)
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  CALL ll2xyz(long,lat,x1,y1,z1)
  if1 = fnxte(j,2,ngrids)
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  CALL ll2xyz(long,lat,x2,y2,z2)
  WRITE(44,'(6e15.7)') x1,y1,z1,x2,y2,z2
ENDDO
CLOSE(44)



! Output gridmap file
IF (flavour == 1) THEN
  WRITE(ygridfile,'(''gridmap_eqcu_'',I10.10,''.dat'')') nfacex
ELSEIF (flavour == 2) THEN
  WRITE(ygridfile,'(''gridmap_tcd_'',I10.10,''.dat'')') nfacex
ELSEIF (flavour == 3) THEN
  WRITE(ygridfile,'(''gridmap_ctcu_'',I10.10,''.dat'')') nfacex
ELSE
  PRINT *,'flavour = ',flavour,'  not recognized. Please pick another.'
  STOP
ENDIF


OPEN(22,FILE=ygridfile,FORM='UNFORMATTED')

! WRITE(22,*) 'GRIDMAP for NGRIDS=',NGRIDS
WRITE(22) ngrids
WRITE(22) nface
WRITE(22) nedge
WRITE(22) nvert
! WRITE(22,*) 'Number of edges of each face - all grids'
WRITE(22) ((neoff(if1,igrid),            &
               if1 = 1, nface(igrid)),   &
               igrid = 1, ngrids)
! WRITE(22,*) 'Number of edges of each vertex - all grids'
WRITE(22) ((neofv(iv1,igrid),            &
               iv1 = 1, nvert(igrid)),   &
               igrid=1, ngrids)
! WRITE(22,*) 'Faces next to each face - all grids'
WRITE(22) (((fnxtf(if1,if2,igrid),       &
               if1 = 1, nface(igrid)),   &
               if2 = 1, 4),              &
               igrid = 1, ngrids)
! WRITE(22,*) 'Edges of each face - all grids'
WRITE(22) (((eoff(if1,ie1,igrid),        &
               if1 = 1, nface(igrid)),   &
               ie1 = 1, 4),              &
               igrid = 1, ngrids)
! WRITE(22,*) 'Vertices of each face - all grids'
WRITE(22) (((voff(if1,iv1,igrid),        &
               if1 = 1, nface(igrid)),   &
               iv1 = 1, 4),              &
               igrid = 1, ngrids)
! WRITE(22,*) 'Faces next to each edge - all grids'
WRITE(22) (((fnxte(ie1,if2,igrid),       &
               ie1 = 1, nedge(igrid)),   &
               if2 = 1, 2),              &
               igrid = 1, ngrids)
! WRITE(22,*) 'Vertices of each edge - all grids'
WRITE(22) (((vofe(ie1,iv2,igrid),        &
               ie1 = 1, nedge(igrid)),   &
               iv2 = 1, 2),              &
               igrid = 1, ngrids)
! WRITE(22,*) 'Faces around each vertex - all grids'
WRITE(22) (((fofv(iv1,if2,igrid),        &
               iv1 = 1, nvert(igrid)),   &
               if2 = 1, 4),              &
               igrid = 1, ngrids)
! WRITE(22,*) 'Edges around each vertex - all grids'
WRITE(22) (((eofv(iv1,IE1,igrid),        &
               iv1 = 1, nvert(igrid)),   &
               ie1 = 1, 4),              &
               igrid = 1, ngrids)
! WRITE(22,*) 'Longitudes of faces - all grids'
WRITE(22) ((flong(if1,igrid),            &
               if1 = 1, nface(igrid)),   &
               igrid = 1, ngrids)
! WRITE(22,*) 'Latitudes of faces - all grids'
WRITE(22) ((flat(if1,igrid),             &
               if1 = 1, nface(igrid)),   &
               igrid = 1, ngrids)
! WRITE(22,*) 'Longitudes of vertices - all grids'
WRITE(22) ((vlong(iv1,igrid),            &
               iv1 = 1, nvert(igrid)),   &
               igrid = 1, ngrids)
! WRITE(22,*) 'Latitudes of vertices - all grids'
WRITE(22) ((vlat(iv1,igrid),             &
               iv1 = 1, nvert(igrid)),   &
               igrid = 1, ngrids)
! WRITE(22,*) 'Areas of faces - all grids'
WRITE(22) ((farea(if1,igrid),            &
              if1 = 1, nface(igrid)),    &
              igrid = 1, ngrids)
! WRITE(22,*) 'Lengths of edges - all grids'
WRITE(22) ((ldist(ie1,igrid),            &
              ie1 = 1, nedge(igrid)),    &
              igrid = 1, ngrids)
! WRITE(22,*) 'Distance between faces across edges - all grids'
WRITE(22) ((ddist(ie1,igrid),            &
              ie1 = 1, nedge(igrid)),    &
              igrid = 1, ngrids)

end subroutine write_grid

! --------------------------------------------------------------------------------------------------

subroutine read_grid(grid)

type(fempsgrid), intent(inout) :: grid


end subroutine read_grid

! --------------------------------------------------------------------------------------------------

subroutine write_pbops(pbops)

type(fempspbops), intent(in) :: pbops

IMPLICIT NONE

INTEGER :: if0, ie0, iv0, igrid, ix, ixx
CHARACTER*10 :: yface
CHARACTER*31 :: ygridfile

! Construct filename
WRITE(yface,'(I10.10)') nfacex
ygridfile = 'gridopermap_'//ygtype//'_'//yface//'.dat'

! Open file for writing
OPEN(changout,FILE=ygridfile,FORM='UNFORMATTED')

! First write ngrids
WRITE(changout) ngrids

! Write numbers of faces, edges and vertices on each grid
WRITE(changout) nface
WRITE(changout) nedge
WRITE(changout) nvert

! Write the numbers of edges of each face and vertex on each grid
WRITE(changout) ((neoff(if0,igrid),         &
                    if0 = 1, nface(igrid)), &
                    igrid = 1, ngrids)
WRITE(changout) ((neofv(iv0,igrid),         &
                    iv0 = 1, nvert(igrid)), &
                    igrid = 1, ngrids)

! Write the connectivity arrays
WRITE(changout) (((fnxtf(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((eoff(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((voff(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((fnxte(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
WRITE(changout) (((vofe(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
WRITE(changout) (((fofv(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((eofv(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)

! Write the geometrical information arrays
WRITE(changout) ((flong(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((flat(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((vlong(iv0,igrid),             &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((vlat(iv0,igrid),              &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((farea(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((varea(iv0,igrid),             &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ldist(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ddist(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)

! Write the sizes of the operator stencils on each grid
WRITE(changout) ((nlsten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nmsten(ie0,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((njsten(iv0,igrid),            &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nhsten(ie0,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nrsten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nrxsten(iv0,igrid),           &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((nwsten(ie0,igrid),            &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
WRITE(changout) ((ntsten(if0,igrid),            &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)

! Write the operator stencils and coefficients
WRITE(changout) (((lsten(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nlsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((msten(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nmsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((jsten(iv0,ix,igrid),         &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((hsten(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rsten(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rxsten(iv0,ix,igrid),        &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nrxsmx),           &
                     igrid = 1, ngrids)
WRITE(changout) (((wsten(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nwsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((tsten(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ntsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((lmass(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nlsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((mmass(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nmsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((jstar(iv0,ix,igrid),         &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((hstar(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rcoeff(if0,ix,igrid),        &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) (((rxcoeff(iv0,ix,igrid),       &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nrxsmx),           &
                     igrid = 1, ngrids)
WRITE(changout) (((wcoeff(ie0,ix,igrid),        &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nwsmx),            &
                     igrid = 1, ngrids)
WRITE(changout) ((((tcoeff(if0,ix,ixx,igrid),   &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ntsmx),            &
                     ixx = 1, ntsmx),           &
                     igrid = 1, ngrids)

! Write the size of the restriction stencil
WRITE(changout) ((ninj(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids-1)

! Write the restriction stencil and coefficients
WRITE(changout) (((injsten(if0,ix,igrid),       &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ninjmx),           &
                     igrid = 1, ngrids-1)
WRITE(changout) (((injwgt(if0,ix,igrid),        &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ninjmx),           &
                     igrid = 1, ngrids-1)

end subroutine write_pbops

! --------------------------------------------------------------------------------------------------

end module femps_inout_mod
