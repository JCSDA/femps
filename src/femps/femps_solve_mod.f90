module femps_solve_mod

use femps_const_mod
use femps_utils_mod
use femps_types_mod
use femps_kinds_mod

implicit none
private

public preliminary

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine preliminary(grid,pbops)

! Preliminary calculations and setting up

implicit none
type(fempsgrid),  intent(inout) :: grid
type(fempspbops), intent(inout) :: pbops

! Read namelist information
! -------------------------
call readnml(grid)
print *,'Namelists read'


! Read in the grid data
! ---------------------
call readgrid(grid,pbops)
print *,'Done readgrid'


! Build a lumped version of the J, M and H matrices
! -------------------------------------------------
call buildjlump(grid,pbops)
call buildmlump(grid,pbops)
call buildhlump(grid,pbops)


! And build a local operator to approximate the inverse of M
! ----------------------------------------------------------
call buildxminv(grid,pbops)

print *,'Lumped J, M and H and approximate M-inverse created'

end subroutine preliminary

! --------------------------------------------------------------------------------------------------

subroutine readnml(grid)

implicit none

type(fempsgrid),  intent(inout) :: grid
integer, parameter :: channml = 20
character*31 :: ygridfile

! Read namelists
! --------------
NAMELIST /rundata/ ygridfile

OPEN(channml,FILE='poissonnml.in',DELIM='APOSTROPHE')
READ(channml,rundata)
CLOSE(channml)

grid%ygridfile = ygridfile

end subroutine readnml

! --------------------------------------------------------------------------------------------------

subroutine buildlap(grid,pbops)

! Extract the diagonal coefficients of the Laplacian operator needed
! for the multigrid Poisson solver at all resolutions

implicit none

type(fempsgrid),  intent(in) :: grid
type(fempspbops), intent(inout) :: pbops

integer :: igrid, if1, if2, if3, ie1, ie2, &
           ixd2, ixm, ixd1, ixl
real*8 :: temp, temp1(grid%nfacex), cd2, cm, cd1, cl


! Under-relaxation parameter. ( *** There is scope to optimize here *** )
! -----------------------------------------------------------------------
do igrid = 1, grid%ngrids
  IF (grid%nefmx == 6) THEN
    pbops%underrel(igrid) = 0.8d0
  ELSEIF (grid%nefmx == 4) THEN
    pbops%underrel(igrid) = 0.8d0
  ELSE
    print *,'Choose a sensible value for underrel in buildlap'
    STOP
  ENDIF
enddo


! Extract diagonal coefficient of Laplacian operator
do igrid = 1, grid%ngrids
  ! Loop over cells
  do if1 = 1, grid%nface(igrid)

    temp = 0.0d0
    ! Loop over edges of if1 involved in Dprimal2 operator
    do ixd2 = 1, grid%neoff(if1,igrid)
      ie1 = grid%eoff(if1,ixd2,igrid)
      cd2 = -pbops%eoffin(if1,ixd2,igrid)
      ! Loop over edges involved in approximate M-inverse operator
      do ixm = 1, pbops%nxminvsten(ie1,igrid)
        ie2 = pbops%xminvsten(ie1,ixm,igrid)
        cm = pbops%xminv(ie1,ixm,igrid)
        ! Loop over cells involved in Ddual1 operator
        cd1 = 1.0d0
        do ixd1 = 1, 2
          if2 = grid%fnxte(ie2,ixd1,igrid)
          cd1 = -cd1
          ! Loop over cells in L operator
          do ixl = 1, pbops%nlsten(if2,igrid)
            if3 = pbops%lsten(if2,ixl,igrid)
            cl = pbops%lmass(if2,ixl,igrid)
            IF (if3 == if1) temp = temp + cd2*cm*cd1*cl
          enddo
        enddo
      enddo
    enddo
    pbops%lapdiag(if1,igrid) = temp

  enddo
enddo

end subroutine buildlap

! --------------------------------------------------------------------------------------------------

subroutine Dprimal1(grid,f,df,igrid,nv,ne)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises pointwise values
! at vertices; df comprises integrals of the derivative
! of f along primal cell edges.

implicit none
type(fempsgrid),  intent(in) :: grid
integer,      intent(in)  :: igrid, nv, ne
real*8,       intent(in)  :: f(nv)
real*8,       intent(out) :: df(ne)

integer :: ie1, iv1, iv2

do ie1 = 1, ne
  iv1 = grid%vofe(ie1,1,igrid)
  iv2 = grid%vofe(ie1,2,igrid)
  df(ie1) = f(iv2) - f(iv1)
enddo

end subroutine Dprimal1

! --------------------------------------------------------------------------------------------------

subroutine Dprimal2(grid,pbops,f,df,igrid,ne,nf)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises integrals along
! primal edges; df comprises integrals of the derivative
! over primal cells.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops

integer,      intent(in)  :: igrid, ne, nf
real*8,       intent(in)  :: f(ne)
real*8,       intent(out) :: df(nf)

integer :: if1, ix, ie1
real*8 :: temp

do if1 = 1, nf
  temp = 0.0d0
  do ix = 1, grid%neoff(if1,igrid)
    ie1 = grid%eoff(if1,ix,igrid)
    temp = temp - f(ie1)*pbops%eoffin(if1,ix,igrid)
  enddo
  df(if1) = temp
enddo

end subroutine Dprimal2

! --------------------------------------------------------------------------------------------------

subroutine Ddual1(grid,f,df,igrid,nf,ne)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises pointwise values
! at face centres; df comprises integrals of the derivative
! of f along dual cell edges.

implicit none
type(fempsgrid), intent(in) :: grid
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf)
real*8,       intent(out) :: df(ne)

integer :: ie1, if1, if2

do ie1 = 1, ne
  if1 = grid%fnxte(ie1,1,igrid)
  if2 = grid%fnxte(ie1,2,igrid)
  df(ie1) = f(if2) - f(if1)
enddo

end subroutine Ddual1

! --------------------------------------------------------------------------------------------------

subroutine Ddual2(grid,pbops,f,df,igrid,ne,nv)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises integrals along
! dual edges; df comprises integrals of the derivative
! over dual cells.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne, nv
real*8,       intent(in)  :: f(ne)
real*8,       intent(out) :: df(nv)

integer :: iv1, ix, ie1
real*8 :: temp

do iv1 = 1, nv
  temp = 0.0d0
  do ix = 1, grid%neofv(iv1,igrid)
    ie1 = grid%eofv(iv1,ix,igrid)
    temp = temp + f(ie1)*pbops%eofvin(iv1,ix,igrid)
  enddo
  df(iv1) = temp
enddo

end subroutine Ddual2

! --------------------------------------------------------------------------------------------------

subroutine buildjlump(grid,pbops)

! Build a lumped version of the j matrix

implicit none

type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(inout) :: pbops

! Two choices for lumped J
! ijlump = 1     Diagonal of J
! ijlump = 2     Exact answer when applied to constant input vector
! ijlump = 3     Exact answer when applied to the discrete representation
!                  of a constant scalar (input vector proportional to varea)

integer, parameter :: ijlump = 3
integer :: igrid, iv1, nv
real*8, allocatable :: p1(:), jp1(:)

do igrid = 1, grid%ngrids

  if (ijlump == 1) then

    do iv1 = 1, grid%nvert(igrid)
      pbops%jlump(iv1,igrid) = pbops%jstar(iv1,1,igrid)
    enddo

  elseif (ijlump == 2) then

    nv = grid%nvert(igrid)
    allocate(p1(nv), jp1(nv))
    p1 = 1.0d0
    call HodgeJ(pbops,p1,jp1,igrid,nv)
    pbops%jlump(1:nv,igrid) = jp1
    deallocate(p1,jp1)

  elseif (ijlump == 3) then

    nv = grid%nvert(igrid)
    allocate(p1(nv), jp1(nv))
    p1 = pbops%varea(1:nv,igrid)
    call HodgeJ(pbops,p1,jp1,igrid,nv)
    pbops%jlump(1:nv,igrid) = jp1/pbops%varea(1:nv,igrid)
    deallocate(p1,jp1)

  else

    print *,'option ijlump = ',ijlump,' not available in subroutine buildjlump'
    stop

  endif

enddo


end subroutine buildjlump

! --------------------------------------------------------------------------------------------------

subroutine buildmlump(grid,pbops)

! Build a lumped version of the M matrix

implicit none

type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(inout) :: pbops

! Three choices for lumped M
! imlump = 1     Diagonal of M
! imlump = 2     Row sum of absolute values of M
! imlump = 3     Best fit to solid body rotation normal to each edge

integer, parameter :: imlump = 3
integer :: igrid, ie1, ie2, ix, iv0, ne, nv, iv1, iv2
real*8 :: temp, slon, slat, clon, clat, a1, a2, a3, num, den, &
          b1, b2, b3, c1, c2, c3, x, y, z
real*8, allocatable :: psi1(:), psi2(:), psi3(:),    &
                       v1(:), v2(:), v3(:), vv(:),   &
                       mv1(:), mv2(:), mv3(:)

do igrid = 1, grid%ngrids

  if (imlump == 1) then

    do ie1 = 1, grid%nedge(igrid)
      pbops%mlump(ie1,igrid) = pbops%mmass(ie1,1,igrid)
    enddo

  elseif (imlump == 2) then

    do ie1 = 1, grid%nedge(igrid)
      temp = 0.0d0
      do ix = 1, pbops%nmsten(ie1,igrid)
        ie2 = pbops%msten(ie1,ix,igrid)
        temp = temp + abs(pbops%mmass(ie1,ix,igrid))
      enddo
      pbops%mlump(ie1,igrid) = temp
    enddo

  elseif (imlump == 3) then

    nv = grid%nvert(igrid)
    ne = grid%nedge(igrid)
    allocate(psi1(nv), psi2(nv), psi3(nv),    &
             v1(ne), v2(ne), v3(ne), vv(ne),  &
             mv1(ne), mv2(ne), mv3(ne))
    ! Set up three solid body rotation fields
    do 	iv0 = 1, nv
      slon = sin(grid%vlong(iv0,igrid))
      clon = cos(grid%vlong(iv0,igrid))
      slat = sin(grid%vlat(iv0,igrid))
      clat = cos(grid%vlat(iv0,igrid))
      psi1(iv0) = clat*clon
      psi2(iv0) = clat*slon
      psi3(iv0) = slat
    enddo
    call Dprimal1(grid,psi1,v1,igrid,nv,ne)
    call Dprimal1(grid,psi2,v2,igrid,nv,ne)
    call Dprimal1(grid,psi3,v3,igrid,nv,ne)
    call massM(pbops,v1,mv1,igrid,ne)
    call massM(pbops,v2,mv2,igrid,ne)
    call massM(pbops,v3,mv3,igrid,ne)
    ! Now loop over edges
    do ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      den = sqrt(a1*a1 + a2*a2 + a3*a3)
      a1 = a1/den
      a2 = a2/den
      a3 = a3/den
      ! Demand that lumped matrix agrees with full M for this v
      num = a1*mv1(ie1) + a2*mv2(ie1) + a3*mv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      pbops%mlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,vv,mv1,mv2,mv3)

  else

    print *,'option imlump = ',imlump,' not available in subroutine buildmlump'
    stop

  endif

enddo

end subroutine buildmlump

! --------------------------------------------------------------------------------------------------

subroutine buildhlump(grid,pbops)

! Build a lumped version of the H matrix

implicit none

type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(inout) :: pbops

! Three choices for lumped H
! ihlump = 1     Diagonal of H
! ihlump = 2     Row sum of absolute values of H
! ihlump = 3     Best fit to global divergent flow on dual grid
! ihlump = 4     Best fit to solid body rotation normal to each edge

integer, parameter :: ihlump = 3
integer :: igrid, ie1, ie2, ix, iv0, if0, ne, nv, nf
real*8 :: temp, slon, slat, clon, clat, a1, a2, a3, num, den
real*8, allocatable :: psi1(:), psi2(:), psi3(:),  &
                       v1(:), v2(:), v3(:),        &
                       hv1(:), hv2(:), hv3(:)

do igrid = 1, grid%ngrids

  if (ihlump == 1) then

    do ie1 = 1, grid%nedge(igrid)
      pbops%hlump(ie1,igrid) = pbops%hstar(ie1,1,igrid)
    enddo

  elseif (ihlump == 2) then

    do ie1 = 1, grid%nedge(igrid)
      temp = 0.0d0
      do ix = 1, pbops%nhsten(ie1,igrid)
        ie2 = pbops%hsten(ie1,ix,igrid)
        temp = temp + abs(pbops%hstar(ie1,ix,igrid))
      enddo
      pbops%hlump(ie1,igrid) = temp
    enddo

  elseif (ihlump == 3) then

    nf = grid%nface(igrid)
    ne = grid%nedge(igrid)
    allocate(psi1(nf), psi2(nf), psi3(nf),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three global divergent fields
    do 	if0 = 1, nf
      slon = sin(grid%flong(if0,igrid))
      clon = cos(grid%flong(if0,igrid))
      slat = sin(grid%flat(if0,igrid))
      clat = cos(grid%flat(if0,igrid))
      psi1(if0) = slat
      psi2(if0) = clat*slon
      psi3(if0) = clat*clon
    enddo
    call Ddual1(grid,psi1,v1,igrid,nf,ne)
    call Ddual1(grid,psi2,v2,igrid,nf,ne)
    call Ddual1(grid,psi3,v3,igrid,nf,ne)
    call HodgeH(pbops,v1,hv1,igrid,ne)
    call HodgeH(pbops,v2,hv2,igrid,ne)
    call HodgeH(pbops,v3,hv3,igrid,ne)
    ! Now loop over edges
    do ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      ! Demand that lumped matrix agrees with full H for this v
      num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      pbops%hlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  elseif (ihlump == 4) then

    nv = grid%nvert(igrid)
    ne = grid%nedge(igrid)
    allocate(psi1(nv), psi2(nv), psi3(nv),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three solid body rotation fields
    do 	iv0 = 1, nv
      slon = sin(grid%vlong(iv0,igrid))
      clon = cos(grid%vlong(iv0,igrid))
      slat = sin(grid%vlat(iv0,igrid))
      clat = cos(grid%vlat(iv0,igrid))
      psi1(iv0) = slat
      psi2(iv0) = clat*slon
      psi3(iv0) = clat*clon
    enddo
    call Dprimal1(grid,psi1,v1,igrid,nv,ne)
    call Dprimal1(grid,psi2,v2,igrid,nv,ne)
    call Dprimal1(grid,psi3,v3,igrid,nv,ne)
    call HodgeH(pbops,v1,hv1,igrid,ne)
    call HodgeH(pbops,v2,hv2,igrid,ne)
    call HodgeH(pbops,v3,hv3,igrid,ne)
    ! Now loop over edges
    do ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      ! Demand that lumped matrix agrees with full H for this v
      num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      pbops%hlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  else

    print *,'option ihlump = ',ihlump,' not available in subroutine buildhlump'
    stop

  endif

enddo

end subroutine buildhlump

! --------------------------------------------------------------------------------------------------

subroutine buildxminv(grid,pbops)

! Determine stencil and coefficients for a local approximation
! to the inverse of M.
!
! Two options are coded. The first is simply the inverse of
! the diagonal mass lumped approximation to M. The second is
! based on a single underrelaxed Jacobi iteration to the
! inverse of M.

implicit none

type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(inout) :: pbops

logical :: llump
integer :: igrid, ie0, ix, ie1
real*8 :: temp, diag
real*8 :: relax

! Underrelaxation coefficient depends on grid
! 0.9 and 1.4 when using mlump in Jacobi
! Check for consistency with massMinv
if (grid%nefmx == 4) then
  ! use sparse approximate inverse on cubed sphere
  llump = .false.
  relax = 0.9d0
else
  ! Use diagonal approximate inverse on hexagonal-icosahedral grid
  llump = .false.
  relax = 1.4d0
endif

if (llump) then

  ! Diagonal operator: stencil size is 1
  pbops%nxmisx = 1
  allocate(pbops%nxminvsten(grid%nedgex,grid%ngrids), pbops%xminvsten(grid%nedgex,pbops%nxmisx,grid%ngrids), &
           pbops%xminv(grid%nedgex,pbops%nxmisx,grid%ngrids))
  do igrid = 1, grid%ngrids
    do ie0 = 1, grid%nedge(igrid)
      ! stencil for edge ie0 is ie0 itself, and coeff is
      ! inverse of the diagonal term of the lumped matrix
      pbops%nxminvsten(ie0,igrid) = 1
      pbops%xminvsten(ie0,1,igrid) = ie0
      pbops%xminv(ie0,1,igrid) = 1.0d0/pbops%mlump(ie0,igrid)
    enddo
  enddo

else

  ! Stencil is the same as for M itself
  pbops%nxmisx = pbops%nmsmx
  allocate(pbops%nxminvsten(grid%nedgex,grid%ngrids), pbops%xminvsten(grid%nedgex,pbops%nxmisx,grid%ngrids), &
           pbops%xminv(grid%nedgex,pbops%nxmisx,grid%ngrids))
  do igrid = 1, grid%ngrids
    do ie0 = 1, grid%nedge(igrid)
      ! stencil for edge ie0 is the same as the stencil for m
      pbops%nxminvsten(ie0,igrid) = pbops%nmsten(ie0,igrid)
      do ix = 1, pbops%nmsten(ie0,igrid)
        ie1 = pbops%msten(ie0,ix,igrid)
        pbops%xminvsten(ie0,ix,igrid) = ie1
        if (ie1 == ie0) then
          diag = 1.0d0 + relax
        else
          diag = 0.0d0
        endif
        temp = pbops%mmass(ie0,ix,igrid)/pbops%mlump(ie1,igrid)
        pbops%xminv(ie0,ix,igrid) = (diag - relax*temp)/pbops%mlump(ie0,igrid)
      enddo
    enddo
  enddo

endif

end subroutine buildxminv

! --------------------------------------------------------------------------------------------------

subroutine massL(pbops,f1,f2,igrid,nf)

! Apply the mass matrix L to field f1 to obtain the result f2

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, nf
real*8,       intent(in)  :: f1(nf)
real*8,       intent(out) :: f2(nf)
integer :: if1, if2, ix
real*8 :: temp

do if1 = 1, nf
  temp = 0.0d0
  do ix = 1, pbops%nlsten(if1,igrid)
    if2 = pbops%lsten(if1,ix,igrid)
    temp = temp + f1(if2)*pbops%lmass(if1,ix,igrid)
  enddo
  f2(if1) = temp
enddo

end subroutine massL

! --------------------------------------------------------------------------------------------------

subroutine massM(pbops,f1,f2,igrid,ne)

! Apply the mass matrix M to field f1 to obtain field f2

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)
integer :: ie1, ie2, ix
real*8 :: temp

do ie1 = 1, ne
  temp = 0.0d0
  do ix = 1, pbops%nmsten(ie1,igrid)
    ie2 = pbops%msten(ie1,ix,igrid)
    temp = temp + f1(ie2)*pbops%mmass(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine massM

! --------------------------------------------------------------------------------------------------

subroutine approxMinv(pbops,f1,f2,igrid,ne)

! Apply an approximate inverse of the mass matrix M
! to field f1 to obtain field f2

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: ie1, ie2, ix
real*8 :: temp

do ie1 = 1, ne
  temp = 0.0d0
  do ix = 1, pbops%nxminvsten(ie1,igrid)
    ie2 = pbops%xminvsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*pbops%xminv(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine approxMinv

! --------------------------------------------------------------------------------------------------

subroutine HodgeJ(pbops,f1,f2,igrid,nv)

! Apply the Hodge star J operator that converts dual face
! integrals f1 to vertex values f2 on grid igrid.

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in) :: igrid, nv
real*8,       intent(in) :: f1(nv)
real*8,       intent(out) :: f2(nv)

integer :: iv1, iv2, ix
real*8 :: temp

do iv1 = 1, nv
  temp = 0.0d0
  do ix = 1, pbops%njsten(iv1,igrid)
    iv2 = pbops%jsten(iv1,ix,igrid)
    temp = temp + f1(iv2)*pbops%jstar(iv1,ix,igrid)
  enddo
  f2(iv1) = temp
enddo

end subroutine HodgeJ

! --------------------------------------------------------------------------------------------------

subroutine HodgeH(pbops,f1,f2,igrid,ne)

! Apply the Hodge star H operator that converts dual edge
! integrals (circulations) f1 to primal edge integrals (fluxes) f2
! on grid igrid.

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: ie1, ie2, ix
real*8 :: temp

do ie1 = 1, ne
  temp = 0.0d0
  do ix = 1, pbops%nhsten(ie1,igrid)
    ie2 = pbops%hsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*pbops%hstar(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine HodgeH

! --------------------------------------------------------------------------------------------------

subroutine massLinv(pbops,f1,f2,igrid,nf,niter)

! Apply the inverse of the mass matrix L to field f1 to obtain
! the result f2.

! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If L is diagonal then there is no need to iterate.

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, nf
real*8,       intent(in)    :: f1(nf)
real*8,       intent(inout) :: f2(nf)

integer :: if1, iter, miter
real*8 :: temp(nf)

miter = ABS(niter)

IF (niter < 0 .OR. pbops%nlsmx == 1) THEN
  ! First guess based on diagonal L
  do if1 = 1, nf
    f2(if1) = f1(if1)/pbops%lmass(if1,1,igrid)
  enddo
ENDIF

IF (pbops%nlsmx > 1) THEN
  ! L is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call massL(pbops,f2,temp,igrid,nf)
    temp = f1 - temp
    do if1 = 1, nf
      f2(if1) = f2(if1) + temp(if1)/pbops%lmass(if1,1,igrid)
    enddo
  enddo
ENDIF

end subroutine massLinv

! --------------------------------------------------------------------------------------------------

subroutine massMinv(grid,pbops,f1,f2,igrid,ne,niter)

! Apply the inverse of the mass matrix M to the field f1
! to obtain the result f2

! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If M is diagonal then there is no need to iterate.

implicit none

type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, ne
real*8,       intent(in)    :: f1(ne)
real*8,       intent(inout) :: f2(ne)

integer :: ie1, iter, miter
real*8 :: temp(ne)
real*8 :: relax

! Underrelaxation coefficient depends on grid
IF (grid%nefmx == 4) THEN
  relax = 0.9d0
ELSE
  relax = 1.4d0
ENDIF

miter = ABS(niter)

IF (niter < 0 .OR. pbops%nmsmx == 1) THEN
  ! First guess based on lumped M
  do ie1 = 1, ne
    f2(ie1) = f1(ie1)/pbops%mlump(ie1,igrid)
  enddo
ENDIF

IF (pbops%nmsmx > 1) THEN
  ! M is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call massM(pbops,f2,temp,igrid,ne)
    temp = relax*(f1 - temp)
    do ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/pbops%mlump(ie1,igrid)
    enddo
  enddo
ENDIF

end subroutine massMinv

! --------------------------------------------------------------------------------------------------

subroutine HodgeJinv(pbops,f1,f2,igrid,nv,niter)

! Apply the inverse Hodge star operator J^{-1} that maps from
! E_p to V_d on grid igrid.
!
! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If J is diagonal then there is no need to iterate.

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, nv
real*8,       intent(in)    :: f1(nv)
real*8,       intent(inout) :: f2(nv)

integer :: iv1, iter, miter
real*8 :: temp(nv)
real*8 :: relax = 1.4d0 ! relax = 1.4 is good for ijlump = 3 on hex and cube grids

miter = ABS(niter)

IF (niter < 0 .OR. pbops%njsmx == 1) THEN
  ! First guess based on lumped J
  do iv1 = 1, nv
    f2(iv1) = f1(iv1)/pbops%jlump(iv1,igrid)
  enddo
ENDIF

IF (pbops%njsmx > 1) THEN
  ! J is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call HodgeJ(pbops,f2,temp,igrid,nv)
    temp = relax*(f1 - temp)
    do iv1 = 1, nv
      f2(iv1) = f2(iv1) + temp(iv1)/pbops%jlump(iv1,igrid)
    enddo
  enddo
ENDIF

end subroutine HodgeJinv

! --------------------------------------------------------------------------------------------------

subroutine HodgeHinv(pbops,f1,f2,igrid,ne,niter)

! Apply the inverse Hodge star operator H^{-1} that maps from
! S_p to S_d on grid igrid.
!
! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If H is diagonal then there is no need to iterate.

implicit none
type(fempspbops), intent(in) :: pbops
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, ne
real*8,       intent(in)    :: f1(ne)
real*8,       intent(inout) :: f2(ne)

integer :: ie1, iter, miter
real*8 :: temp(ne)
real*8 :: relax = 1.4d0 ! relax = 1.4 is good for ihlump = 3 on hex and cube grids

miter = ABS(niter)

IF (niter < 0 .OR. pbops%nhsmx == 1) THEN
  ! First guess based on diagonal H
  do ie1 = 1, ne
    f2(ie1) = f1(ie1)/pbops%hlump(ie1,igrid)
  enddo
ENDIF

IF (pbops%nhsmx > 1) THEN
  ! H is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call HodgeH(pbops,f2,temp,igrid,ne)
    temp = relax*(f1 - temp)
    do ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/pbops%hlump(ie1,igrid)
    enddo
  enddo
ENDIF

end subroutine HodgeHinv

! --------------------------------------------------------------------------------------------------

subroutine operW_original(grid,pbops,f1,f2,igrid,ne)

! Apply the W operator:
! given fluxes f1 across primal edges, construct
! the rotated fluxes across dual edges f2, on grid igrid.

! This is the original formulation, building W from R
! a la TRiSK. It probably requires an MPI reduce operation
! so is likely to be inefficient.

implicit none

type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: if1, ie1, ie2, ix, ix1, ix2, ixv, ne1
real*8 :: ss, w

! Initialize to zero
f2 = 0.0d0

! Loop over faces
do if1 = 1, grid%nface(igrid)
  ne1 = grid%neoff(if1,igrid)
  ! For each edge of this face
  do ix1 = 1, ne1
    ss = -0.5
    ie1 = grid%eoff(if1,ix1,igrid)
    ! Find the contribution to f2 from every other
    ! edge of this face
    do ix = 0, ne1 - 2
      ixv = MOD(ix1 + ix - 1,ne1) + 1
      ix2 = MOD(ix1 + ix,ne1) + 1
      ie2 = grid%eoff(if1,ix2,igrid)
      ss = ss + pbops%rcoeff(if1,ixv,igrid)
      w = -ss*pbops%eoffin(if1,ix1,igrid)*pbops%eoffin(if1,ix2,igrid)
      f2(ie1) = f2(ie1) + w*f1(ie2)
    enddo
  enddo
enddo

end subroutine operW_original

! --------------------------------------------------------------------------------------------------

subroutine operR_original(grid,pbops,f1,f2,igrid,nf,nv)

! Apply the R operator:
! map from V_p to E_p

! This is the original formulation. It loops over `source'
! entities rather than target entities and so will require
! an MPI reduce.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, nf, nv
real*8,       intent(in)  :: f1(nf)
real*8,       intent(out) :: f2(nv)

integer :: if1, iv1, ix1, ne1

! Initialize to zero
f2 = 0.0d0

! Loop over faces
do if1 = 1, grid%nface(igrid)
  ne1 = grid%neoff(if1,igrid)
  ! Share out this face's contributions to its surrounding vertices
  do ix1 = 1, ne1
    iv1 = pbops%rsten(if1,ix1,igrid)
    f2(iv1) = f2(iv1) + f1(if1)*pbops%rcoeff(if1,ix1,igrid)
  enddo
enddo

end subroutine operR_original

! --------------------------------------------------------------------------------------------------

subroutine operW(grid,pbops,f1,f2,igrid,ne)

! Apply the W operator:
! given edge integrals of normal components f1 on primal edges,
! construct edge integrals of normal components of perpendicular
! field f2, on grid igrid.

! This formulation uses pre-build stencil and coefficients to
! avoid the need for MPI reduce. It is mathematically equivalent
! to the original formulation.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: ie0, ne1, ix1, ie1
real*8 :: temp

! Loop over vertices
do ie0 = 1, grid%nedge(igrid)
  ne1 = pbops%nwsten(ie0,igrid)
  ! Collect contributions from stencil
  temp = 0.0d0
  do ix1 = 1, ne1
    ie1 = pbops%wsten(ie0,ix1,igrid)
    temp = temp + f1(ie1)*pbops%wcoeff(ie0,ix1,igrid)
  enddo
  f2(ie0) = temp
enddo

end subroutine operW

! --------------------------------------------------------------------------------------------------

subroutine operR(grid,pbops,f1,f2,igrid,nf,nv)

! Apply the R operator:
! given face integrals f1 on primal faces, map to dual cell
! integrals f2, on grid igrid.

! This formulation stores the coefficients in the transpose of
! the original formulation to avoid an MPI reduce. It is
! mathematically equivalent to the original formulation.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, nf, nv
real*8,       intent(in)  :: f1(nf)
real*8,       intent(out) :: f2(nv)

integer :: iv0, if1, ix1, ne1
real*8 :: temp

! Loop over vertices
do iv0 = 1, grid%nvert(igrid)
  ne1 = pbops%nrxsten(iv0,igrid)
  ! Collect contributions from surrounding faces
  temp = 0.0d0
  do ix1 = 1, ne1
    if1 = pbops%rxsten(iv0,ix1,igrid)
    temp = temp + f1(if1)*pbops%rxcoeff(iv0,ix1,igrid)
  enddo
  f2(iv0) = temp
enddo

end subroutine operR

! --------------------------------------------------------------------------------------------------

subroutine operT(grid,pbops,f1,f2,igrid,ne,nf)

! Apply the T operator:
! compute cell integrals of 2 x kinetic energy from edge integrals
! of normal fluxes

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne, nf
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(nf)

integer :: if1, ix1, ix2, ne1, ie1, ie2
real*8 :: temp

! Loop over faces
do if1 = 1, grid%nface(igrid)
  ne1 = pbops%ntsten(if1,igrid)
  temp = 0.0d0
  ! Loop over all pairs of edges of this cell
  do ix1 = 1, ne1
    ie1 = pbops%tsten(if1,ix1,igrid)
    do ix2 = 1, ne1
      ie2 = pbops%tsten(if1,ix2,igrid)
      temp = temp + pbops%tcoeff(if1,ix1,ix2,igrid)*f1(ie1)*f1(ie2)
    enddo
  enddo
  f2(if1) = temp
enddo

end subroutine operT

! --------------------------------------------------------------------------------------------------

subroutine restrict(grid,pbops,f1,nf1,f2,nf2,igrid)

! To perform the restriction operation needed for a multigrid solver.
! Restrict field f1 from grid igrid + 1 to grid igrid and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms).

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: nf1, nf2, igrid
real*8,       intent(in)  :: f1(nf1)
real*8,       intent(out) :: f2(nf2)

integer :: if1, if2, ix
real*8 :: wgt

! Safety check
IF (nf2 .ne. grid%nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine restrict'
  STOP
ENDIF

do if2 = 1, nf2
  f2(if2) = 0.0d0
  do ix = 1, pbops%ninj(if2,igrid)
    if1 = pbops%injsten(if2,ix,igrid)
    wgt = pbops%injwgt(if2,ix,igrid)
    f2(if2) = f2(if2) + wgt*f1(if1)
  enddo
enddo

end subroutine restrict

! --------------------------------------------------------------------------------------------------

subroutine prolong(grid,pbops,f2,nf2,f1,nf1,igrid)

! To perform the prolongation operation needed for a multigrid solver.
! Prolong field f2 from grid igrid to grid igrid + 1 and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms); so f2 must be converted to point values
! (zero-form) first and f1 must be converted back to an area-integral
! (2-form) at the end. The prolongation operator is the adjoint of
! the restriction operator, so uses the same stencil and weights.

! *** Currently this operator uses a loop over source entities,
! implying a need for an MPI reduce. The stencil and weights could
! be stored in transpose form to allow a loop over target entities
! and hence avoid an MPI reduce ***

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: nf1, nf2, igrid
real*8,       intent(in)  :: f2(nf2)
real*8,       intent(out) :: f1(nf1)
integer :: if1, if2, ix, igridp
real*8 :: wgt, f2if2, temp1(nf1), temp2(nf2)

! Safety check
IF (nf2 .ne. grid%nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine prolong'
  STOP
ENDIF

igridp = igrid + 1
temp2(1:nf2) = f2(1:nf2)/grid%farea(1:nf2,igrid)
temp1 = 0.0d0
do if2 = 1, nf2
  f2if2 = temp2(if2)
  do ix = 1, pbops%ninj(if2,igrid)
    if1 = pbops%injsten(if2,ix,igrid)
    wgt = pbops%injwgt(if2,ix,igrid)
    temp1(if1) = temp1(if1) + wgt*f2if2
  enddo
enddo
f1(1:nf1) = temp1(1:nf1)*grid%farea(1:nf1,igridp)

end subroutine prolong

! --------------------------------------------------------------------------------------------------

subroutine laplace(grid,pbops,f,hf,igrid,nf,ne)

! To apply the Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf)
real*8,       intent(out) :: hf(nf)

real*8 :: temp1(nf), temp2(ne), temp3(ne)
integer :: niter

call massL(pbops,f,temp1,igrid,nf)
call Ddual1(grid,temp1,temp2,igrid,nf,ne)
niter = -20
call massMinv(grid,pbops,temp2,temp3,igrid,ne,niter)
call Dprimal2(grid,pbops,temp3,hf,igrid,ne,nf)


end subroutine laplace

! --------------------------------------------------------------------------------------------------

subroutine xlaplace(grid,pbops,f,hf,igrid,nf,ne)

! To apply the APPROXIMATE Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf)
real*8,       intent(out) :: hf(nf)

real*8 :: temp1(nf), temp2(ne), temp3(ne)

call massL(pbops,f,temp1,igrid,nf)
call Ddual1(grid,temp1,temp2,igrid,nf,ne)
call approxMinv(pbops,temp2,temp3,igrid,ne)
call Dprimal2(grid,pbops,temp3,hf,igrid,ne,nf)

end subroutine xlaplace

! --------------------------------------------------------------------------------------------------

subroutine residual(grid,pbops,f,rhs,res,igrid,nf,ne)

! Compute the residual res in the approximate Poisson equation on grid igrid
! when f is the input field and rhs is the right hand side. Note that
! f, rhs and res are area integrals (2-forms).

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf), rhs(nf)
real*8,       intent(out) :: res(nf)

call xlaplace(grid,pbops,f,res,igrid,nf,ne)
res = rhs - res
!print *,'     residual: ',res(1:5)

end subroutine residual

! --------------------------------------------------------------------------------------------------

subroutine relax(grid,pbops,f,rhs,igrid,nf,ne,niter)

! To carry out niter Jacobi relaxation iterations for the multigrid
! solver on grid igrid

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)    :: igrid, nf, ne, niter
real*8,       intent(in)    :: rhs(nf)
real*8,       intent(inout) :: f(nf)

real*8, allocatable :: res(:), finc(:)
real*8 :: u
integer :: iter

allocate(res(nf), finc(nf))

u = pbops%underrel(igrid)
do iter = 1, niter
  call residual(grid,pbops,f,rhs,res,igrid,nf,ne)
  finc = res/pbops%lapdiag(1:nf,igrid)
  f = f + u*finc
enddo

deallocate(res, finc)

end subroutine relax

! --------------------------------------------------------------------------------------------------

subroutine fullmgsolve(grid,pbops,phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 xminv Ddual1 L phi = RHS
!
! using full multigrid algorithm.
! Coefficients are contained in module laplace.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: ng
real*8,       intent(in)  :: rr(grid%nfacex)
real*8,       intent(out) :: phi(grid%nfacex)

! Numbers of iterations on coarsest grid and other grids
integer, parameter :: niterc = 10, niter = 2, npass = 1
integer :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
real*8, allocatable :: ff(:,:), rf(:,:)
real*8 :: temp1(grid%nfacex)

! Allocate space on all grids
allocate(ff(grid%nfacex,grid%ngrids),rf(grid%nfacex,grid%ngrids))

! One pass should be enough. Warn user if npass is set to
! some other value for testing
IF (npass > 1) PRINT *,'mgsolve: npass = ',npass

! Initialize solution to zero
phi = 0.0d0

! For diagnostics
!nf1 = grid%nface(grid%ngrids)
!ne1 = grid%nedge(grid%ngrids)
!call residual(grid,pbops,phi,rr,temp1,grid%ngrids,nf1,ne1)
!print *,'Pass ',0,'  RMS residual = ',SQRT(SUM(temp1*temp1)/nf1)

do ipass = 1, npass

  ! Initialize rhs as residual using latest estimate
  IF (ipass == 1) THEN
    ! No need to do the calculation
    rf(:,grid%ngrids) = rr(:)
  ELSE
    nf1 = grid%nface(grid%ngrids)
    ne1 = grid%nedge(grid%ngrids)
    call residual(grid,pbops,phi,rr,rf(1,grid%ngrids),grid%ngrids,nf1,ne1)
  ENDIF

  ! Initialize correction to solution on all grids to zero
  ff = 0.0d0

  ! Restrict right hand side to each grid in the hierarchy
  do igrid = grid%ngrids-1, grid%ngrids-ng+1, -1
    igridp = igrid + 1
    nf1 = grid%nface(igridp)
    nf2 = grid%nface(igrid)
    call restrict(grid,pbops,rf(1,igridp),nf1,rf(1,igrid),nf2,igrid)
  enddo

  ! Iterate to convergence on coarsest grid
  igrid = grid%ngrids-ng+1
  nf1 = grid%nface(igrid)
  ne1 = grid%nedge(igrid)
  ff(1:nf1,igrid) = 0.0d0
  call relax(grid,pbops,ff(1,igrid),rf(1,igrid),igrid,nf1,ne1,niterc)

  ! Sequence of growing V-cycles
  do igridp = grid%ngrids-ng+2, grid%ngrids

    igrid = igridp - 1
    nf1 = grid%nface(igridp)
    ne1 = grid%nedge(igridp)
    nf2 = grid%nface(igrid)
    ne2 = grid%nedge(igrid)

    ! Prolong solution to grid igridp
    ! and execute one V-cycle starting from grid igridp

    ! Prolong
    call prolong(grid,pbops,ff(1,igrid),nf2,ff(1,igridp),nf1,igrid)

    ! Descending part of V-cycle
    do jgrid = igrid, grid%ngrids-ng+1, -1

      jgridp = jgrid + 1
      nf1 = grid%nface(jgridp)
      ne1 = grid%nedge(jgridp)
      nf2 = grid%nface(jgrid)
      ne2 = grid%nedge(jgrid)

      ! Relax on grid jgridp
      call relax(grid,pbops,ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

      ! Calculate residual on jgridp
      call residual(grid,pbops,ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)

      ! Restrict residual to jgrid
      call restrict(grid,pbops,temp1,nf1,rf(1,jgrid),nf2,jgrid)

      ! Set correction first guess to zero on grid jgrid-1
      ff(1:nf2,jgrid) = 0.0d0

    enddo

    ! Iterate to convergence on coarsest grid
    jgrid = grid%ngrids-ng+1
    nf1 = grid%nface(jgrid)
    ne1 = grid%nedge(jgrid)
    ff(1:nf1,jgrid) = 0.0d0
    call relax(grid,pbops,ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

    ! Ascending part of V-cycle
    do jgrid = grid%ngrids-ng+1, igrid

      jgridp = jgrid + 1
      igrid = igrid - 1
      nf1 = grid%nface(jgridp)
      ne1 = grid%nedge(jgridp)
      nf2 = grid%nface(jgrid)
      ne2 = grid%nedge(jgrid)

      ! Prolong correction to jgridp
      call prolong(grid,pbops,ff(1,jgrid),nf2,temp1,nf1,jgrid)

      ! Add correction to solution on jgridp
      ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

      ! Relax on grid jgridp
      call relax(grid,pbops,ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

    enddo

  enddo

  ! Add correction to phi
  phi = phi + ff(:,grid%ngrids)

  ! For diagnostics
  !nf1 = grid%nface(grid%ngrids)
  !ne1 = grid%nedge(grid%ngrids)
  !call residual(grid,pbops,phi,rr,temp1,grid%ngrids,nf1,ne1)
  !print *,'      RMS residual in fullmgsolve = ',SQRT(SUM(temp1*temp1)/nf1)

enddo


deallocate(ff,rf)

end subroutine fullmgsolve

! --------------------------------------------------------------------------------------------------

subroutine mgsolve(grid,pbops,phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 xminv Ddual1 L phi = RHS
!
! using a single V-cycle multigrid algorithm.
! Coefficients are contained in module laplace.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: ng
real*8,       intent(in)  :: rr(grid%nfacex)
real*8,       intent(out) :: phi(grid%nfacex)

! Numbers of iterations on coarsest grid and other grids
integer, parameter :: niterc = 10, niter = 2
integer :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
real*8, allocatable :: ff(:,:), rf(:,:)
real*8 :: temp1(grid%nfacex)

! Allocate space on all grids
allocate(ff(grid%nfacex,grid%ngrids),rf(grid%nfacex,grid%ngrids))

! Initialize solution to zero
phi = 0.0d0

! Initialize rhs on finest grid
rf(:,grid%ngrids) = rr(:)

! Initialize correction to solution on all grids to zero
ff = 0.0d0


! Descending part of V-cycle
do jgrid = grid%ngrids-1, grid%ngrids-ng+1, -1

  jgridp = jgrid + 1
  nf1 = grid%nface(jgridp)
  ne1 = grid%nedge(jgridp)
  nf2 = grid%nface(jgrid)
  ne2 = grid%nedge(jgrid)

  ! Relax on grid jgridp
  call relax(grid,pbops,ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

  ! Calculate residual on jgridp
  call residual(grid,pbops,ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)

  ! Restrict residual to jgrid
  call restrict(grid,pbops,temp1,nf1,rf(1,jgrid),nf2,jgrid)

  ! Set correction first guess to zero on grid jgrid
  ff(1:nf2,jgrid) = 0.0d0

enddo

! Iterate to convergence on coarsest grid
jgrid = grid%ngrids-ng+1
nf1 = grid%nface(jgrid)
ne1 = grid%nedge(jgrid)
ff(1:nf1,jgrid) = 0.0d0
call relax(grid,pbops,ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

! Ascending part of V-cycle
do jgrid = grid%ngrids-ng+1, grid%ngrids-1

  jgridp = jgrid + 1
  nf1 = grid%nface(jgridp)
  ne1 = grid%nedge(jgridp)
  nf2 = grid%nface(jgrid)
  ne2 = grid%nedge(jgrid)

  ! Prolong correction to jgridp
  call prolong(grid,pbops,ff(1,jgrid),nf2,temp1,nf1,jgrid)

  ! Add correction to solution on jgridp
  ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

  ! Relax on grid jgridp
  call relax(grid,pbops,ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

enddo

! Add correction to phi
phi = phi + ff(:,grid%ngrids)

!  ! For diagnostics
!  nf1 = grid%nface(grid%ngrids)
!  ne1 = grid%nedge(grid%ngrids)
!  call residual(grid,pbops,phi,rr,temp1,grid%ngrids,nf1,ne1)
!  print *,'     RMS residual in mgsolve = ',SQRT(SUM(temp1*temp1)/nf1)

deallocate(ff,rf)

end subroutine mgsolve

! --------------------------------------------------------------------------------------------------

subroutine readgrid(grid,pbops)

! To allocate array space for the grid information in module grid
! and to read the information from file

implicit none
type(fempsgrid), intent(inout) :: grid
type(fempspbops), intent(inout) :: pbops

integer :: if0, ie0, iv0, igrid, ix, ixx, if1, if2, iv1, iv2, ie1
integer, parameter :: changrid = 25         ! grid information

! Open file for reading
OPEN(changrid,FILE=grid%ygridfile,FORM='UNFORMATTED')

! First read grid%ngrids
READ(changrid) grid%ngrids


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
         grid%farea(grid%nfacex,grid%ngrids), pbops%varea(grid%nvertx,grid%ngrids), &
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
READ(changrid) ((pbops%varea(iv0,igrid),              &
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
pbops%varea = pbops%varea*rearth*rearth
grid%ldist = grid%ldist*rearth
grid%ddist = grid%ddist*rearth

! Determine smallest face area on each grid
do igrid = 1, grid%ngrids
  grid%fareamin(igrid) = MINVAL(grid%farea(1:grid%nface(igrid),igrid))
enddo


! Allocate arrays for size of operator stencils
allocate(pbops%nlsten(grid%nfacex,grid%ngrids), pbops%nmsten(grid%nedgex,grid%ngrids), &
         pbops%njsten(grid%nvertx,grid%ngrids), pbops%nhsten(grid%nedgex,grid%ngrids), &
         pbops%nrsten(grid%nfacex,grid%ngrids), pbops%nrxsten(grid%nvertx,grid%ngrids), &
         pbops%nwsten(grid%nedgex,grid%ngrids), pbops%ntsten(grid%nfacex,grid%ngrids))

! Read the sizes of the operator stencils on each grid
pbops%nlsten = 0
pbops%nmsten = 0
pbops%njsten = 0
pbops%nhsten = 0
pbops%nrsten = 0
pbops%nrxsten = 0
pbops%nwsten = 0
pbops%ntsten = 0
READ(changrid) ((pbops%nlsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((pbops%nmsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((pbops%njsten(iv0,igrid),             &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((pbops%nhsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((pbops%nrsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((pbops%nrxsten(iv0,igrid),            &
                     iv0 = 1, grid%nvert(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((pbops%nwsten(ie0,igrid),             &
                     ie0 = 1, grid%nedge(igrid)),    &
                     igrid = 1, grid%ngrids)
READ(changrid) ((pbops%ntsten(if0,igrid),             &
                     if0 = 1, grid%nface(igrid)),    &
                     igrid = 1, grid%ngrids)

! Find maximum values in order to allocate subsequent arrays
pbops%nlsmx = MAXVAL(pbops%nlsten)
pbops%nmsmx = MAXVAL(pbops%nmsten)
pbops%njsmx = MAXVAL(pbops%njsten)
pbops%nhsmx = MAXVAL(pbops%nhsten)
pbops%nrsmx = MAXVAL(pbops%nrsten)
pbops%nrxsmx = MAXVAL(pbops%nrxsten)
pbops%nwsmx = MAXVAL(pbops%nwsten)
pbops%ntsmx = MAXVAL(pbops%ntsten)

PRINT *,'Maximum stencil sizes:'
PRINT *,'massL ...  ',pbops%nlsmx
PRINT *,'massM ...  ',pbops%nmsmx
PRINT *,'HodgeJ ... ',pbops%njsmx
PRINT *,'HodgeH ... ',pbops%nhsmx
PRINT *,'operR ...  ',pbops%nrxsmx
PRINT *,'operW ...  ',pbops%nwsmx
PRINT *,'operT ...  ',pbops%ntsmx
PRINT *,' '


! Allocate arrays for operator stencils and coefficients
allocate(pbops%lsten(grid%nfacex,pbops%nlsmx,grid%ngrids), pbops%msten(grid%nedgex,pbops%nmsmx,grid%ngrids), &
         pbops%jsten(grid%nvertx,pbops%njsmx,grid%ngrids), pbops%hsten(grid%nedgex,pbops%nhsmx,grid%ngrids), &
         pbops%rsten(grid%nfacex,pbops%nrsmx,grid%ngrids), pbops%rxsten(grid%nvertx,pbops%nrxsmx,grid%ngrids), &
         pbops%wsten(grid%nedgex,pbops%nwsmx,grid%ngrids), pbops%tsten(grid%nfacex,pbops%ntsmx,grid%ngrids))
allocate(pbops%lmass(grid%nfacex,pbops%nlsmx,grid%ngrids), pbops%mmass(grid%nedgex,pbops%nmsmx,grid%ngrids), &
         pbops%jstar(grid%nvertx,pbops%njsmx,grid%ngrids), pbops%hstar(grid%nedgex,pbops%nhsmx,grid%ngrids), &
         pbops%rcoeff(grid%nfacex,pbops%nrsmx,grid%ngrids), pbops%rxcoeff(grid%nvertx,pbops%nrxsmx,grid%ngrids), &
         pbops%wcoeff(grid%nedgex,pbops%nwsmx,grid%ngrids), pbops%tcoeff(grid%nfacex,pbops%ntsmx,pbops%ntsmx,grid%ngrids), &
         pbops%jlump(grid%nvertx,grid%ngrids), pbops%mlump(grid%nedgex,grid%ngrids), pbops%hlump(grid%nedgex,grid%ngrids))

! Read the operator stencils and coefficients
READ(changrid) (((pbops%lsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, pbops%nlsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%msten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, pbops%nmsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%jsten(iv0,ix,igrid),          &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, pbops%njsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%hsten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, pbops%nhsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%rsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, pbops%nrsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%rxsten(iv0,ix,igrid),         &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, pbops%nrxsmx),           &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%wsten(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, pbops%nwsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%tsten(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, pbops%ntsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%lmass(if0,ix,igrid),          &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, pbops%nlsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%mmass(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, pbops%nmsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%jstar(iv0,ix,igrid),          &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, pbops%njsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%hstar(ie0,ix,igrid),          &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, pbops%nhsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%rcoeff(if0,ix,igrid),         &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, pbops%nrsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%rxcoeff(iv0,ix,igrid),        &
                     iv0 = 1, grid%nvert(igrid)),    &
                     ix = 1, pbops%nrxsmx),           &
                     igrid = 1, grid%ngrids)
READ(changrid) (((pbops%wcoeff(ie0,ix,igrid),         &
                     ie0 = 1, grid%nedge(igrid)),    &
                     ix = 1, pbops%nwsmx),            &
                     igrid = 1, grid%ngrids)
READ(changrid) ((((pbops%tcoeff(if0,ix,ixx,igrid),    &
                     if0 = 1, grid%nface(igrid)),    &
                     ix = 1, pbops%ntsmx),            &
                     ixx = 1, pbops%ntsmx),           &
                     igrid = 1, grid%ngrids)

! Dimensionalize
pbops%lmass = pbops%lmass/(rearth*rearth)

! Construct the tables eoffin and eofvin
allocate(pbops%eoffin(grid%nfacex,grid%nefmx,grid%ngrids), pbops%eofvin(grid%nvertx,grid%nevmx,grid%ngrids))
do igrid = 1, grid%ngrids

  do if1 = 1, grid%nface(igrid)
    do ix = 1, grid%neoff(if1,igrid)
      ie1 = grid%eoff(if1,ix,igrid)
      if2 = grid%fnxte(ie1,1,igrid)
      IF (if1 == if2) THEN
        ! This edge points out of face if1
	      pbops%eoffin(if1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points into face if1
	      pbops%eoffin(if1,ix,igrid) = 1.0d0
      ENDIF
    enddo
  enddo

  do iv1 = 1, grid%nvert(igrid)
    do ix = 1, grid%neofv(iv1,igrid)
      ie1 = grid%eofv(iv1,ix,igrid)
      iv2 = grid%vofe(ie1,1,igrid)
      IF (iv1 == iv2) THEN
        ! This edge points away from vertex iv1
	      pbops%eofvin(iv1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points towards vertex iv1
	      pbops%eofvin(iv1,ix,igrid) = 1.0d0
      ENDIF
    enddo
  enddo

enddo


! Allocate array for size of restriction stencil
allocate(pbops%ninj(grid%nfacex,grid%ngrids-1))

! Read the size of the restriction stencil on each grid
pbops%ninj = 0
READ(changrid) ((pbops%ninj(if0,igrid),              &
                    if0 = 1, grid%nface(igrid)),    &
                    igrid = 1, grid%ngrids-1)

! Find maximum value in order to allocate subsequent arrays
pbops%ninjmx = MAXVAL(pbops%ninj)

! Allocate arrays for restriction stencils and weights
allocate(pbops%injsten(grid%nfacex,pbops%ninjmx,grid%ngrids-1))
allocate(pbops%injwgt(grid%nfacex,pbops%ninjmx,grid%ngrids-1))

! Read the restriction stencil and weights
READ(changrid) (((pbops%injsten(if0,ix,igrid),       &
                    if0 = 1, grid%nface(igrid)),    &
                    ix = 1, pbops%ninjmx),           &
                    igrid = 1, grid%ngrids-1)
READ(changrid) (((pbops%injwgt(if0,ix,igrid),        &
                    if0 = 1, grid%nface(igrid)),    &
                    ix = 1, pbops%ninjmx),           &
                    igrid = 1, grid%ngrids-1)

end subroutine readgrid

! --------------------------------------------------------------------------------------------------

subroutine centroid(grid,if0,long,lat,igrid)

! Find the centroid of cell if0 on grid igrid

implicit none
type(fempsgrid), intent(in) :: grid
integer,      intent(in)  :: if0, igrid
real*8,       intent(out) :: long, lat

integer :: ixe, ie1, iv1, iv2
real*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          xc, yc, zc, a, aby3, mag

! Coordinates of `centre' of face (i.e. dual vertex)
long1 = grid%flong(if0,igrid)
lat1 = grid%flat(if0,igrid)
call ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the face
! Hence find area of face and centroid
xc = 0.0d0
yc = 0.0d0
zc = 0.0d0
do ixe = 1, grid%neoff(if0,igrid)
  ie1 = grid%eoff(if0,ixe,igrid)
  iv1 = grid%vofe(ie1,1,igrid)
  iv2 = grid%vofe(ie1,2,igrid)
  long1 = grid%vlong(iv1,igrid)
  lat1 = grid%vlat(iv1,igrid)
  call ll2xyz(long1,lat1,x1,y1,z1)
  long1 = grid%vlong(iv2,igrid)
  lat1 = grid%vlat(iv2,igrid)
  call ll2xyz(long1,lat1,x2,y2,z2)
  call starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
  aby3 = a/3.0d0
  xc = xc + (x0 + x1 + x2)*aby3
  yc = yc + (y0 + y1 + y2)*aby3
  zc = zc + (z0 + z1 + z2)*aby3
enddo
mag = SQRT(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag
call xyz2ll(xc,yc,zc,long,lat)

end subroutine centroid

! --------------------------------------------------------------------------------------------------

subroutine dual_centroid(grid,iv0,long,lat,igrid)

! Find the centroid of dual cell iv0 on grid igrid

implicit none
type(fempsgrid), intent(in) :: grid
integer,      intent(in)  :: iv0, igrid
real*8,       intent(out) :: long, lat

integer :: ixe, ie1, iv1, iv2
real*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          xc, yc, zc, a, aby3, mag

! Coordinates of `centre' of dual cell (i.e. vertex)
long1 = grid%vlong(iv0,igrid)
lat1 = grid%vlat(iv0,igrid)
call ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the dual cell
! Hence find area of dual cell and centroid
xc = 0.0d0
yc = 0.0d0
zc = 0.0d0
do ixe = 1, grid%neofv(iv0,igrid)
  ie1 = grid%eofv(iv0,ixe,igrid)
  iv1 = grid%fnxte(ie1,1,igrid)
  iv2 = grid%fnxte(ie1,2,igrid)
  long1 = grid%flong(iv1,igrid)
  lat1 = grid%flat(iv1,igrid)
  call ll2xyz(long1,lat1,x1,y1,z1)
  long1 = grid%flong(iv2,igrid)
  lat1 = grid%flat(iv2,igrid)
  call ll2xyz(long1,lat1,x2,y2,z2)
  call starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
  aby3 = a/3.0d0
  xc = xc + (x0 + x1 + x2)*aby3
  yc = yc + (y0 + y1 + y2)*aby3
  zc = zc + (z0 + z1 + z2)*aby3
enddo
mag = SQRT(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag
call xyz2ll(xc,yc,zc,long,lat)

end subroutine dual_centroid

! --------------------------------------------------------------------------------------------------

subroutine testpoisson(grid,pbops)

! To test the solution of thePoisson problem

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(inout) :: pbops

! Number of passes
integer :: npass = 10
integer :: nf, ne, nv, if1, iv1, ipass, nprt
real*8, allocatable :: psi0(:), zeta(:), psi(:), ff1(:), ff2(:), ff3(:), ff4(:), temp1(:), temp2(:)
real*8 :: long, lat, psibar

nf = grid%nface(grid%ngrids)
ne = grid%nedge(grid%ngrids)
nv = grid%nvert(grid%ngrids)

print *,' '
print *,'--------------------------'
print *,' '
print *,'Testing mgsolve '
print *,' '

! Number of values to print for testing
nprt = 5

allocate(psi0(nf),zeta(nf),psi(nf),ff1(nf),ff2(nf),ff3(nf),ff4(nf),temp1(ne),temp2(ne))


! Build coefficients used in Laplacian operator on all grids
call buildlap(grid,pbops)


! Set up test data
! Large-scale part
do if1 = 1, nf
  long = grid%flong(if1,grid%ngrids)
  lat = grid%flat(if1,grid%ngrids)
  ! psi0(if1) = SIN(lat)
  psi0(if1) = COS(lat)*SIN(long)
enddo
! Plus small-scale part
psi0(10) = 10.0d0*psi0(10)
! Remove global mean (to ensure unambiguous result)
psibar = SUM(psi0*grid%farea(:,grid%ngrids))/SUM(grid%farea(:,grid%ngrids))
psi0 = psi0 - psibar
! Convert to area integrals
ff1 = psi0*grid%farea(:,grid%ngrids)
print *,'Original field ff1 =     ',ff1(1:nprt)
print *,' '

! Calculate laplacian
call laplace(grid,pbops,ff1,zeta,grid%ngrids,nf,ne)


! Now invert laplacian to check we get back to where we started

! Initialize result to zero
! Note psi will be stream function times grid cell area
psi = 0.0d0
temp2 = 0.0d0

! Iterate several passes
do ipass = 1, npass

  print *,'Pass ',ipass

  ! Compute residual based on latest estimate
  call massL(pbops,psi,ff2,grid%ngrids,nf)

  call Ddual1(grid,ff2,temp1,grid%ngrids,nf,ne)

  ! Improve the estimate temp2 that we obtained in the previous pass
  call massMinv(grid,pbops,temp1,temp2,grid%ngrids,ne,4)

  call Dprimal2(grid,pbops,temp2,ff3,grid%ngrids,ne,nf)

  ! Residual
  ff4 = zeta - ff3

  ! Now solve the approximate Poisson equation
  ! D2 xminv D1bar L psi' = residual
  call mgsolve(grid,pbops,ff3,ff4,grid%ngrids)

  ! And increment best estimate
  ! *** We could think about adding beta*ff3 here and tuning beta for optimal convergence ***
  psi = psi + ff3

  ! Remove global mean (to ensure unambiguous result)
  psibar = SUM(psi)/SUM(grid%farea(:,grid%ngrids))
  psi = psi - psibar*grid%farea(:,grid%ngrids)

  print *,'Original field ff1      = ',ff1(1:nprt)
  print *,'Soln of Poisson eqn psi = ', psi(1:nprt)
  ff4 = (ff1-psi)/grid%farea(:,grid%ngrids)
  print *,'RMS err in global problem = ',sqrt(sum(ff4*ff4)/nf)
  print *,' '

enddo

end subroutine testpoisson

! --------------------------------------------------------------------------------------------------

end module femps_solve_mod
