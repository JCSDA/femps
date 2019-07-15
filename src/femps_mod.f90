module femps_mod

use femps_constants_mod
use femps_utils_mod

implicit none
private

public femps

! --------------------------------------------------------------------------------------------------

type femps

! Grid file
character*31 :: ygridfile = 'gridopermap_cube_0000013824.dat'

real*8, allocatable :: lapdiag(:,:), underrel(:)

contains

 procedure :: preliminary
 procedure :: delete
 procedure :: readnml
 procedure :: allocateall
 procedure :: buildlap
 procedure :: Dprimal1
 procedure :: Dprimal2
 procedure :: Ddual1
 procedure :: Ddual2
 procedure :: buildjlump
 procedure :: buildmlump
 procedure :: buildhlump
 procedure :: buildxminv
 procedure :: massL
 procedure :: massM
 procedure :: approxMinv
 procedure :: HodgeJ
 procedure :: HodgeH
 procedure :: massLinv
 procedure :: massMinv
 procedure :: HodgeJinv
 procedure :: HodgeHinv
 procedure :: operW_original
 procedure :: operR_original
 procedure :: operW
 procedure :: operR
 procedure :: operT
 procedure :: restrict
 procedure :: prolong
 procedure :: laplace
 procedure :: xlaplace
 procedure :: residual
 procedure :: relax
 procedure :: fullmgsolve
 procedure :: mgsolve
 procedure :: readgrid
 procedure :: centroid
 procedure :: dual_centroid
 procedure :: testpoisson 

end type femps

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine preliminary(self,grid,pbobs)

! Preliminary calculations and setting up

implicit none
class(femps),     intent(inout) :: self
type(fempsgrid),  intent(inout) :: grid
type(fempspbobs), intent(inout) :: pbobs

! Read namelist information
! -------------------------
call self%readnml()
print *,'Namelists read'


! Read in the grid data
! ---------------------
call self%readgrid()
print *,'Done readgrid'


! Allocate array space now that resolution is known
! -------------------------------------------------
call self%allocateall()
print *,'Done allocateall'


! Build a lumped version of the J, M and H matrices
! -------------------------------------------------
call self%buildjlump()
call self%buildmlump()
call self%buildhlump()


! And build a local operator to approximate the inverse of M
! ----------------------------------------------------------
call self%buildxminv()

print *,'Lumped J, M and H and approximate M-inverse created'

end subroutine preliminary

! --------------------------------------------------------------------------------------------------

subroutine readnml(self)

implicit none
class(femps), intent(inout) :: self

integer, parameter :: channml = 20
character*31 :: ygridfile

! Read namelists
! --------------
NAMELIST /rundata/ ygridfile

OPEN(channml,FILE='poissonnml.in',DELIM='APOSTROPHE')
READ(channml,rundata)
CLOSE(channml)

self%ygridfile = ygridfile

end subroutine readnml

! --------------------------------------------------------------------------------------------------

subroutine allocateall(self)

implicit none
class(femps), intent(inout) :: self

! Arrays in module laplacian
! --------------------------
allocate(self%lapdiag(self%nfacex,self%ngrids))
allocate(self%underrel(self%ngrids))

end subroutine allocateall

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(femps), intent(inout) :: self

if (allocated(self%lapdiag   )) deallocate(self%lapdiag   )
if (allocated(self%underrel  )) deallocate(self%underrel  )

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine buildlap(self)

! Extract the diagonal coefficients of the Laplacian operator needed
! for the multigrid Poisson solver at all resolutions

implicit none
class(femps), intent(inout) :: self

integer :: igrid, if1, if2, if3, ie1, ie2, &
           ixd2, ixm, ixd1, ixl
real*8 :: temp, temp1(self%nfacex), cd2, cm, cd1, cl


! Under-relaxation parameter. ( *** There is scope to optimize here *** )
! -----------------------------------------------------------------------
do igrid = 1, self%ngrids
  IF (self%nefmx == 6) THEN
    self%underrel(igrid) = 0.8d0
  ELSEIF (self%nefmx == 4) THEN
    self%underrel(igrid) = 0.8d0
  ELSE
    print *,'Choose a sensible value for self%underrel in buildlap'
    STOP
  ENDIF
enddo


! Extract diagonal coefficient of Laplacian operator
do igrid = 1, self%ngrids
  ! Loop over cells
  do if1 = 1, self%nface(igrid)

    temp = 0.0d0
    ! Loop over edges of if1 involved in Dprimal2 operator
    do ixd2 = 1, self%neoff(if1,igrid)
      ie1 = self%eoff(if1,ixd2,igrid)
      cd2 = -self%eoffin(if1,ixd2,igrid)
      ! Loop over edges involved in approximate M-inverse operator
      do ixm = 1, self%nxminvsten(ie1,igrid)
        ie2 = self%xminvsten(ie1,ixm,igrid)
        cm = self%xminv(ie1,ixm,igrid)
        ! Loop over cells involved in Ddual1 operator
        cd1 = 1.0d0
        do ixd1 = 1, 2
          if2 = self%fnxte(ie2,ixd1,igrid)
          cd1 = -cd1
          ! Loop over cells in L operator
          do ixl = 1, self%nlsten(if2,igrid)
            if3 = self%lsten(if2,ixl,igrid)
            cl = self%lmass(if2,ixl,igrid)
            IF (if3 == if1) temp = temp + cd2*cm*cd1*cl
          enddo
        enddo
      enddo
    enddo
    self%lapdiag(if1,igrid) = temp

  enddo
enddo

end subroutine buildlap

! --------------------------------------------------------------------------------------------------

subroutine Dprimal1(self,f,df,igrid,nv,ne)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises pointwise values
! at vertices; df comprises integrals of the derivative
! of f along primal cell edges.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nv, ne
real*8,       intent(in)  :: f(nv)
real*8,       intent(out) :: df(ne)

integer :: ie1, iv1, iv2

do ie1 = 1, ne
  iv1 = self%vofe(ie1,1,igrid)
  iv2 = self%vofe(ie1,2,igrid)
  df(ie1) = f(iv2) - f(iv1)
enddo

end subroutine Dprimal1

! --------------------------------------------------------------------------------------------------

subroutine Dprimal2(self,f,df,igrid,ne,nf)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises integrals along
! primal edges; df comprises integrals of the derivative
! over primal cells.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne, nf
real*8,       intent(in)  :: f(ne)
real*8,       intent(out) :: df(nf)

integer :: if1, ix, ie1
real*8 :: temp

do if1 = 1, nf
  temp = 0.0d0
  do ix = 1, self%neoff(if1,igrid)
    ie1 = self%eoff(if1,ix,igrid)
    temp = temp - f(ie1)*self%eoffin(if1,ix,igrid)
  enddo
  df(if1) = temp
enddo

end subroutine Dprimal2

! --------------------------------------------------------------------------------------------------

subroutine Ddual1(self,f,df,igrid,nf,ne)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises pointwise values
! at face centres; df comprises integrals of the derivative
! of f along dual cell edges.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf)
real*8,       intent(out) :: df(ne)

integer :: ie1, if1, if2

do ie1 = 1, ne
  if1 = self%fnxte(ie1,1,igrid)
  if2 = self%fnxte(ie1,2,igrid)
  df(ie1) = f(if2) - f(if1)
enddo

end subroutine Ddual1

! --------------------------------------------------------------------------------------------------

subroutine Ddual2(self,f,df,igrid,ne,nv)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises integrals along
! dual edges; df comprises integrals of the derivative
! over dual cells.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne, nv
real*8,       intent(in)  :: f(ne)
real*8,       intent(out) :: df(nv)

integer :: iv1, ix, ie1
real*8 :: temp

do iv1 = 1, nv
  temp = 0.0d0
  do ix = 1, self%neofv(iv1,igrid)
    ie1 = self%eofv(iv1,ix,igrid)
    temp = temp + f(ie1)*self%eofvin(iv1,ix,igrid)
  enddo
  df(iv1) = temp
enddo

end subroutine Ddual2

! --------------------------------------------------------------------------------------------------

subroutine buildjlump(self)

! Build a lumped version of the j matrix

implicit none
class(femps), intent(inout) :: self

! Two choices for lumped J
! ijlump = 1     Diagonal of J
! ijlump = 2     Exact answer when applied to constant input vector
! ijlump = 3     Exact answer when applied to the discrete representation
!                  of a constant scalar (input vector proportional to self%varea)

integer, parameter :: ijlump = 3
integer :: igrid, iv1, nv
real*8, allocatable :: p1(:), jp1(:)

do igrid = 1, self%ngrids

  if (ijlump == 1) then

    do iv1 = 1, self%nvert(igrid)
      self%jlump(iv1,igrid) = self%jstar(iv1,1,igrid)
    enddo

  elseif (ijlump == 2) then

    nv = self%nvert(igrid)
    allocate(p1(nv), jp1(nv))
    p1 = 1.0d0
    call self%HodgeJ(p1,jp1,igrid,nv)
    self%jlump(1:nv,igrid) = jp1
    deallocate(p1,jp1)

  elseif (ijlump == 3) then

    nv = self%nvert(igrid)
    allocate(p1(nv), jp1(nv))
    p1 = self%varea(1:nv,igrid)
    call self%HodgeJ(p1,jp1,igrid,nv)
    self%jlump(1:nv,igrid) = jp1/self%varea(1:nv,igrid)
    deallocate(p1,jp1)

  else

    print *,'option ijlump = ',ijlump,' not available in subroutine buildjlump'
    stop

  endif

enddo


end subroutine buildjlump

! --------------------------------------------------------------------------------------------------

subroutine buildmlump(self)

! Build a lumped version of the M matrix

implicit none
class(femps), intent(inout) :: self

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

do igrid = 1, self%ngrids

  if (imlump == 1) then

    do ie1 = 1, self%nedge(igrid)
      self%mlump(ie1,igrid) = self%mmass(ie1,1,igrid)
    enddo

  elseif (imlump == 2) then

    do ie1 = 1, self%nedge(igrid)
      temp = 0.0d0
      do ix = 1, self%nmsten(ie1,igrid)
        ie2 = self%msten(ie1,ix,igrid)
        temp = temp + abs(self%mmass(ie1,ix,igrid))
      enddo
      self%mlump(ie1,igrid) = temp
    enddo

  elseif (imlump == 3) then

    nv = self%nvert(igrid)
    ne = self%nedge(igrid)
    allocate(psi1(nv), psi2(nv), psi3(nv),    &
             v1(ne), v2(ne), v3(ne), vv(ne),  &
             mv1(ne), mv2(ne), mv3(ne))
    ! Set up three solid body rotation fields
    do 	iv0 = 1, nv
      slon = sin(self%vlong(iv0,igrid))
      clon = cos(self%vlong(iv0,igrid))
      slat = sin(self%vlat(iv0,igrid))
      clat = cos(self%vlat(iv0,igrid))
      psi1(iv0) = clat*clon
      psi2(iv0) = clat*slon
      psi3(iv0) = slat
    enddo
    call self%Dprimal1(psi1,v1,igrid,nv,ne)
    call self%Dprimal1(psi2,v2,igrid,nv,ne)
    call self%Dprimal1(psi3,v3,igrid,nv,ne)
    call self%massM(v1,mv1,igrid,ne)
    call self%massM(v2,mv2,igrid,ne)
    call self%massM(v3,mv3,igrid,ne)
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
      self%mlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,vv,mv1,mv2,mv3)

  else

    print *,'option imlump = ',imlump,' not available in subroutine buildmlump'
    stop

  endif

enddo

end subroutine buildmlump

! --------------------------------------------------------------------------------------------------

subroutine buildhlump(self)

! Build a lumped version of the H matrix

implicit none
class(femps), intent(inout) :: self

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

do igrid = 1, self%ngrids

  if (ihlump == 1) then

    do ie1 = 1, self%nedge(igrid)
      self%hlump(ie1,igrid) = self%hstar(ie1,1,igrid)
    enddo

  elseif (ihlump == 2) then

    do ie1 = 1, self%nedge(igrid)
      temp = 0.0d0
      do ix = 1, self%nhsten(ie1,igrid)
        ie2 = self%hsten(ie1,ix,igrid)
        temp = temp + abs(self%hstar(ie1,ix,igrid))
      enddo
      self%hlump(ie1,igrid) = temp
    enddo

  elseif (ihlump == 3) then

    nf = self%nface(igrid)
    ne = self%nedge(igrid)
    allocate(psi1(nf), psi2(nf), psi3(nf),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three global divergent fields
    do 	if0 = 1, nf
      slon = sin(self%flong(if0,igrid))
      clon = cos(self%flong(if0,igrid))
      slat = sin(self%flat(if0,igrid))
      clat = cos(self%flat(if0,igrid))
      psi1(if0) = slat
      psi2(if0) = clat*slon
      psi3(if0) = clat*clon
    enddo
    call self%Ddual1(psi1,v1,igrid,nf,ne)
    call self%Ddual1(psi2,v2,igrid,nf,ne)
    call self%Ddual1(psi3,v3,igrid,nf,ne)
    call self%HodgeH(v1,hv1,igrid,ne)
    call self%HodgeH(v2,hv2,igrid,ne)
    call self%HodgeH(v3,hv3,igrid,ne)
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
      self%hlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  elseif (ihlump == 4) then

    nv = self%nvert(igrid)
    ne = self%nedge(igrid)
    allocate(psi1(nv), psi2(nv), psi3(nv),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three solid body rotation fields
    do 	iv0 = 1, nv
      slon = sin(self%vlong(iv0,igrid))
      clon = cos(self%vlong(iv0,igrid))
      slat = sin(self%vlat(iv0,igrid))
      clat = cos(self%vlat(iv0,igrid))
      psi1(iv0) = slat
      psi2(iv0) = clat*slon
      psi3(iv0) = clat*clon
    enddo
    call self%Dprimal1(psi1,v1,igrid,nv,ne)
    call self%Dprimal1(psi2,v2,igrid,nv,ne)
    call self%Dprimal1(psi3,v3,igrid,nv,ne)
    call self%HodgeH(v1,hv1,igrid,ne)
    call self%HodgeH(v2,hv2,igrid,ne)
    call self%HodgeH(v3,hv3,igrid,ne)
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
      self%hlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  else

    print *,'option ihlump = ',ihlump,' not available in subroutine buildhlump'
    stop

  endif

enddo

end subroutine buildhlump

! --------------------------------------------------------------------------------------------------

subroutine buildxminv(self)

! Determine stencil and coefficients for a local approximation
! to the inverse of M.
!
! Two options are coded. The first is simply the inverse of
! the diagonal mass lumped approximation to M. The second is
! based on a single underrelaxed Jacobi iteration to the
! inverse of M.

implicit none
class(femps), intent(inout) :: self

logical :: llump
integer :: igrid, ie0, ix, ie1
real*8 :: temp, diag
real*8 :: relax

! Underrelaxation coefficient depends on grid
! 0.9 and 1.4 when using self%mlump in Jacobi
! Check for consistency with massMinv
if (self%nefmx == 4) then
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
  self%nxmisx = 1
  allocate(self%nxminvsten(self%nedgex,self%ngrids), self%xminvsten(self%nedgex,self%nxmisx,self%ngrids), &
           self%xminv(self%nedgex,self%nxmisx,self%ngrids))
  do igrid = 1, self%ngrids
    do ie0 = 1, self%nedge(igrid)
      ! stencil for edge ie0 is ie0 itself, and coeff is
      ! inverse of the diagonal term of the lumped matrix
      self%nxminvsten(ie0,igrid) = 1
      self%xminvsten(ie0,1,igrid) = ie0
      self%xminv(ie0,1,igrid) = 1.0d0/self%mlump(ie0,igrid)
    enddo
  enddo

else

  ! Stencil is the same as for M itself
  self%nxmisx = self%nmsmx
  allocate(self%nxminvsten(self%nedgex,self%ngrids), self%xminvsten(self%nedgex,self%nxmisx,self%ngrids), &
           self%xminv(self%nedgex,self%nxmisx,self%ngrids))
  do igrid = 1, self%ngrids
    do ie0 = 1, self%nedge(igrid)
      ! stencil for edge ie0 is the same as the stencil for m
      self%nxminvsten(ie0,igrid) = self%nmsten(ie0,igrid)
      do ix = 1, self%nmsten(ie0,igrid)
        ie1 = self%msten(ie0,ix,igrid)
        self%xminvsten(ie0,ix,igrid) = ie1
        if (ie1 == ie0) then
          diag = 1.0d0 + relax
        else
          diag = 0.0d0
        endif
        temp = self%mmass(ie0,ix,igrid)/self%mlump(ie1,igrid)
        self%xminv(ie0,ix,igrid) = (diag - relax*temp)/self%mlump(ie0,igrid)
      enddo
    enddo
  enddo

endif

end subroutine buildxminv

! --------------------------------------------------------------------------------------------------

subroutine massL(self,f1,f2,igrid,nf)

! Apply the mass matrix L to field f1 to obtain the result f2

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nf
real*8,       intent(in)  :: f1(nf)
real*8,       intent(out) :: f2(nf)
integer :: if1, if2, ix
real*8 :: temp

do if1 = 1, nf
  temp = 0.0d0
  do ix = 1, self%nlsten(if1,igrid)
    if2 = self%lsten(if1,ix,igrid)
    temp = temp + f1(if2)*self%lmass(if1,ix,igrid)
  enddo
  f2(if1) = temp
enddo

end subroutine massL

! --------------------------------------------------------------------------------------------------

subroutine massM(self,f1,f2,igrid,ne)

! Apply the mass matrix M to field f1 to obtain field f2

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)
integer :: ie1, ie2, ix
real*8 :: temp

do ie1 = 1, ne
  temp = 0.0d0
  do ix = 1, self%nmsten(ie1,igrid)
    ie2 = self%msten(ie1,ix,igrid)
    temp = temp + f1(ie2)*self%mmass(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine massM

! --------------------------------------------------------------------------------------------------

subroutine approxMinv(self,f1,f2,igrid,ne)

! Apply an approximate inverse of the mass matrix M
! to field f1 to obtain field f2

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: ie1, ie2, ix
real*8 :: temp

do ie1 = 1, ne
  temp = 0.0d0
  do ix = 1, self%nxminvsten(ie1,igrid)
    ie2 = self%xminvsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*self%xminv(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine approxMinv

! --------------------------------------------------------------------------------------------------

subroutine HodgeJ(self,f1,f2,igrid,nv)

! Apply the Hodge star J operator that converts dual face
! integrals f1 to vertex values f2 on grid igrid.

implicit none
class(femps), intent(in) :: self
integer,      intent(in) :: igrid, nv
real*8,       intent(in) :: f1(nv)
real*8,       intent(out) :: f2(nv)

integer :: iv1, iv2, ix
real*8 :: temp

do iv1 = 1, nv
  temp = 0.0d0
  do ix = 1, self%njsten(iv1,igrid)
    iv2 = self%jsten(iv1,ix,igrid)
    temp = temp + f1(iv2)*self%jstar(iv1,ix,igrid)
  enddo
  f2(iv1) = temp
enddo

end subroutine HodgeJ

! --------------------------------------------------------------------------------------------------

subroutine HodgeH(self,f1,f2,igrid,ne)

! Apply the Hodge star H operator that converts dual edge
! integrals (circulations) f1 to primal edge integrals (fluxes) f2
! on grid igrid.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: ie1, ie2, ix
real*8 :: temp

do ie1 = 1, ne
  temp = 0.0d0
  do ix = 1, self%nhsten(ie1,igrid)
    ie2 = self%hsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*self%hstar(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine HodgeH

! --------------------------------------------------------------------------------------------------

subroutine massLinv(self,f1,f2,igrid,nf,niter)

! Apply the inverse of the mass matrix L to field f1 to obtain
! the result f2.

! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If L is diagonal then there is no need to iterate.

implicit none
class(femps), intent(in)    :: self
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, nf
real*8,       intent(in)    :: f1(nf)
real*8,       intent(inout) :: f2(nf)

integer :: if1, iter, miter
real*8 :: temp(nf)

miter = ABS(niter)

IF (niter < 0 .OR. self%nlsmx == 1) THEN
  ! First guess based on diagonal L
  do if1 = 1, nf
    f2(if1) = f1(if1)/self%lmass(if1,1,igrid)
  enddo
ENDIF

IF (self%nlsmx > 1) THEN
  ! L is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call self%massL(f2,temp,igrid,nf)
    temp = f1 - temp
    do if1 = 1, nf
      f2(if1) = f2(if1) + temp(if1)/self%lmass(if1,1,igrid)
    enddo
  enddo
ENDIF

end subroutine massLinv

! --------------------------------------------------------------------------------------------------

subroutine massMinv(self,f1,f2,igrid,ne,niter)

! Apply the inverse of the mass matrix M to the field f1
! to obtain the result f2

! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If M is diagonal then there is no need to iterate.

implicit none
class(femps), intent(in)    :: self
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, ne
real*8,       intent(in)    :: f1(ne)
real*8,       intent(inout) :: f2(ne)

integer :: ie1, iter, miter
real*8 :: temp(ne)
real*8 :: relax

! Underrelaxation coefficient depends on grid
IF (self%nefmx == 4) THEN
  relax = 0.9d0
ELSE
  relax = 1.4d0
ENDIF

miter = ABS(niter)

IF (niter < 0 .OR. self%nmsmx == 1) THEN
  ! First guess based on lumped M
  do ie1 = 1, ne
    f2(ie1) = f1(ie1)/self%mlump(ie1,igrid)
  enddo
ENDIF

IF (self%nmsmx > 1) THEN
  ! M is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call self%massM(f2,temp,igrid,ne)
    temp = relax*(f1 - temp)
    do ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/self%mlump(ie1,igrid)
    enddo
  enddo
ENDIF

end subroutine massMinv

! --------------------------------------------------------------------------------------------------

subroutine HodgeJinv(self,f1,f2,igrid,nv,niter)

! Apply the inverse Hodge star operator J^{-1} that maps from
! E_p to V_d on grid igrid.
!
! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If J is diagonal then there is no need to iterate.

implicit none
class(femps), intent(in)    :: self
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, nv
real*8,       intent(in)    :: f1(nv)
real*8,       intent(inout) :: f2(nv)

integer :: iv1, iter, miter
real*8 :: temp(nv)
real*8 :: relax = 1.4d0 ! relax = 1.4 is good for ijlump = 3 on hex and cube grids 

miter = ABS(niter)

IF (niter < 0 .OR. self%njsmx == 1) THEN
  ! First guess based on lumped J
  do iv1 = 1, nv
    f2(iv1) = f1(iv1)/self%jlump(iv1,igrid)
  enddo
ENDIF

IF (self%njsmx > 1) THEN
  ! J is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call self%HodgeJ(f2,temp,igrid,nv)
    temp = relax*(f1 - temp)
    do iv1 = 1, nv
      f2(iv1) = f2(iv1) + temp(iv1)/self%jlump(iv1,igrid)
    enddo
  enddo
ENDIF

end subroutine HodgeJinv

! --------------------------------------------------------------------------------------------------

subroutine HodgeHinv(self,f1,f2,igrid,ne,niter)

! Apply the inverse Hodge star operator H^{-1} that maps from
! S_p to S_d on grid igrid.
!
! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If H is diagonal then there is no need to iterate.

implicit none
class(femps), intent(in)    :: self
integer,      intent(in)    :: niter
integer,      intent(in)    :: igrid, ne
real*8,       intent(in)    :: f1(ne)
real*8,       intent(inout) :: f2(ne)

integer :: ie1, iter, miter
real*8 :: temp(ne)
real*8 :: relax = 1.4d0 ! relax = 1.4 is good for ihlump = 3 on hex and cube grids 

miter = ABS(niter)

IF (niter < 0 .OR. self%nhsmx == 1) THEN
  ! First guess based on diagonal H
  do ie1 = 1, ne
    f2(ie1) = f1(ie1)/self%hlump(ie1,igrid)
  enddo
ENDIF

IF (self%nhsmx > 1) THEN
  ! H is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call self%HodgeH(f2,temp,igrid,ne)
    temp = relax*(f1 - temp)
    do ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/self%hlump(ie1,igrid)
    enddo
  enddo
ENDIF

end subroutine HodgeHinv

! --------------------------------------------------------------------------------------------------

subroutine operW_original(self,f1,f2,igrid,ne)

! Apply the W operator:
! given fluxes f1 across primal edges, construct
! the rotated fluxes across dual edges f2, on grid igrid.

! This is the original formulation, building W from R
! a la TRiSK. It probably requires an MPI reduce operation
! so is likely to be inefficient.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: if1, ie1, ie2, ix, ix1, ix2, ixv, ne1
real*8 :: ss, w

! Initialize to zero
f2 = 0.0d0

! Loop over faces
do if1 = 1, self%nface(igrid)
  ne1 = self%neoff(if1,igrid)
  ! For each edge of this face
  do ix1 = 1, ne1
    ss = -0.5
    ie1 = self%eoff(if1,ix1,igrid)
    ! Find the contribution to f2 from every other
    ! edge of this face
    do ix = 0, ne1 - 2
      ixv = MOD(ix1 + ix - 1,ne1) + 1
      ix2 = MOD(ix1 + ix,ne1) + 1
      ie2 = self%eoff(if1,ix2,igrid)
      ss = ss + self%rcoeff(if1,ixv,igrid)
      w = -ss*self%eoffin(if1,ix1,igrid)*self%eoffin(if1,ix2,igrid)
      f2(ie1) = f2(ie1) + w*f1(ie2)
    enddo
  enddo
enddo

end subroutine operW_original

! --------------------------------------------------------------------------------------------------

subroutine operR_original(self,f1,f2,igrid,nf,nv)

! Apply the R operator:
! map from V_p to E_p

! This is the original formulation. It loops over `source'
! entities rather than target entities and so will require
! an MPI reduce.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nf, nv
real*8,       intent(in)  :: f1(nf)
real*8,       intent(out) :: f2(nv)

integer :: if1, iv1, ix1, ne1

! Initialize to zero
f2 = 0.0d0

! Loop over faces
do if1 = 1, self%nface(igrid)
  ne1 = self%neoff(if1,igrid)
  ! Share out this face's contributions to its surrounding vertices
  do ix1 = 1, ne1
    iv1 = self%rsten(if1,ix1,igrid)
    f2(iv1) = f2(iv1) + f1(if1)*self%rcoeff(if1,ix1,igrid)
  enddo
enddo

end subroutine operR_original

! --------------------------------------------------------------------------------------------------

subroutine operW(self,f1,f2,igrid,ne)

! Apply the W operator:
! given edge integrals of normal components f1 on primal edges,
! construct edge integrals of normal components of perpendicular
! field f2, on grid igrid.

! This formulation uses pre-build stencil and coefficients to
! avoid the need for MPI reduce. It is mathematically equivalent
! to the original formulation.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(ne)

integer :: ie0, ne1, ix1, ie1
real*8 :: temp

! Loop over vertices
do ie0 = 1, self%nedge(igrid)
  ne1 = self%nwsten(ie0,igrid)
  ! Collect contributions from stencil
  temp = 0.0d0
  do ix1 = 1, ne1
    ie1 = self%wsten(ie0,ix1,igrid)
    temp = temp + f1(ie1)*self%wcoeff(ie0,ix1,igrid)
  enddo
  f2(ie0) = temp
enddo

end subroutine operW

! --------------------------------------------------------------------------------------------------

subroutine operR(self,f1,f2,igrid,nf,nv)

! Apply the R operator:
! given face integrals f1 on primal faces, map to dual cell
! integrals f2, on grid igrid.

! This formulation stores the coefficients in the transpose of
! the original formulation to avoid an MPI reduce. It is
! mathematically equivalent to the original formulation.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nf, nv
real*8,       intent(in)  :: f1(nf)
real*8,       intent(out) :: f2(nv)

integer :: iv0, if1, ix1, ne1
real*8 :: temp

! Loop over vertices
do iv0 = 1, self%nvert(igrid)
  ne1 = self%nrxsten(iv0,igrid)
  ! Collect contributions from surrounding faces
  temp = 0.0d0
  do ix1 = 1, ne1
    if1 = self%rxsten(iv0,ix1,igrid)
    temp = temp + f1(if1)*self%rxcoeff(iv0,ix1,igrid)
  enddo
  f2(iv0) = temp
enddo

end subroutine operR

! --------------------------------------------------------------------------------------------------

subroutine operT(self,f1,f2,igrid,ne,nf)

! Apply the T operator:
! compute cell integrals of 2 x kinetic energy from edge integrals
! of normal fluxes

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, ne, nf
real*8,       intent(in)  :: f1(ne)
real*8,       intent(out) :: f2(nf)

integer :: if1, ix1, ix2, ne1, ie1, ie2
real*8 :: temp

! Loop over faces
do if1 = 1, self%nface(igrid)
  ne1 = self%ntsten(if1,igrid)
  temp = 0.0d0
  ! Loop over all pairs of edges of this cell
  do ix1 = 1, ne1
    ie1 = self%tsten(if1,ix1,igrid)
    do ix2 = 1, ne1
      ie2 = self%tsten(if1,ix2,igrid)
      temp = temp + self%tcoeff(if1,ix1,ix2,igrid)*f1(ie1)*f1(ie2)
    enddo
  enddo
  f2(if1) = temp
enddo

end subroutine operT

! --------------------------------------------------------------------------------------------------

subroutine restrict(self,f1,nf1,f2,nf2,igrid)

! To perform the restriction operation needed for a multigrid solver.
! Restrict field f1 from grid igrid + 1 to grid igrid and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms).

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: nf1, nf2, igrid
real*8,       intent(in)  :: f1(nf1)
real*8,       intent(out) :: f2(nf2)

integer :: if1, if2, ix
real*8 :: wgt

! Safety check
IF (nf2 .ne. self%nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine restrict'
  STOP
ENDIF

do if2 = 1, nf2
  f2(if2) = 0.0d0
  do ix = 1, self%ninj(if2,igrid)
    if1 = self%injsten(if2,ix,igrid)
    wgt = self%injwgt(if2,ix,igrid)
    f2(if2) = f2(if2) + wgt*f1(if1)
  enddo
enddo

end subroutine restrict

! --------------------------------------------------------------------------------------------------

subroutine prolong(self,f2,nf2,f1,nf1,igrid)

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
class(femps), intent(in)  :: self
integer,      intent(in)  :: nf1, nf2, igrid
real*8,       intent(in)  :: f2(nf2)
real*8,       intent(out) :: f1(nf1)
integer :: if1, if2, ix, igridp
real*8 :: wgt, f2if2, temp1(nf1), temp2(nf2)

! Safety check
IF (nf2 .ne. self%nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine prolong'
  STOP
ENDIF

igridp = igrid + 1
temp2(1:nf2) = f2(1:nf2)/self%farea(1:nf2,igrid)
temp1 = 0.0d0
do if2 = 1, nf2
  f2if2 = temp2(if2)
  do ix = 1, self%ninj(if2,igrid)
    if1 = self%injsten(if2,ix,igrid)
    wgt = self%injwgt(if2,ix,igrid)
    temp1(if1) = temp1(if1) + wgt*f2if2
  enddo
enddo
f1(1:nf1) = temp1(1:nf1)*self%farea(1:nf1,igridp)

end subroutine prolong

! --------------------------------------------------------------------------------------------------

subroutine laplace(self,f,hf,igrid,nf,ne)

! To apply the Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf)
real*8,       intent(out) :: hf(nf)

real*8 :: temp1(nf), temp2(ne), temp3(ne)
integer :: niter

call self%massL(f,temp1,igrid,nf)
call self%Ddual1(temp1,temp2,igrid,nf,ne)
niter = -20
call self%massMinv(temp2,temp3,igrid,ne,niter)
call self%Dprimal2(temp3,hf,igrid,ne,nf)


end subroutine laplace

! --------------------------------------------------------------------------------------------------

subroutine xlaplace(self,f,hf,igrid,nf,ne)

! To apply the APPROXIMATE Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf)
real*8,       intent(out) :: hf(nf)

real*8 :: temp1(nf), temp2(ne), temp3(ne)

call self%massL(f,temp1,igrid,nf)
call self%Ddual1(temp1,temp2,igrid,nf,ne)
call self%approxMinv(temp2,temp3,igrid,ne)
call self%Dprimal2(temp3,hf,igrid,ne,nf)

end subroutine xlaplace

! --------------------------------------------------------------------------------------------------

subroutine residual(self,f,rhs,res,igrid,nf,ne)

! Compute the residual res in the approximate Poisson equation on grid igrid
! when f is the input field and rhs is the right hand side. Note that
! f, rhs and res are area integrals (2-forms).

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: igrid, nf, ne
real*8,       intent(in)  :: f(nf), rhs(nf)
real*8,       intent(out) :: res(nf)

call self%xlaplace(f,res,igrid,nf,ne)
res = rhs - res
!print *,'     residual: ',res(1:5)

end subroutine residual

! --------------------------------------------------------------------------------------------------

subroutine relax(self,f,rhs,igrid,nf,ne,niter)

! To carry out niter Jacobi relaxation iterations for the multigrid
! solver on grid igrid

implicit none
class(femps), intent(in)    :: self
integer,      intent(in)    :: igrid, nf, ne, niter
real*8,       intent(in)    :: rhs(nf)
real*8,       intent(inout) :: f(nf)

real*8, allocatable :: res(:), finc(:)
real*8 :: u
integer :: iter

allocate(res(nf), finc(nf))

u = self%underrel(igrid)
do iter = 1, niter
  call self%residual(f,rhs,res,igrid,nf,ne)
  finc = res/self%lapdiag(1:nf,igrid)
  f = f + u*finc
enddo

deallocate(res, finc)

end subroutine relax

! --------------------------------------------------------------------------------------------------

subroutine fullmgsolve(self,phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 self%xMinv Ddual1 L phi = RHS
!
! using full multigrid algorithm.
! Coefficients are contained in module laplace.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: ng
real*8,       intent(in)  :: rr(self%nfacex)
real*8,       intent(out) :: phi(self%nfacex)

! Numbers of iterations on coarsest grid and other grids
integer, parameter :: niterc = 10, niter = 2, npass = 1
integer :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
real*8, allocatable :: ff(:,:), rf(:,:)
real*8 :: temp1(self%nfacex)

! Allocate space on all grids
allocate(ff(self%nfacex,self%ngrids),rf(self%nfacex,self%ngrids))

! One pass should be enough. Warn user if npass is set to
! some other value for testing
IF (npass > 1) PRINT *,'mgsolve: npass = ',npass

! Initialize solution to zero
phi = 0.0d0

! For diagnostics
!nf1 = self%nface(self%ngrids)
!ne1 = self%nedge(self%ngrids)
!call self%residual(phi,rr,temp1,self%ngrids,nf1,ne1)
!print *,'Pass ',0,'  RMS residual = ',SQRT(SUM(temp1*temp1)/nf1)

do ipass = 1, npass

  ! Initialize rhs as residual using latest estimate
  IF (ipass == 1) THEN
    ! No need to do the calculation
    rf(:,self%ngrids) = rr(:)
  ELSE
    nf1 = self%nface(self%ngrids)
    ne1 = self%nedge(self%ngrids)
    call self%residual(phi,rr,rf(1,self%ngrids),self%ngrids,nf1,ne1)
  ENDIF

  ! Initialize correction to solution on all grids to zero
  ff = 0.0d0

  ! Restrict right hand side to each grid in the hierarchy
  do igrid = self%ngrids-1, self%ngrids-ng+1, -1
    igridp = igrid + 1
    nf1 = self%nface(igridp)
    nf2 = self%nface(igrid)
    call self%restrict(rf(1,igridp),nf1,rf(1,igrid),nf2,igrid)
  enddo

  ! Iterate to convergence on coarsest grid
  igrid = self%ngrids-ng+1
  nf1 = self%nface(igrid)
  ne1 = self%nedge(igrid)
  ff(1:nf1,igrid) = 0.0d0
  call self%relax(ff(1,igrid),rf(1,igrid),igrid,nf1,ne1,niterc)

  ! Sequence of growing V-cycles
  do igridp = self%ngrids-ng+2, self%ngrids

    igrid = igridp - 1
    nf1 = self%nface(igridp)
    ne1 = self%nedge(igridp)
    nf2 = self%nface(igrid)
    ne2 = self%nedge(igrid)

    ! Prolong solution to grid igridp
    ! and execute one V-cycle starting from grid igridp

    ! Prolong
    call self%prolong(ff(1,igrid),nf2,ff(1,igridp),nf1,igrid)

    ! Descending part of V-cycle
    do jgrid = igrid, self%ngrids-ng+1, -1
    
      jgridp = jgrid + 1
      nf1 = self%nface(jgridp)
      ne1 = self%nedge(jgridp)
      nf2 = self%nface(jgrid)
      ne2 = self%nedge(jgrid)

      ! Relax on grid jgridp
      call self%relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

      ! Calculate residual on jgridp
      call self%residual(ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)
  
      ! Restrict residual to jgrid
      call self%restrict(temp1,nf1,rf(1,jgrid),nf2,jgrid)

      ! Set correction first guess to zero on grid jgrid-1
      ff(1:nf2,jgrid) = 0.0d0

    enddo

    ! Iterate to convergence on coarsest grid
    jgrid = self%ngrids-ng+1
    nf1 = self%nface(jgrid)
    ne1 = self%nedge(jgrid)
    ff(1:nf1,jgrid) = 0.0d0
    call self%relax(ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

    ! Ascending part of V-cycle
    do jgrid = self%ngrids-ng+1, igrid

      jgridp = jgrid + 1
      igrid = igrid - 1
      nf1 = self%nface(jgridp)
      ne1 = self%nedge(jgridp)
      nf2 = self%nface(jgrid)
      ne2 = self%nedge(jgrid)

      ! Prolong correction to jgridp
      call self%prolong(ff(1,jgrid),nf2,temp1,nf1,jgrid)

      ! Add correction to solution on jgridp
      ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

      ! Relax on grid jgridp
      call self%relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

    enddo

  enddo

  ! Add correction to phi
  phi = phi + ff(:,self%ngrids)

  ! For diagnostics
  !nf1 = self%nface(self%ngrids)
  !ne1 = self%nedge(self%ngrids)
  !call self%residual(phi,rr,temp1,self%ngrids,nf1,ne1)
  !print *,'      RMS residual in fullmgsolve = ',SQRT(SUM(temp1*temp1)/nf1)

enddo


deallocate(ff,rf)

end subroutine fullmgsolve

! --------------------------------------------------------------------------------------------------

subroutine mgsolve(self,phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 self%xMinv Ddual1 L phi = RHS
!
! using a single V-cycle multigrid algorithm.
! Coefficients are contained in module laplace.

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: ng
real*8,       intent(in)  :: rr(self%nfacex)
real*8,       intent(out) :: phi(self%nfacex)

! Numbers of iterations on coarsest grid and other grids
integer, parameter :: niterc = 10, niter = 2
integer :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
real*8, allocatable :: ff(:,:), rf(:,:)
real*8 :: temp1(self%nfacex)

! Allocate space on all grids
allocate(ff(self%nfacex,self%ngrids),rf(self%nfacex,self%ngrids))

! Initialize solution to zero
phi = 0.0d0

! Initialize rhs on finest grid
rf(:,self%ngrids) = rr(:)

! Initialize correction to solution on all grids to zero
ff = 0.0d0


! Descending part of V-cycle
do jgrid = self%ngrids-1, self%ngrids-ng+1, -1
    
  jgridp = jgrid + 1
  nf1 = self%nface(jgridp)
  ne1 = self%nedge(jgridp)
  nf2 = self%nface(jgrid)
  ne2 = self%nedge(jgrid)

  ! Relax on grid jgridp
  call self%relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

  ! Calculate residual on jgridp
  call self%residual(ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)
  
  ! Restrict residual to jgrid
  call self%restrict(temp1,nf1,rf(1,jgrid),nf2,jgrid)

  ! Set correction first guess to zero on grid jgrid
  ff(1:nf2,jgrid) = 0.0d0

enddo

! Iterate to convergence on coarsest grid
jgrid = self%ngrids-ng+1
nf1 = self%nface(jgrid)
ne1 = self%nedge(jgrid)
ff(1:nf1,jgrid) = 0.0d0
call self%relax(ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

! Ascending part of V-cycle
do jgrid = self%ngrids-ng+1, self%ngrids-1

  jgridp = jgrid + 1
  nf1 = self%nface(jgridp)
  ne1 = self%nedge(jgridp)
  nf2 = self%nface(jgrid)
  ne2 = self%nedge(jgrid)

  ! Prolong correction to jgridp
  call self%prolong(ff(1,jgrid),nf2,temp1,nf1,jgrid)

  ! Add correction to solution on jgridp
  ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

  ! Relax on grid jgridp
  call self%relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

enddo

! Add correction to phi
phi = phi + ff(:,self%ngrids)

!  ! For diagnostics
!  nf1 = self%nface(self%ngrids)
!  ne1 = self%nedge(self%ngrids)
!  call residual(phi,rr,temp1,self%ngrids,nf1,ne1)
!  print *,'     RMS residual in mgsolve = ',SQRT(SUM(temp1*temp1)/nf1)

deallocate(ff,rf)

end subroutine mgsolve

! --------------------------------------------------------------------------------------------------

subroutine readgrid(self)

! To allocate array space for the grid information in module grid
! and to read the information from file

implicit none
class(femps), intent(inout)  :: self

integer :: if0, ie0, iv0, igrid, ix, ixx, if1, if2, iv1, iv2, ie1
integer, parameter :: changrid = 25         ! grid information

! Open file for reading
OPEN(changrid,FILE=self%ygridfile,FORM='UNFORMATTED')

! First read self%ngrids
READ(changrid) self%ngrids


! Allocate self%nface, self%nedge, self%nvert
allocate(self%nface(self%ngrids), self%nedge(self%ngrids), self%nvert(self%ngrids))

! Read numbers of faces, edges and vertices on each grid
READ(changrid) self%nface
READ(changrid) self%nedge
READ(changrid) self%nvert

! Find maximum values in order to allocated subsequent arrays
self%nfacex = MAXVAL(self%nface)
self%nedgex = MAXVAL(self%nedge)
self%nvertx = MAXVAL(self%nvert)

! Allocate self%neoff, self%neofv
allocate(self%neoff(self%nfacex,self%ngrids), self%neofv(self%nvertx,self%ngrids))

! Read the numbers of edges of each face and vertex on each grid
self%neoff = 0
self%neofv = 0
READ(changrid) ((self%neoff(if0,igrid),          &
                    if0 = 1, self%nface(igrid)), &
                    igrid = 1, self%ngrids)
READ(changrid) ((self%neofv(iv0,igrid),          &
                    iv0 = 1, self%nvert(igrid)), &
                    igrid = 1, self%ngrids)

! Find maximum values in order to allocate subsequent arrays
self%nefmx = MAXVAL(self%neoff)
self%nevmx = MAXVAL(self%neofv)


! Allocate connectivity arrays arrays
allocate(self%fnxtf(self%nfacex,self%nefmx,self%ngrids), self%eoff(self%nfacex,self%nefmx,self%ngrids), &
         self%voff(self%nfacex,self%nefmx,self%ngrids),  self%fnxte(self%nedgex,2,self%ngrids),    &
         self%vofe(self%nedgex,2,self%ngrids),      self%fofv(self%nvertx,self%nevmx,self%ngrids), &
         self%eofv(self%nvertx,self%nevmx,self%ngrids))

! Read the connectivity arrays
READ(changrid) (((self%fnxtf(if0,ix,igrid),          &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%nefmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%eoff(if0,ix,igrid),           &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%nefmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%voff(if0,ix,igrid),           &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%nefmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%fnxte(ie0,ix,igrid),          &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%vofe(ie0,ix,igrid),           &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%fofv(iv0,ix,igrid),           &
                     iv0 = 1, self%nvert(igrid)),    &
                     ix = 1, self%nevmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%eofv(iv0,ix,igrid),           &
                     iv0 = 1, self%nvert(igrid)),    &
                     ix = 1, self%nevmx),            &
                     igrid = 1, self%ngrids)


! Allocate the geometrical information arrays
allocate(self%flong(self%nfacex,self%ngrids), self%flat(self%nfacex,self%ngrids),  &
         self%vlong(self%nvertx,self%ngrids), self%vlat(self%nvertx,self%ngrids),  &
         self%farea(self%nfacex,self%ngrids), self%varea(self%nvertx,self%ngrids), &
         self%ldist(self%nedgex,self%ngrids), self%ddist(self%nedgex,self%ngrids), &
         self%fareamin(self%ngrids))

! Read the geometrical information arrays
READ(changrid) ((self%flong(if0,igrid),              &
                     if0 = 1, self%nface(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%flat(if0,igrid),               &
                     if0 = 1, self%nface(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%vlong(iv0,igrid),              &
                     iv0 = 1, self%nvert(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%vlat(iv0,igrid),               &
                     iv0 = 1, self%nvert(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%farea(if0,igrid),              &
                     if0 = 1, self%nface(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%varea(iv0,igrid),              &
                     iv0 = 1, self%nvert(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%ldist(ie0,igrid),              &
                     ie0 = 1, self%nedge(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%ddist(ie0,igrid),              &
                     ie0 = 1, self%nedge(igrid)),    &
                     igrid = 1, self%ngrids)

! Dimensionalize
self%farea = self%farea*rearth*rearth
self%varea = self%varea*rearth*rearth
self%ldist = self%ldist*rearth
self%ddist = self%ddist*rearth

! Determine smallest face area on each grid
do igrid = 1, self%ngrids
  self%fareamin(igrid) = MINVAL(self%farea(1:self%nface(igrid),igrid))
enddo


! Allocate arrays for size of operator stencils
allocate(self%nlsten(self%nfacex,self%ngrids), self%nmsten(self%nedgex,self%ngrids), &
         self%njsten(self%nvertx,self%ngrids), self%nhsten(self%nedgex,self%ngrids), &
         self%nrsten(self%nfacex,self%ngrids), self%nrxsten(self%nvertx,self%ngrids), &
         self%nwsten(self%nedgex,self%ngrids), self%ntsten(self%nfacex,self%ngrids))

! Read the sizes of the operator stencils on each grid
self%nlsten = 0
self%nmsten = 0
self%njsten = 0
self%nhsten = 0
self%nrsten = 0
self%nrxsten = 0
self%nwsten = 0
self%ntsten = 0
READ(changrid) ((self%nlsten(if0,igrid),             &
                     if0 = 1, self%nface(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%nmsten(ie0,igrid),             &
                     ie0 = 1, self%nedge(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%njsten(iv0,igrid),             &
                     iv0 = 1, self%nvert(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%nhsten(ie0,igrid),             &
                     ie0 = 1, self%nedge(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%nrsten(if0,igrid),             &
                     if0 = 1, self%nface(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%nrxsten(iv0,igrid),            &
                     iv0 = 1, self%nvert(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%nwsten(ie0,igrid),             &
                     ie0 = 1, self%nedge(igrid)),    &
                     igrid = 1, self%ngrids)
READ(changrid) ((self%ntsten(if0,igrid),             &
                     if0 = 1, self%nface(igrid)),    &
                     igrid = 1, self%ngrids)

! Find maximum values in order to allocate subsequent arrays
self%nlsmx = MAXVAL(self%nlsten)
self%nmsmx = MAXVAL(self%nmsten)
self%njsmx = MAXVAL(self%njsten)
self%nhsmx = MAXVAL(self%nhsten)
self%nrsmx = MAXVAL(self%nrsten)
self%nrxsmx = MAXVAL(self%nrxsten)
self%nwsmx = MAXVAL(self%nwsten)
self%ntsmx = MAXVAL(self%ntsten)

PRINT *,'Maximum stencil sizes:'
PRINT *,'massL ...  ',self%nlsmx
PRINT *,'massM ...  ',self%nmsmx
PRINT *,'HodgeJ ... ',self%njsmx
PRINT *,'HodgeH ... ',self%nhsmx
PRINT *,'operR ...  ',self%nrxsmx
PRINT *,'operW ...  ',self%nwsmx
PRINT *,'operT ...  ',self%ntsmx
PRINT *,' '


! Allocate arrays for operator stencils and coefficients
allocate(self%lsten(self%nfacex,self%nlsmx,self%ngrids), self%msten(self%nedgex,self%nmsmx,self%ngrids), &
         self%jsten(self%nvertx,self%njsmx,self%ngrids), self%hsten(self%nedgex,self%nhsmx,self%ngrids), &
         self%rsten(self%nfacex,self%nrsmx,self%ngrids), self%rxsten(self%nvertx,self%nrxsmx,self%ngrids), &
         self%wsten(self%nedgex,self%nwsmx,self%ngrids), self%tsten(self%nfacex,self%ntsmx,self%ngrids))
allocate(self%lmass(self%nfacex,self%nlsmx,self%ngrids), self%mmass(self%nedgex,self%nmsmx,self%ngrids), &
         self%jstar(self%nvertx,self%njsmx,self%ngrids), self%hstar(self%nedgex,self%nhsmx,self%ngrids), &
         self%rcoeff(self%nfacex,self%nrsmx,self%ngrids), self%rxcoeff(self%nvertx,self%nrxsmx,self%ngrids), &
         self%wcoeff(self%nedgex,self%nwsmx,self%ngrids), self%tcoeff(self%nfacex,self%ntsmx,self%ntsmx,self%ngrids), &
         self%jlump(self%nvertx,self%ngrids), self%mlump(self%nedgex,self%ngrids), self%hlump(self%nedgex,self%ngrids))

! Read the operator stencils and coefficients
READ(changrid) (((self%lsten(if0,ix,igrid),          &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%nlsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%msten(ie0,ix,igrid),          &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, self%nmsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%jsten(iv0,ix,igrid),          &
                     iv0 = 1, self%nvert(igrid)),    &
                     ix = 1, self%njsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%hsten(ie0,ix,igrid),          &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, self%nhsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%rsten(if0,ix,igrid),          &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%nrsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%rxsten(iv0,ix,igrid),         &
                     iv0 = 1, self%nvert(igrid)),    &
                     ix = 1, self%nrxsmx),           &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%wsten(ie0,ix,igrid),          &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, self%nwsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%tsten(if0,ix,igrid),          &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%ntsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%lmass(if0,ix,igrid),          &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%nlsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%mmass(ie0,ix,igrid),          &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, self%nmsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%jstar(iv0,ix,igrid),          &
                     iv0 = 1, self%nvert(igrid)),    &
                     ix = 1, self%njsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%hstar(ie0,ix,igrid),          &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, self%nhsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%rcoeff(if0,ix,igrid),         &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%nrsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%rxcoeff(iv0,ix,igrid),        &
                     iv0 = 1, self%nvert(igrid)),    &
                     ix = 1, self%nrxsmx),           &
                     igrid = 1, self%ngrids)
READ(changrid) (((self%wcoeff(ie0,ix,igrid),         &
                     ie0 = 1, self%nedge(igrid)),    &
                     ix = 1, self%nwsmx),            &
                     igrid = 1, self%ngrids)
READ(changrid) ((((self%tcoeff(if0,ix,ixx,igrid),    &
                     if0 = 1, self%nface(igrid)),    &
                     ix = 1, self%ntsmx),            &
                     ixx = 1, self%ntsmx),           &
                     igrid = 1, self%ngrids)

! Dimensionalize
self%lmass = self%lmass/(rearth*rearth)

! Construct the tables self%eoffin and self%eofvin
allocate(self%eoffin(self%nfacex,self%nefmx,self%ngrids), self%eofvin(self%nvertx,self%nevmx,self%ngrids))
do igrid = 1, self%ngrids

  do if1 = 1, self%nface(igrid)
    do ix = 1, self%neoff(if1,igrid)
      ie1 = self%eoff(if1,ix,igrid)
      if2 = self%fnxte(ie1,1,igrid)
      IF (if1 == if2) THEN
        ! This edge points out of face if1
	self%eoffin(if1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points into face if1
	self%eoffin(if1,ix,igrid) = 1.0d0
      ENDIF
    enddo
  enddo

  do iv1 = 1, self%nvert(igrid)
    do ix = 1, self%neofv(iv1,igrid)
      ie1 = self%eofv(iv1,ix,igrid)
      iv2 = self%vofe(ie1,1,igrid)
      IF (iv1 == iv2) THEN
        ! This edge points away from vertex iv1
	self%eofvin(iv1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points towards vertex iv1
	self%eofvin(iv1,ix,igrid) = 1.0d0
      ENDIF
    enddo
  enddo

enddo


! Allocate array for size of restriction stencil
allocate(self%ninj(self%nfacex,self%ngrids-1))

! Read the size of the restriction stencil on each grid
self%ninj = 0
READ(changrid) ((self%ninj(if0,igrid),              &
                    if0 = 1, self%nface(igrid)),    &
                    igrid = 1, self%ngrids-1)

! Find maximum value in order to allocate subsequent arrays
self%ninjmx = MAXVAL(self%ninj)

! Allocate arrays for restriction stencils and weights
allocate(self%injsten(self%nfacex,self%ninjmx,self%ngrids-1))
allocate(self%injwgt(self%nfacex,self%ninjmx,self%ngrids-1))

! Read the restriction stencil and weights
READ(changrid) (((self%injsten(if0,ix,igrid),       &
                    if0 = 1, self%nface(igrid)),    &
                    ix = 1, self%ninjmx),           &
                    igrid = 1, self%ngrids-1)
READ(changrid) (((self%injwgt(if0,ix,igrid),        &
                    if0 = 1, self%nface(igrid)),    &
                    ix = 1, self%ninjmx),           &
                    igrid = 1, self%ngrids-1)

end subroutine readgrid

! --------------------------------------------------------------------------------------------------

subroutine centroid(self,if0,long,lat,igrid)

! Find the centroid of cell if0 on grid igrid

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: if0, igrid
real*8,       intent(out) :: long, lat

integer :: ixe, ie1, iv1, iv2
real*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          xc, yc, zc, a, aby3, mag

! Coordinates of `centre' of face (i.e. dual vertex)
long1 = self%flong(if0,igrid)
lat1 = self%flat(if0,igrid)
call ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the face
! Hence find area of face and centroid
xc = 0.0d0
yc = 0.0d0
zc = 0.0d0
do ixe = 1, self%neoff(if0,igrid)
  ie1 = self%eoff(if0,ixe,igrid)
  iv1 = self%vofe(ie1,1,igrid)
  iv2 = self%vofe(ie1,2,igrid)
  long1 = self%vlong(iv1,igrid)
  lat1 = self%vlat(iv1,igrid)
  call ll2xyz(long1,lat1,x1,y1,z1)
  long1 = self%vlong(iv2,igrid)
  lat1 = self%vlat(iv2,igrid)
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

subroutine dual_centroid(self,iv0,long,lat,igrid)

! Find the centroid of dual cell iv0 on grid igrid

implicit none
class(femps), intent(in)  :: self
integer,      intent(in)  :: iv0, igrid
real*8,       intent(out) :: long, lat

integer :: ixe, ie1, iv1, iv2
real*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          xc, yc, zc, a, aby3, mag

! Coordinates of `centre' of dual cell (i.e. vertex)
long1 = self%vlong(iv0,igrid)
lat1 = self%vlat(iv0,igrid)
call ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the dual cell
! Hence find area of dual cell and centroid
xc = 0.0d0
yc = 0.0d0
zc = 0.0d0
do ixe = 1, self%neofv(iv0,igrid)
  ie1 = self%eofv(iv0,ixe,igrid)
  iv1 = self%fnxte(ie1,1,igrid)
  iv2 = self%fnxte(ie1,2,igrid)
  long1 = self%flong(iv1,igrid)
  lat1 = self%flat(iv1,igrid)
  call ll2xyz(long1,lat1,x1,y1,z1)
  long1 = self%flong(iv2,igrid)
  lat1 = self%flat(iv2,igrid)
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

subroutine testpoisson(self)

! To test the solution of thePoisson problem

implicit none
class(femps), intent(inout) :: self

! Number of passes
integer :: npass = 10
integer :: nf, ne, nv, if1, iv1, ipass, nprt
real*8, allocatable :: psi0(:), zeta(:), psi(:), ff1(:), ff2(:), ff3(:), ff4(:), temp1(:), temp2(:)
real*8 :: long, lat, psibar

nf = self%nface(self%ngrids)
ne = self%nedge(self%ngrids)
nv = self%nvert(self%ngrids)

print *,' '
print *,'--------------------------'
print *,' '
print *,'Testing mgsolve '
print *,' '

! Number of values to print for testing
nprt = 5

allocate(psi0(nf),zeta(nf),psi(nf),ff1(nf),ff2(nf),ff3(nf),ff4(nf),temp1(ne),temp2(ne))


! Build coefficients used in Laplacian operator on all grids
call self%buildlap


! Set up test data
! Large-scale part
do if1 = 1, nf
  long = self%flong(if1,self%ngrids)
  lat = self%flat(if1,self%ngrids)
  ! psi0(if1) = SIN(lat)
  psi0(if1) = COS(lat)*SIN(long)
enddo
! Plus small-scale part
psi0(10) = 10.0d0*psi0(10)
! Remove global mean (to ensure unambiguous result)
psibar = SUM(psi0*self%farea(:,self%ngrids))/SUM(self%farea(:,self%ngrids))
psi0 = psi0 - psibar
! Convert to area integrals
ff1 = psi0*self%farea(:,self%ngrids)
print *,'Original field ff1 =     ',ff1(1:nprt)
print *,' '

! Calculate laplacian
call self%laplace(ff1,zeta,self%ngrids,nf,ne)


! Now invert laplacian to check we get back to where we started

! Initialize result to zero
! Note psi will be stream function times grid cell area
psi = 0.0d0
temp2 = 0.0d0

! Iterate several passes
do ipass = 1, npass

  print *,'Pass ',ipass

  ! Compute residual based on latest estimate
  call self%massL(psi,ff2,self%ngrids,nf)

  call self%Ddual1(ff2,temp1,self%ngrids,nf,ne)

  ! Improve the estimate temp2 that we obtained in the previous pass
  call self%massMinv(temp1,temp2,self%ngrids,ne,4)

  call self%Dprimal2(temp2,ff3,self%ngrids,ne,nf)

  ! Residual
  ff4 = zeta - ff3

  ! Now solve the approximate Poisson equation
  ! D2 self%xMinv D1bar L psi' = residual
  call self%mgsolve(ff3,ff4,self%ngrids)

  ! And increment best estimate
  ! *** We could think about adding beta*ff3 here and tuning beta for optimal convergence ***
  psi = psi + ff3

  ! Remove global mean (to ensure unambiguous result)
  psibar = SUM(psi)/SUM(self%farea(:,self%ngrids))
  psi = psi - psibar*self%farea(:,self%ngrids)

  print *,'Original field ff1      = ',ff1(1:nprt)
  print *,'Soln of Poisson eqn psi = ', psi(1:nprt)
  ff4 = (ff1-psi)/self%farea(:,self%ngrids)
  print *,'RMS err in global problem = ',sqrt(sum(ff4*ff4)/nf)
  print *,' '

enddo

end subroutine testpoisson

! --------------------------------------------------------------------------------------------------

end module femps_mod
