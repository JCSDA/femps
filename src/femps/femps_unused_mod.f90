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

subroutine Ddual2(grid,pbops,f,df,igrid,ne,nv)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises integrals along
! dual edges; df comprises integrals of the derivative
! over dual cells.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne, nv
real(kind=kind_real),       intent(in)  :: f(ne)
real(kind=kind_real),       intent(out) :: df(nv)

integer :: iv1, ix, ie1
real(kind=kind_real) :: temp

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
real(kind=kind_real),       intent(in)    :: f1(nf)
real(kind=kind_real),       intent(inout) :: f2(nf)

integer :: if1, iter, miter
real(kind=kind_real) :: temp(nf)

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
real(kind=kind_real),       intent(in)    :: f1(nv)
real(kind=kind_real),       intent(inout) :: f2(nv)

integer :: iv1, iter, miter
real(kind=kind_real) :: temp(nv)
real(kind=kind_real) :: relax = 1.4d0 ! relax = 1.4 is good for ijlump = 3 on hex and cube grids

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

subroutine operT(grid,pbops,f1,f2,igrid,ne,nf)

! Apply the T operator:
! compute cell integrals of 2 x kinetic energy from edge integrals
! of normal fluxes

implicit none
type(fempsgrid), intent(in) :: grid
type(fempspbops), intent(in) :: pbops
integer,      intent(in)  :: igrid, ne, nf
real(kind=kind_real),       intent(in)  :: f1(ne)
real(kind=kind_real),       intent(out) :: f2(nf)

integer :: if1, ix1, ix2, ne1, ie1, ie2
real(kind=kind_real) :: temp

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
