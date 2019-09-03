MODULE grid

IMPLICIT NONE

! COUNTING

! Number of grids in multigrid hierarchy
INTEGER :: ngrids

! nface        Number of faces on each grid
! nedge        Number of edges on each grid
! nvert        Number of vertices on each edge 
INTEGER, ALLOCATABLE :: nface(:), nedge(:), nvert(:)
INTEGER :: nfacex, nedgex, nvertx

! neoff        Number of edges and vertices of each face on each grid
! neofv        Number of edges of each vertex on each grid
INTEGER, ALLOCATABLE :: neoff(:,:), neofv(:,:)
INTEGER :: nefmx, nevmx


! CONNECTIVITY

! fnxtf        Faces next to each face on each grid
! eoff         Edges of each face on each grid
! voff         Vertices of each face on each grid
! fnxte        Faces either side of each edge on each grid
! vofe         Vertices at the ends of each edge on each grid
! fofv         Faces around each vertex on each grid
! eofv         Edges incident on each vertex on each grid
INTEGER, ALLOCATABLE :: fnxtf(:,:,:), eoff(:,:,:), voff(:,:,:), &
                        fnxte(:,:,:), vofe(:,:,:), &
                        fofv(:,:,:), eofv(:,:,:)

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


! eoffin(f,j,:)   Indicates whether the normal at the j'th edge is
!                 inward or outward relative to face f.
! eofvin(v,j,:)   Indicates whether the tangent at the j'th edge is
!                 inward or outward relative to vertex v.
INTEGER, ALLOCATABLE :: eoffin(:,:,:), eofvin(:,:,:)


! COORDINATES AND GEOMETRICAL INFORMATION

! flong        Longitude of faces on each grid
! flat         Latitude of faces on each grid
! vlong        Longitude of vertices on each grid
! vlat         Latitude of vertices on each grid
! farea        Area of faces on each grid
! varea        Area of dual faces on each grid
! ldist        Primal edge length, i.e. distance between neighbouring vertices
! ddist        Dual edge length, i.e. distance between neighbouring face centres
! fareamin     Minimum face areaon each grid
REAL*8, ALLOCATABLE :: flong(:,:), flat(:,:), &
                       vlong(:,:), vlat(:,:), &
                       farea(:,:), varea(:,:), &
                       ldist(:,:), ddist(:,:), &
                       fareamin(:)

! Conventions
!
! Latitude and longitude in radians.
! Area and lengths are for the unit sphere.


! HODGE STAR, MASS MATRIX, AND RELATED OPERATORS

! nlsten        Number of faces in stencil for L mass matrix
! lsten         Stencil for L mass matrix
! lmass         Coefficients for L mass matrix
! nmsten        Number of faces in stencil for M mass matrix
! msten         Stencil for M mass matrix
! mmass         Coefficients for M mass matrix
! njsten        Number of vertices in stencil for J operator
! jsten         Stencil for J operator
! jstar         Coefficients for J operator
! nhsten        Number of edges in stencil for H operator
! hsten         Stencil for H operator
! hstar         Coefficients for H operator
! nrsten        Number of vertices in stencil for R operator (= neoff)
! rsten         Stencil for R operator (= voff)
! rcoeff        Coefficients for R operator
! nrxsten       Number of faces in stencil for R transpose operator (= neofv)
! rxsten        Stencil for R transpose operator (= fofv)
! rxcoeff       Coefficients for R transpose operator
! nwsten        Number of edges in stencil for W operator
! wsten         Stencil for W operator
! wcoeff        Coefficients for W operator
! ntsten        Number of edges in stencel for T operator
! tsten         Stencil for T operator
! tcoeff        Coefficients for T operator
! jlump         Coefficients of lumped J matrix
! mlump         Coefficients of lumped M matrix
! hlump         Coefficients of lumped H matrix
! nxminvten     Number of edges in stencil for approximate inverse of M
! xminvsten     Stencil for approximate inverse of M
! xminv         Coefficients for approximate inverse of M
! velcoeff      Coefficients to reconstruct velocity vector in cells
INTEGER, ALLOCATABLE :: nlsten(:,:), nmsten(:,:), njsten(:,:), &
                        nhsten(:,:), nrsten(:,:), nrxsten(:,:), &
                        nwsten(:,:), ntsten(:,:), nxminvsten(:,:)
INTEGER, ALLOCATABLE :: lsten(:,:,:), msten(:,:,:), jsten(:,:,:), &
                        hsten(:,:,:), rsten(:,:,:), rxsten(:,:,:), &
                        wsten(:,:,:), tsten(:,:,:), xminvsten(:,:,:)
REAL*8, ALLOCATABLE :: lmass(:,:,:), mmass(:,:,:), jstar(:,:,:), &
                       hstar(:,:,:), rcoeff(:,:,:), rxcoeff(:,:,:), &
                       wcoeff(:,:,:), tcoeff(:,:,:,:), jlump(:,:), &
                       mlump(:,:), hlump(:,:), xminv(:,:,:), &
                       velcoeff(:,:,:)
INTEGER :: nlsmx, nmsmx, njsmx, nhsmx, nrsmx, nrxsmx, nwsmx, ntsmx, nxmisx


! RESTRICTION AND PROLONGATION OPERATORS FOR MULTIGRID

! ninj           Number of faces in stencil for restriction operator
! injsten        Stencil for restriction operator
! injwgt         Weights for restriction operator
INTEGER, ALLOCATABLE :: ninj(:,:), injsten(:,:,:)
REAL*8, ALLOCATABLE :: injwgt(:,:,:)
INTEGER :: ninjmx


END MODULE grid

! ================================================================

MODULE runtype

! Default values are set here; they may be overridden
! via the namelist

IMPLICIT NONE

! File containing grid information
! CHARACTER*31 :: ygridfile = 'gridopermap_hex_0000010242.dat'
CHARACTER*31 :: ygridfile = 'gridopermap_cube_0000013824.dat'


END MODULE runtype

! ================================================================

MODULE constants

! Various physical and geometrical constants

IMPLICIT NONE

! Pi
REAL*8, PARAMETER :: pi = 3.14159265358979323d0, piby2 = 0.5d0*pi

! Earth's radius
REAL*8 :: rearth = 6371220.0d0


END MODULE constants

! ================================================================

MODULE laplacecoeff

! Reference state values on the grid hierarchy needed for
! multigrid poisson solver

USE grid
IMPLICIT NONE

! lapdiag        Diagonal coefficient of Laplacian operator
!                on all grids
! underrel       Under-relaxation parameter on all grids
REAL*8, ALLOCATABLE :: lapdiag(:,:), underrel(:)


END MODULE laplacecoeff

! ================================================================

MODULE channels

! Tidy list of all I/O channels in one place to avoid accidental
! overuse of any channel number

IMPLICIT NONE

INTEGER, PARAMETER :: channml = 20          ! For reading namelists
INTEGER, PARAMETER :: changrid = 25         ! Grid information

END MODULE channels

! ================================================================

PROGRAM fempoisson

IMPLICIT NONE

! ----------------------------------------------------------------

CALL preliminary
print *,'Done preliminary'

! Test the Poisson solver
CALL testpoisson
print *,'Done testpoisson'


! ----------------------------------------------------------------

END PROGRAM fempoisson

! ================================================================

SUBROUTINE preliminary

! Preliminary calculations and setting up

USE constants
USE grid
USE runtype

IMPLICIT NONE

! ----------------------------------------------------------------

! Read namelist information
CALL readnml
print *,'Namelists read'

! ----------------------------------------------------------------

! Read in the grid data
CALL readgrid
print *,'Done readgrid'

! ----------------------------------------------------------------

! Allocate array space now that resolution is known
CALL allocateall
print *,'Done allocateall'

! ----------------------------------------------------------------

! Build a lumped version of the J, M and H matrices
CALL buildjlump
CALL buildmlump
CALL buildhlump

! And build a local operator to approximate the inverse of M
CALL buildxminv

PRINT *,'Lumped J, M and H and approximate M-inverse created'

! ----------------------------------------------------------------

END SUBROUTINE preliminary

! ================================================================

SUBROUTINE readnml

! Read namelists

USE runtype
USE channels
IMPLICIT NONE


NAMELIST /rundata/ ygridfile

OPEN(channml,FILE='poissonnml.in',DELIM='APOSTROPHE')
READ(channml,rundata)
CLOSE(channml)


END SUBROUTINE readnml

! ================================================================

SUBROUTINE allocateall

! Allocate array space that will be needed throughout the code

USE constants
USE laplacecoeff

IMPLICIT NONE

! ----------------------------------------------------------------

! Arrays in module laplacian
ALLOCATE(lapdiag(nfacex,ngrids))
ALLOCATE(underrel(ngrids))

! ----------------------------------------------------------------

END SUBROUTINE allocateall

! ================================================================

SUBROUTINE buildlap

! Extract the diagonal coefficients of the Laplacian operator needed
! for the multigrid Poisson solver at all resolutions

USE grid
USE laplacecoeff
USE constants

IMPLICIT NONE
INTEGER :: igrid, if1, if2, if3, ie1, ie2, &
           ixd2, ixm, ixd1, ixl
REAL*8 :: temp, temp1(nfacex), cd2, cm, cd1, cl

! ---------------------------------------------------------------

! Under-relaxation parameter.
! *** There is scope to optimize here ***
DO igrid = 1, ngrids
  IF (nefmx == 6) THEN
    underrel(igrid) = 0.8d0
  ELSEIF (nefmx == 4) THEN
    underrel(igrid) = 0.8d0
  ELSE
    print *,'Choose a sensible value for underrel in buildlap'
    STOP
  ENDIF
ENDDO


! Extract diagonal coefficient of Laplacian operator
DO igrid = 1, ngrids
  ! Loop over cells
  DO if1 = 1, nface(igrid)

    temp = 0.0d0
    ! Loop over edges of if1 involved in Dprimal2 operator
    DO ixd2 = 1, neoff(if1,igrid)
      ie1 = eoff(if1,ixd2,igrid)
      cd2 = -eoffin(if1,ixd2,igrid)
      ! Loop over edges involved in approximate M-inverse operator
      DO ixm = 1, nxminvsten(ie1,igrid)
        ie2 = xminvsten(ie1,ixm,igrid)
        cm = xminv(ie1,ixm,igrid)
        ! Loop over cells involved in Ddual1 operator
        cd1 = 1.0d0
        DO ixd1 = 1, 2
          if2 = fnxte(ie2,ixd1,igrid)
          cd1 = -cd1
          ! Loop over cells in L operator
          DO ixl = 1, nlsten(if2,igrid)
            if3 = lsten(if2,ixl,igrid)
            cl = lmass(if2,ixl,igrid)
            IF (if3 == if1) temp = temp + cd2*cm*cd1*cl
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    lapdiag(if1,igrid) = temp

  ENDDO
ENDDO


! ---------------------------------------------------------------

END SUBROUTINE buildlap

! ================================================================

SUBROUTINE Dprimal1(f,df,igrid,nv,ne)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises pointwise values
! at vertices; df comprises integrals of the derivative
! of f along primal cell edges.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nv, ne
REAL*8, INTENT(IN) :: f(nv)
REAL*8, INTENT(OUT) :: df(ne)
INTEGER :: ie1, iv1, iv2

! ----------------------------------------------------------------

DO ie1 = 1, ne
  iv1 = vofe(ie1,1,igrid)
  iv2 = vofe(ie1,2,igrid)
  df(ie1) = f(iv2) - f(iv1)
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Dprimal1

! ================================================================

SUBROUTINE Dprimal2(f,df,igrid,ne,nf)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises integrals along
! primal edges; df comprises integrals of the derivative
! over primal cells.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne, nf
REAL*8, INTENT(IN) :: f(ne)
REAL*8, INTENT(OUT) :: df(nf)
INTEGER :: if1, ix, ie1
REAL*8 :: temp

! ----------------------------------------------------------------

DO if1 = 1, nf
  temp = 0.0d0
  DO ix = 1, neoff(if1,igrid)
    ie1 = eoff(if1,ix,igrid)
    temp = temp - f(ie1)*eoffin(if1,ix,igrid)
  ENDDO
  df(if1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Dprimal2

! ================================================================

SUBROUTINE Ddual1(f,df,igrid,nf,ne)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises pointwise values
! at face centres; df comprises integrals of the derivative
! of f along dual cell edges.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
REAL*8, INTENT(IN) :: f(nf)
REAL*8, INTENT(OUT) :: df(ne)
INTEGER :: ie1, if1, if2

! ----------------------------------------------------------------

DO ie1 = 1, ne
  if1 = fnxte(ie1,1,igrid)
  if2 = fnxte(ie1,2,igrid)
  df(ie1) = f(if2) - f(if1)
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Ddual1

! ================================================================

SUBROUTINE Ddual2(f,df,igrid,ne,nv)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises integrals along
! dual edges; df comprises integrals of the derivative
! over dual cells.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne, nv
REAL*8, INTENT(IN) :: f(ne)
REAL*8, INTENT(OUT) :: df(nv)
INTEGER :: iv1, ix, ie1
REAL*8 :: temp

! ----------------------------------------------------------------

DO iv1 = 1, nv
  temp = 0.0d0
  DO ix = 1, neofv(iv1,igrid)
    ie1 = eofv(iv1,ix,igrid)
    temp = temp + f(ie1)*eofvin(iv1,ix,igrid)
  ENDDO
  df(iv1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE Ddual2

! ================================================================

SUBROUTINE buildjlump

! Build a lumped version of the j matrix

USE grid
IMPLICIT NONE

! Two choices for lumped J
! ijlump = 1     Diagonal of J
! ijlump = 2     Exact answer when applied to constant input vector
! ijlump = 3     Exact answer when applied to the discrete representation
!                  of a constant scalar (input vector proportional to varea)

INTEGER, PARAMETER :: ijlump = 3
INTEGER :: igrid, iv1, nv
REAL*8, ALLOCATABLE :: p1(:), jp1(:)

! ----------------------------------------------------------------

DO igrid = 1, ngrids

  IF (ijlump == 1) THEN

    DO iv1 = 1, nvert(igrid)
      jlump(iv1,igrid) = jstar(iv1,1,igrid)
    ENDDO

  ELSEIF (ijlump == 2) THEN

    nv = nvert(igrid)
    ALLOCATE(p1(nv), jp1(nv))
    p1 = 1.0d0
    CALL HodgeJ(p1,jp1,igrid,nv)
    jlump(1:nv,igrid) = jp1
    DEALLOCATE(p1,jp1)

  ELSEIF (ijlump == 3) THEN

    nv = nvert(igrid)
    ALLOCATE(p1(nv), jp1(nv))
    p1 = varea(1:nv,igrid)
    CALL HodgeJ(p1,jp1,igrid,nv)
    jlump(1:nv,igrid) = jp1/varea(1:nv,igrid)
    DEALLOCATE(p1,jp1)

  ELSE

    PRINT *,'Option ijlump = ',ijlump,' not available in subroutine buildjlump'
    STOP

  ENDIF

ENDDO

! ----------------------------------------------------------------

END SUBROUTINE buildjlump

! ================================================================

SUBROUTINE buildmlump

! Build a lumped version of the M matrix

USE grid
IMPLICIT NONE

! Three choices for lumped M
! imlump = 1     Diagonal of M
! imlump = 2     Row sum of absolute values of M
! imlump = 3     Best fit to solid body rotation normal to each edge

INTEGER, PARAMETER :: imlump = 3
INTEGER :: igrid, ie1, ie2, ix, iv0, ne, nv, iv1, iv2
REAL*8 :: temp, slon, slat, clon, clat, a1, a2, a3, num, den, &
          b1, b2, b3, c1, c2, c3, x, y, z
REAL*8, ALLOCATABLE :: psi1(:), psi2(:), psi3(:),    &
                       v1(:), v2(:), v3(:), vv(:),   &
                       mv1(:), mv2(:), mv3(:)

! ----------------------------------------------------------------

DO igrid = 1, ngrids

  IF (imlump == 1) THEN

    DO ie1 = 1, nedge(igrid)
      mlump(ie1,igrid) = mmass(ie1,1,igrid)
    ENDDO

  ELSEIF (imlump == 2) THEN

    DO ie1 = 1, nedge(igrid)
      temp = 0.0d0
      DO ix = 1, nmsten(ie1,igrid)
        ie2 = msten(ie1,ix,igrid)
        temp = temp + ABS(mmass(ie1,ix,igrid))
      ENDDO
      mlump(ie1,igrid) = temp
    ENDDO

  ELSEIF (imlump == 3) THEN

    nv = nvert(igrid)
    ne = nedge(igrid)
    ALLOCATE(psi1(nv), psi2(nv), psi3(nv),    &
             v1(ne), v2(ne), v3(ne), vv(ne),  &
             mv1(ne), mv2(ne), mv3(ne))
    ! Set up three solid body rotation fields
    DO 	iv0 = 1, nv
      slon = SIN(vlong(iv0,igrid))
      clon = COS(vlong(iv0,igrid))
      slat = SIN(vlat(iv0,igrid))
      clat = COS(vlat(iv0,igrid))
      psi1(iv0) = clat*clon
      psi2(iv0) = clat*slon
      psi3(iv0) = slat
    ENDDO
    CALL Dprimal1(psi1,v1,igrid,nv,ne)
    CALL Dprimal1(psi2,v2,igrid,nv,ne)
    CALL Dprimal1(psi3,v3,igrid,nv,ne)
    CALL massM(v1,mv1,igrid,ne)
    CALL massM(v2,mv2,igrid,ne)
    CALL massM(v3,mv3,igrid,ne)
    ! Now loop over edges
    DO ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      den = SQRT(a1*a1 + a2*a2 + a3*a3)
      a1 = a1/den
      a2 = a2/den
      a3 = a3/den
      ! Demand that lumped matrix agrees with full M for this v
      num = a1*mv1(ie1) + a2*mv2(ie1) + a3*mv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      mlump(ie1,igrid) = num/den
    ENDDO
    DEALLOCATE(psi1,psi2,psi3,v1,v2,v3,vv,mv1,mv2,mv3)

  ELSE

    PRINT *,'Option imlump = ',imlump,' not available in subroutine buildmlump'
    STOP

  ENDIF

ENDDO

! ----------------------------------------------------------------

END SUBROUTINE buildmlump

! ================================================================

SUBROUTINE buildhlump

! Build a lumped version of the H matrix

USE grid
IMPLICIT NONE

! Three choices for lumped H
! ihlump = 1     Diagonal of H
! ihlump = 2     Row sum of absolute values of H
! ihlump = 3     Best fit to global divergent flow on dual grid
! ihlump = 4     Best fit to solid body rotation normal to each edge

INTEGER, PARAMETER :: ihlump = 3
INTEGER :: igrid, ie1, ie2, ix, iv0, if0, ne, nv, nf
REAL*8 :: temp, slon, slat, clon, clat, a1, a2, a3, num, den
REAL*8, ALLOCATABLE :: psi1(:), psi2(:), psi3(:),  &
                       v1(:), v2(:), v3(:),        &
                       hv1(:), hv2(:), hv3(:)

! ----------------------------------------------------------------

DO igrid = 1, ngrids

  IF (ihlump == 1) THEN

    DO ie1 = 1, nedge(igrid)
      hlump(ie1,igrid) = hstar(ie1,1,igrid)
    ENDDO

  ELSEIF (ihlump == 2) THEN

    DO ie1 = 1, nedge(igrid)
      temp = 0.0d0
      DO ix = 1, nhsten(ie1,igrid)
        ie2 = hsten(ie1,ix,igrid)
        temp = temp + ABS(hstar(ie1,ix,igrid))
      ENDDO
      hlump(ie1,igrid) = temp
    ENDDO

  ELSEIF (ihlump == 3) THEN

    nf = nface(igrid)
    ne = nedge(igrid)
    ALLOCATE(psi1(nf), psi2(nf), psi3(nf),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three global divergent fields
    DO 	if0 = 1, nf
      slon = SIN(flong(if0,igrid))
      clon = COS(flong(if0,igrid))
      slat = SIN(flat(if0,igrid))
      clat = COS(flat(if0,igrid))
      psi1(if0) = slat
      psi2(if0) = clat*slon
      psi3(if0) = clat*clon
    ENDDO
    CALL Ddual1(psi1,v1,igrid,nf,ne)
    CALL Ddual1(psi2,v2,igrid,nf,ne)
    CALL Ddual1(psi3,v3,igrid,nf,ne)
    CALL HodgeH(v1,hv1,igrid,ne)
    CALL HodgeH(v2,hv2,igrid,ne)
    CALL HodgeH(v3,hv3,igrid,ne)
    ! Now loop over edges
    DO ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      ! Demand that lumped matrix agrees with full H for this v
      num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      hlump(ie1,igrid) = num/den
    ENDDO
    DEALLOCATE(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  ELSEIF (ihlump == 4) THEN

    nv = nvert(igrid)
    ne = nedge(igrid)
    ALLOCATE(psi1(nv), psi2(nv), psi3(nv),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three solid body rotation fields
    DO 	iv0 = 1, nv
      slon = SIN(vlong(iv0,igrid))
      clon = COS(vlong(iv0,igrid))
      slat = SIN(vlat(iv0,igrid))
      clat = COS(vlat(iv0,igrid))
      psi1(iv0) = slat
      psi2(iv0) = clat*slon
      psi3(iv0) = clat*clon
    ENDDO
    CALL Dprimal1(psi1,v1,igrid,nv,ne)
    CALL Dprimal1(psi2,v2,igrid,nv,ne)
    CALL Dprimal1(psi3,v3,igrid,nv,ne)
    CALL HodgeH(v1,hv1,igrid,ne)
    CALL HodgeH(v2,hv2,igrid,ne)
    CALL HodgeH(v3,hv3,igrid,ne)
    ! Now loop over edges
    DO ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      ! Demand that lumped matrix agrees with full H for this v
      num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      hlump(ie1,igrid) = num/den
    ENDDO
    DEALLOCATE(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  ELSE

    PRINT *,'Option ihlump = ',ihlump,' not available in subroutine buildhlump'
    STOP

  ENDIF

ENDDO

! ----------------------------------------------------------------

END SUBROUTINE buildhlump

! ================================================================

SUBROUTINE buildxminv

! Determine stencil and coefficients for a local approximation
! to the inverse of M.
!
! Two options are coded. The first is simply the inverse of
! the diagonal mass lumped approximation to M. The second is
! based on a single underrelaxed Jacobi iteration to the
! inverse of M.

USE grid
IMPLICIT NONE

LOGICAL :: llump
INTEGER :: igrid, ie0, ix, ie1
REAL*8 :: temp, diag
REAL*8 :: relax

! ----------------------------------------------------------------

! Underrelaxation coefficient depends on grid
! 0.9 and 1.4 when using mlump in Jacobi
! Check for consistency with massMinv
IF (nefmx == 4) THEN
  ! Use sparse approximate inverse on cubed sphere
  llump = .false.
  relax = 0.9d0
ELSE
  ! Use diagonal approximate inverse on hexagonal-icosahedral grid
  llump = .false.
  relax = 1.4d0
ENDIF


IF (llump) THEN

  ! Diagonal operator: stencil size is 1
  nxmisx = 1
  ALLOCATE(nxminvsten(nedgex,ngrids), xminvsten(nedgex,nxmisx,ngrids), &
           xminv(nedgex,nxmisx,ngrids))
  DO igrid = 1, ngrids
    DO ie0 = 1, nedge(igrid)
      ! Stencil for edge ie0 is ie0 itself, and coeff is
      ! inverse of the diagonal term of the lumped matrix
      nxminvsten(ie0,igrid) = 1
      xminvsten(ie0,1,igrid) = ie0
      xminv(ie0,1,igrid) = 1.0d0/mlump(ie0,igrid)
    ENDDO
  ENDDO

ELSE

  ! Stencil is the same as for M itself
  nxmisx = nmsmx
  ALLOCATE(nxminvsten(nedgex,ngrids), xminvsten(nedgex,nxmisx,ngrids), &
           xminv(nedgex,nxmisx,ngrids))
  DO igrid = 1, ngrids
    DO ie0 = 1, nedge(igrid)
      ! Stencil for edge ie0 is the same as the stencil for M
      nxminvsten(ie0,igrid) = nmsten(ie0,igrid)
      DO ix = 1, nmsten(ie0,igrid)
        ie1 = msten(ie0,ix,igrid)
        xminvsten(ie0,ix,igrid) = ie1
        IF (ie1 == ie0) THEN
          diag = 1.0d0 + relax
        ELSE
          diag = 0.0d0
        ENDIF
        temp = mmass(ie0,ix,igrid)/mlump(ie1,igrid)
        xminv(ie0,ix,igrid) = (diag - relax*temp)/mlump(ie0,igrid)
      ENDDO
    ENDDO
  ENDDO

ENDIF


! ----------------------------------------------------------------

END SUBROUTINE buildxminv

! ================================================================

SUBROUTINE massL(f1,f2,igrid,nf)

! Apply the mass matrix L to field f1 to obtain the result f2

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(OUT) :: f2(nf)
INTEGER :: if1, if2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO if1 = 1, nf
  temp = 0.0d0
  DO ix = 1, nlsten(if1,igrid)
    if2 = lsten(if1,ix,igrid)
    temp = temp + f1(if2)*lmass(if1,ix,igrid)
  ENDDO
  f2(if1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE massL

! ================================================================

SUBROUTINE massM(f1,f2,igrid,ne)

! Apply the mass matrix M to field f1 to obtain field f2

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: ie1, ie2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO ie1 = 1, ne
  temp = 0.0d0
  DO ix = 1, nmsten(ie1,igrid)
    ie2 = msten(ie1,ix,igrid)
    temp = temp + f1(ie2)*mmass(ie1,ix,igrid)
  ENDDO
  f2(ie1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE massM

! ================================================================

SUBROUTINE approxMinv(f1,f2,igrid,ne)

! Apply an approximate inverse of the mass matrix M
! to field f1 to obtain field f2

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: ie1, ie2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO ie1 = 1, ne
  temp = 0.0d0
  DO ix = 1, nxminvsten(ie1,igrid)
    ie2 = xminvsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*xminv(ie1,ix,igrid)
  ENDDO
  f2(ie1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE approxMinv

! ================================================================

SUBROUTINE HodgeJ(f1,f2,igrid,nv)

! Apply the Hodge star J operator that converts dual face
! integrals f1 to vertex values f2 on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nv
REAL*8, INTENT(IN) :: f1(nv)
REAL*8, INTENT(OUT) :: f2(nv)
INTEGER :: iv1, iv2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO iv1 = 1, nv
  temp = 0.0d0
  DO ix = 1, njsten(iv1,igrid)
    iv2 = jsten(iv1,ix,igrid)
    temp = temp + f1(iv2)*jstar(iv1,ix,igrid)
  ENDDO
  f2(iv1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE HodgeJ

! ================================================================

SUBROUTINE HodgeH(f1,f2,igrid,ne)

! Apply the Hodge star H operator that converts dual edge
! integrals (circulations) f1 to primal edge integrals (fluxes) f2
! on grid igrid.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: ie1, ie2, ix
REAL*8 :: temp

! ----------------------------------------------------------------

DO ie1 = 1, ne
  temp = 0.0d0
  DO ix = 1, nhsten(ie1,igrid)
    ie2 = hsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*hstar(ie1,ix,igrid)
  ENDDO
  f2(ie1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE HodgeH

! ================================================================

SUBROUTINE massLinv(f1,f2,igrid,nf,niter)

! Apply the inverse of the mass matrix L to field f1 to obtain
! the result f2.

! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If L is diagonal then there is no need to iterate.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: niter
INTEGER, INTENT(IN) :: igrid, nf
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(INOUT) :: f2(nf)
INTEGER :: if1, iter, miter
REAL*8 :: temp(nf)

! ----------------------------------------------------------------

miter = ABS(niter)

IF (niter < 0 .OR. nlsmx == 1) THEN
  ! First guess based on diagonal L
  DO if1 = 1, nf
    f2(if1) = f1(if1)/lmass(if1,1,igrid)
  ENDDO
ENDIF

IF (nlsmx > 1) THEN
  ! L is not diagonal, so use Jacobi iteration to invert
  DO iter = 1, miter
    CALL massL(f2,temp,igrid,nf)
    temp = f1 - temp
    DO if1 = 1, nf
      f2(if1) = f2(if1) + temp(if1)/lmass(if1,1,igrid)
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------

END SUBROUTINE massLinv

! ================================================================

SUBROUTINE massMinv(f1,f2,igrid,ne,niter)

! Apply the inverse of the mass matrix M to the field f1
! to obtain the result f2

! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If M is diagonal then there is no need to iterate.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: niter
INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(INOUT) :: f2(ne)
INTEGER :: ie1, iter, miter
REAL*8 :: temp(ne)
REAL*8 :: relax

! ----------------------------------------------------------------

! Underrelaxation coefficient depends on grid
IF (nefmx == 4) THEN
  relax = 0.9d0
ELSE
  relax = 1.4d0
ENDIF

miter = ABS(niter)

IF (niter < 0 .OR. nmsmx == 1) THEN
  ! First guess based on lumped M
  DO ie1 = 1, ne
    f2(ie1) = f1(ie1)/mlump(ie1,igrid)
  ENDDO
ENDIF

IF (nmsmx > 1) THEN
  ! M is not diagonal, so use Jacobi iteration to invert
  DO iter = 1, miter
    CALL massM(f2,temp,igrid,ne)
    temp = relax*(f1 - temp)
    DO ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/mlump(ie1,igrid)
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------

END SUBROUTINE massMinv

! ================================================================

SUBROUTINE HodgeJinv(f1,f2,igrid,nv,niter)

! Apply the inverse Hodge star operator J^{-1} that maps from
! E_p to V_d on grid igrid.
!
! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If J is diagonal then there is no need to iterate.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: niter
INTEGER, INTENT(IN) :: igrid, nv
REAL*8, INTENT(IN) :: f1(nv)
REAL*8, INTENT(INOUT) :: f2(nv)
INTEGER :: iv1, iter, miter
REAL*8 :: temp(nv)
REAL*8 :: relax = 1.4d0 ! relax = 1.4 is good for ijlump = 3 on hex and cube grids 

! ----------------------------------------------------------------

miter = ABS(niter)

IF (niter < 0 .OR. njsmx == 1) THEN
  ! First guess based on lumped J
  DO iv1 = 1, nv
    f2(iv1) = f1(iv1)/jlump(iv1,igrid)
  ENDDO
ENDIF

IF (njsmx > 1) THEN
  ! J is not diagonal, so use Jacobi iteration to invert
  DO iter = 1, miter
    CALL HodgeJ(f2,temp,igrid,nv)
    temp = relax*(f1 - temp)
    DO iv1 = 1, nv
      f2(iv1) = f2(iv1) + temp(iv1)/jlump(iv1,igrid)
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------

END SUBROUTINE HodgeJinv

! ================================================================

SUBROUTINE HodgeHinv(f1,f2,igrid,ne,niter)

! Apply the inverse Hodge star operator H^{-1} that maps from
! S_p to S_d on grid igrid.
!
! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If H is diagonal then there is no need to iterate.

USE grid
IMPLICIT NONE

INTEGER,INTENT(IN) :: niter
INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(INOUT) :: f2(ne)
INTEGER :: ie1, iter, miter
REAL*8 :: temp(ne)
REAL*8 :: relax = 1.4d0 ! relax = 1.4 is good for ihlump = 3 on hex and cube grids 

! ----------------------------------------------------------------

miter = ABS(niter)

IF (niter < 0 .OR. nhsmx == 1) THEN
  ! First guess based on diagonal H
  DO ie1 = 1, ne
    f2(ie1) = f1(ie1)/hlump(ie1,igrid)
  ENDDO
ENDIF

IF (nhsmx > 1) THEN
  ! H is not diagonal, so use Jacobi iteration to invert
  DO iter = 1, miter
    CALL HodgeH(f2,temp,igrid,ne)
    temp = relax*(f1 - temp)
    DO ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/hlump(ie1,igrid)
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------

END SUBROUTINE HodgeHinv

! ================================================================

SUBROUTINE operW_original(f1,f2,igrid,ne)

! Apply the W operator:
! given fluxes f1 across primal edges, construct
! the rotated fluxes across dual edges f2, on grid igrid.

! This is the original formulation, building W from R
! a la TRiSK. It probably requires an MPI reduce operation
! so is likely to be inefficient.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: if1, ie1, ie2, ix, ix1, ix2, ixv, ne1
REAL*8 :: ss, w

! ----------------------------------------------------------------

! Initialize to zero
f2 = 0.0d0

! Loop over faces
DO if1 = 1, nface(igrid)
  ne1 = neoff(if1,igrid)
  ! For each edge of this face
  DO ix1 = 1, ne1
    ss = -0.5
    ie1 = eoff(if1,ix1,igrid)
    ! Find the contribution to f2 from every other
    ! edge of this face
    DO ix = 0, ne1 - 2
      ixv = MOD(ix1 + ix - 1,ne1) + 1
      ix2 = MOD(ix1 + ix,ne1) + 1
      ie2 = eoff(if1,ix2,igrid)
      ss = ss + rcoeff(if1,ixv,igrid)
      w = -ss*eoffin(if1,ix1,igrid)*eoffin(if1,ix2,igrid)
      f2(ie1) = f2(ie1) + w*f1(ie2)
    ENDDO
  ENDDO
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE operW_original

! ================================================================

SUBROUTINE operR_original(f1,f2,igrid,nf,nv)

! Apply the R operator:
! map from V_p to E_p

! This is the original formulation. It loops over `source'
! entities rather than target entities and so will require
! an MPI reduce.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, nv
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(OUT) :: f2(nv)
INTEGER :: if1, iv1, ix1, ne1

! ----------------------------------------------------------------

! Initialize to zero
f2 = 0.0d0

! Loop over faces
DO if1 = 1, nface(igrid)
  ne1 = neoff(if1,igrid)
  ! Share out this face's contributions to its surrounding vertices
  DO ix1 = 1, ne1
    iv1 = rsten(if1,ix1,igrid)
    f2(iv1) = f2(iv1) + f1(if1)*rcoeff(if1,ix1,igrid)
  ENDDO
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE operR_original

! ================================================================

SUBROUTINE operW(f1,f2,igrid,ne)

! Apply the W operator:
! given edge integrals of normal components f1 on primal edges,
! construct edge integrals of normal components of perpendicular
! field f2, on grid igrid.

! This formulation uses pre-build stencil and coefficients to
! avoid the need for MPI reduce. It is mathematically equivalent
! to the original formulation.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(ne)
INTEGER :: ie0, ne1, ix1, ie1
REAL*8 :: temp

! ----------------------------------------------------------------

! Loop over vertices
DO ie0 = 1, nedge(igrid)
  ne1 = nwsten(ie0,igrid)
  ! Collect contributions from stencil
  temp = 0.0d0
  DO ix1 = 1, ne1
    ie1 = wsten(ie0,ix1,igrid)
    temp = temp + f1(ie1)*wcoeff(ie0,ix1,igrid)
  ENDDO
  f2(ie0) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE operW

! ================================================================

SUBROUTINE operR(f1,f2,igrid,nf,nv)

! Apply the R operator:
! given face integrals f1 on primal faces, map to dual cell
! integrals f2, on grid igrid.

! This formulation stores the coefficients in the transpose of
! the original formulation to avoid an MPI reduce. It is
! mathematically equivalent to the original formulation.

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, nv
REAL*8, INTENT(IN) :: f1(nf)
REAL*8, INTENT(OUT) :: f2(nv)
INTEGER :: iv0, if1, ix1, ne1
REAL*8 :: temp

! ----------------------------------------------------------------

! Loop over vertices
DO iv0 = 1, nvert(igrid)
  ne1 = nrxsten(iv0,igrid)
  ! Collect contributions from surrounding faces
  temp = 0.0d0
  DO ix1 = 1, ne1
    if1 = rxsten(iv0,ix1,igrid)
    temp = temp + f1(if1)*rxcoeff(iv0,ix1,igrid)
  ENDDO
  f2(iv0) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE operR

! ================================================================

SUBROUTINE operT(f1,f2,igrid,ne,nf)

! Apply the T operator:
! compute cell integrals of 2 x kinetic energy from edge integrals
! of normal fluxes

USE grid
IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, ne, nf
REAL*8, INTENT(IN) :: f1(ne)
REAL*8, INTENT(OUT) :: f2(nf)
INTEGER :: if1, ix1, ix2, ne1, ie1, ie2
REAL*8 :: temp

! ----------------------------------------------------------------

! Loop over faces
DO if1 = 1, nface(igrid)
  ne1 = ntsten(if1,igrid)
  temp = 0.0d0
  ! Loop over all pairs of edges of this cell
  DO ix1 = 1, ne1
    ie1 = tsten(if1,ix1,igrid)
    DO ix2 = 1, ne1
      ie2 = tsten(if1,ix2,igrid)
      temp = temp + tcoeff(if1,ix1,ix2,igrid)*f1(ie1)*f1(ie2)
    ENDDO
  ENDDO
  f2(if1) = temp
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE operT

! ================================================================

SUBROUTINE restrict(f1,nf1,f2,nf2,igrid)

! To perform the restriction operation needed for a multigrid solver.
! Restrict field f1 from grid igrid + 1 to grid igrid and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms).

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: nf1, nf2, igrid
REAL*8, INTENT(IN) :: f1(nf1)
REAL*8, INTENT(OUT) :: f2(nf2)
INTEGER :: if1, if2, ix
REAL*8 :: wgt

! Safety check
IF (nf2 .ne. nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine restrict'
  STOP
ENDIF

DO if2 = 1, nf2
  f2(if2) = 0.0d0
  DO ix = 1, ninj(if2,igrid)
    if1 = injsten(if2,ix,igrid)
    wgt = injwgt(if2,ix,igrid)
    f2(if2) = f2(if2) + wgt*f1(if1)
  ENDDO
ENDDO

! ----------------------------------------------------------------

END SUBROUTINE restrict

! ================================================================

SUBROUTINE prolong(f2,nf2,f1,nf1,igrid)

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

USE grid

IMPLICIT NONE
INTEGER, INTENT(IN) :: nf1, nf2, igrid
REAL*8, INTENT(IN) :: f2(nf2)
REAL*8, INTENT(OUT) :: f1(nf1)
INTEGER :: if1, if2, ix, igridp
REAL*8 :: wgt, f2if2, temp1(nf1), temp2(nf2)

! Safety check
IF (nf2 .ne. nface(igrid)) THEN
  PRINT *,'Wrong size array in subroutine prolong'
  STOP
ENDIF

igridp = igrid + 1
temp2(1:nf2) = f2(1:nf2)/farea(1:nf2,igrid)
temp1 = 0.0d0
DO if2 = 1, nf2
  f2if2 = temp2(if2)
  DO ix = 1, ninj(if2,igrid)
    if1 = injsten(if2,ix,igrid)
    wgt = injwgt(if2,ix,igrid)
    temp1(if1) = temp1(if1) + wgt*f2if2
  ENDDO
ENDDO
f1(1:nf1) = temp1(1:nf1)*farea(1:nf1,igridp)

! ----------------------------------------------------------------

END SUBROUTINE prolong

! ================================================================

SUBROUTINE laplace(f,hf,igrid,nf,ne)

! To apply the Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

USE laplacecoeff

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
INTEGER :: niter
REAL*8, INTENT(IN) :: f(nf)
REAL*8, INTENT(OUT) :: hf(nf)
REAL*8 :: temp1(nf), temp2(ne), temp3(ne)


CALL massL(f,temp1,igrid,nf)
CALL Ddual1(temp1,temp2,igrid,nf,ne)
niter = -20
CALL massMinv(temp2,temp3,igrid,ne,niter)
CALL Dprimal2(temp3,hf,igrid,ne,nf)


END SUBROUTINE laplace

! ================================================================

SUBROUTINE xlaplace(f,hf,igrid,nf,ne)

! To apply the APPROXIMATE Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

USE laplacecoeff

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
REAL*8, INTENT(IN) :: f(nf)
REAL*8, INTENT(OUT) :: hf(nf)
REAL*8 :: temp1(nf), temp2(ne), temp3(ne)


CALL massL(f,temp1,igrid,nf)
CALL Ddual1(temp1,temp2,igrid,nf,ne)
CALL approxMinv(temp2,temp3,igrid,ne)
CALL Dprimal2(temp3,hf,igrid,ne,nf)


END SUBROUTINE xlaplace

! ================================================================

SUBROUTINE residual(f,rhs,res,igrid,nf,ne)

! Compute the residual res in the approximate Poisson equation on grid igrid
! when f is the input field and rhs is the right hand side. Note that
! f, rhs and res are area integrals (2-forms).

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne
REAL*8, INTENT(IN) :: f(nf), rhs(nf)
REAL*8, INTENT(OUT) :: res(nf)


CALL xlaplace(f,res,igrid,nf,ne)
res = rhs - res
!print *,'     residual: ',res(1:5)

END SUBROUTINE residual

! ================================================================

SUBROUTINE relax(f,rhs,igrid,nf,ne,niter)

! To carry out niter Jacobi relaxation iterations for the multigrid
! solver on grid igrid

USE laplacecoeff

IMPLICIT NONE

INTEGER, INTENT(IN) :: igrid, nf, ne, niter
REAL*8, INTENT(IN) :: rhs(nf)
REAL*8, INTENT(INOUT) :: f(nf)
REAL*8, ALLOCATABLE :: res(:), finc(:)
REAL*8 :: u
INTEGER :: iter

ALLOCATE(res(nf), finc(nf))

u = underrel(igrid)
DO iter = 1, niter
  CALL residual(f,rhs,res,igrid,nf,ne)
  finc = res/lapdiag(1:nf,igrid)
  f = f + u*finc
ENDDO

DEALLOCATE(res, finc)

END SUBROUTINE relax

! ================================================================

SUBROUTINE fullmgsolve(phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 xMinv Ddual1 L phi = RHS
!
! using full multigrid algorithm.
! Coefficients are contained in module laplace.

USE grid

IMPLICIT NONE

! Numbers of iterations on coarsest grid and other grids
INTEGER, PARAMETER :: niterc = 10, niter = 2, npass = 1

INTEGER, INTENT(IN) :: ng
REAL*8, INTENT(IN) :: rr(nfacex)
REAL*8, INTENT(OUT) :: phi(nfacex)

INTEGER :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
REAL*8, ALLOCATABLE :: ff(:,:), rf(:,:)
REAL*8 :: temp1(nfacex)

! ------------------------------------------------------------------------

! Allocate space on all grids
ALLOCATE(ff(nfacex,ngrids),rf(nfacex,ngrids))

! One pass should be enough. Warn user if npass is set to
! some other value for testing
IF (npass > 1) PRINT *,'mgsolve: npass = ',npass

! ------------------------------------------------------------------------

! Initialize solution to zero
phi = 0.0d0

! For diagnostics
!nf1 = nface(ngrids)
!ne1 = nedge(ngrids)
!CALL residual(phi,rr,temp1,ngrids,nf1,ne1)
!print *,'Pass ',0,'  RMS residual = ',SQRT(SUM(temp1*temp1)/nf1)

DO ipass = 1, npass

  ! Initialize rhs as residual using latest estimate
  IF (ipass == 1) THEN
    ! No need to do the calculation
    rf(:,ngrids) = rr(:)
  ELSE
    nf1 = nface(ngrids)
    ne1 = nedge(ngrids)
    CALL residual(phi,rr,rf(1,ngrids),ngrids,nf1,ne1)
  ENDIF

  ! Initialize correction to solution on all grids to zero
  ff = 0.0d0

  ! Restrict right hand side to each grid in the hierarchy
  DO igrid = ngrids-1, ngrids-ng+1, -1
    igridp = igrid + 1
    nf1 = nface(igridp)
    nf2 = nface(igrid)
    CALL restrict(rf(1,igridp),nf1,rf(1,igrid),nf2,igrid)
  ENDDO

  ! Iterate to convergence on coarsest grid
  igrid = ngrids-ng+1
  nf1 = nface(igrid)
  ne1 = nedge(igrid)
  ff(1:nf1,igrid) = 0.0d0
  CALL relax(ff(1,igrid),rf(1,igrid),igrid,nf1,ne1,niterc)

  ! Sequence of growing V-cycles
  DO igridp = ngrids-ng+2, ngrids

    igrid = igridp - 1
    nf1 = nface(igridp)
    ne1 = nedge(igridp)
    nf2 = nface(igrid)
    ne2 = nedge(igrid)

    ! Prolong solution to grid igridp
    ! and execute one V-cycle starting from grid igridp

    ! Prolong
    CALL prolong(ff(1,igrid),nf2,ff(1,igridp),nf1,igrid)

    ! Descending part of V-cycle
    DO jgrid = igrid, ngrids-ng+1, -1
    
      jgridp = jgrid + 1
      nf1 = nface(jgridp)
      ne1 = nedge(jgridp)
      nf2 = nface(jgrid)
      ne2 = nedge(jgrid)

      ! Relax on grid jgridp
      CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

      ! Calculate residual on jgridp
      CALL residual(ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)
  
      ! Restrict residual to jgrid
      CALL restrict(temp1,nf1,rf(1,jgrid),nf2,jgrid)

      ! Set correction first guess to zero on grid jgrid-1
      ff(1:nf2,jgrid) = 0.0d0

    ENDDO

    ! Iterate to convergence on coarsest grid
    jgrid = ngrids-ng+1
    nf1 = nface(jgrid)
    ne1 = nedge(jgrid)
    ff(1:nf1,jgrid) = 0.0d0
    CALL relax(ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

    ! Ascending part of V-cycle
    DO jgrid = ngrids-ng+1, igrid

      jgridp = jgrid + 1
      igrid = igrid - 1
      nf1 = nface(jgridp)
      ne1 = nedge(jgridp)
      nf2 = nface(jgrid)
      ne2 = nedge(jgrid)

      ! Prolong correction to jgridp
      CALL prolong(ff(1,jgrid),nf2,temp1,nf1,jgrid)

      ! Add correction to solution on jgridp
      ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

      ! Relax on grid jgridp
      CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

    ENDDO

  ENDDO

  ! Add correction to phi
  phi = phi + ff(:,ngrids)

  ! For diagnostics
  !nf1 = nface(ngrids)
  !ne1 = nedge(ngrids)
  !CALL residual(phi,rr,temp1,ngrids,nf1,ne1)
  !print *,'      RMS residual in fullmgsolve = ',SQRT(SUM(temp1*temp1)/nf1)

ENDDO


DEALLOCATE(ff,rf)

! ----------------------------------------------------------------

END SUBROUTINE fullmgsolve

! ================================================================

SUBROUTINE mgsolve(phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 xMinv Ddual1 L phi = RHS
!
! using a single V-cycle multigrid algorithm.
! Coefficients are contained in module laplace.

USE grid

IMPLICIT NONE

! Numbers of iterations on coarsest grid and other grids
INTEGER, PARAMETER :: niterc = 10, niter = 2

INTEGER, INTENT(IN) :: ng
REAL*8, INTENT(IN) :: rr(nfacex)
REAL*8, INTENT(OUT) :: phi(nfacex)

INTEGER :: ipass, nf1, nf2, ne1, ne2, igrid, igridp, jgrid, jgridp, iter
REAL*8, ALLOCATABLE :: ff(:,:), rf(:,:)
REAL*8 :: temp1(nfacex)

! ------------------------------------------------------------------------

! Allocate space on all grids
ALLOCATE(ff(nfacex,ngrids),rf(nfacex,ngrids))

! ------------------------------------------------------------------------

! Initialize solution to zero
phi = 0.0d0

! Initialize rhs on finest grid
rf(:,ngrids) = rr(:)

! Initialize correction to solution on all grids to zero
ff = 0.0d0


! Descending part of V-cycle
DO jgrid = ngrids-1, ngrids-ng+1, -1
    
  jgridp = jgrid + 1
  nf1 = nface(jgridp)
  ne1 = nedge(jgridp)
  nf2 = nface(jgrid)
  ne2 = nedge(jgrid)

  ! Relax on grid jgridp
  CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

  ! Calculate residual on jgridp
  CALL residual(ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1,ne1)
  
  ! Restrict residual to jgrid
  CALL restrict(temp1,nf1,rf(1,jgrid),nf2,jgrid)

  ! Set correction first guess to zero on grid jgrid
  ff(1:nf2,jgrid) = 0.0d0

ENDDO

! Iterate to convergence on coarsest grid
jgrid = ngrids-ng+1
nf1 = nface(jgrid)
ne1 = nedge(jgrid)
ff(1:nf1,jgrid) = 0.0d0
CALL relax(ff(1,jgrid),rf(1,jgrid),jgrid,nf1,ne1,niterc)

! Ascending part of V-cycle
DO jgrid = ngrids-ng+1, ngrids-1

  jgridp = jgrid + 1
  nf1 = nface(jgridp)
  ne1 = nedge(jgridp)
  nf2 = nface(jgrid)
  ne2 = nedge(jgrid)

  ! Prolong correction to jgridp
  CALL prolong(ff(1,jgrid),nf2,temp1,nf1,jgrid)

  ! Add correction to solution on jgridp
  ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

  ! Relax on grid jgridp
  CALL relax(ff(1,jgridp),rf(1,jgridp),jgridp,nf1,ne1,niter)

ENDDO


! Add correction to phi
phi = phi + ff(:,ngrids)

!  ! For diagnostics
!  nf1 = nface(ngrids)
!  ne1 = nedge(ngrids)
!  CALL residual(phi,rr,temp1,ngrids,nf1,ne1)
!  print *,'     RMS residual in mgsolve = ',SQRT(SUM(temp1*temp1)/nf1)

DEALLOCATE(ff,rf)

! ----------------------------------------------------------------

END SUBROUTINE mgsolve

! ================================================================

SUBROUTINE readgrid

! To allocate array space for the grid information in module grid
! and to read the information from file

USE runtype
USE constants
USE grid
USE channels

IMPLICIT NONE

INTEGER :: if0, ie0, iv0, igrid, ix, ixx, if1, if2, iv1, iv2, ie1

! ----------------------------------------------------------------

! Open file for reading
OPEN(changrid,FILE=ygridfile,FORM='UNFORMATTED')

! First read ngrids
READ(changrid) ngrids


! Allocate nface, nedge, nvert
ALLOCATE(nface(ngrids), nedge(ngrids), nvert(ngrids))

! Read numbers of faces, edges and vertices on each grid
READ(changrid) nface
READ(changrid) nedge
READ(changrid) nvert

! Find maximum values in order to allocated subsequent arrays
nfacex = MAXVAL(nface)
nedgex = MAXVAL(nedge)
nvertx = MAXVAL(nvert)

! Allocate neoff, neofv
ALLOCATE(neoff(nfacex,ngrids), neofv(nvertx,ngrids))

! Read the numbers of edges of each face and vertex on each grid
neoff = 0
neofv = 0
READ(changrid) ((neoff(if0,igrid),          &
                    if0 = 1, nface(igrid)), &
                    igrid = 1, ngrids)
READ(changrid) ((neofv(iv0,igrid),          &
                    iv0 = 1, nvert(igrid)), &
                    igrid = 1, ngrids)

! Find maximum values in order to allocate subsequent arrays
nefmx = MAXVAL(neoff)
nevmx = MAXVAL(neofv)


! Allocate connectivity arrays arrays
ALLOCATE(fnxtf(nfacex,nefmx,ngrids), eoff(nfacex,nefmx,ngrids), &
         voff(nfacex,nefmx,ngrids),  fnxte(nedgex,2,ngrids),    &
         vofe(nedgex,2,ngrids),      fofv(nvertx,nevmx,ngrids), &
         eofv(nvertx,nevmx,ngrids))

! Read the connectivity arrays
READ(changrid) (((fnxtf(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((eoff(if0,ix,igrid),           &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((voff(if0,ix,igrid),           &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nefmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((fnxte(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
READ(changrid) (((vofe(ie0,ix,igrid),           &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, 2),                &
                     igrid = 1, ngrids)
READ(changrid) (((fofv(iv0,ix,igrid),           &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((eofv(iv0,ix,igrid),           &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nevmx),            &
                     igrid = 1, ngrids)


! Allocate the geometrical information arrays
ALLOCATE(flong(nfacex,ngrids), flat(nfacex,ngrids),  &
         vlong(nvertx,ngrids), vlat(nvertx,ngrids),  &
         farea(nfacex,ngrids), varea(nvertx,ngrids), &
         ldist(nedgex,ngrids), ddist(nedgex,ngrids), &
         fareamin(ngrids))

! Read the geometrical information arrays
READ(changrid) ((flong(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((flat(if0,igrid),               &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((vlong(iv0,igrid),              &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((vlat(iv0,igrid),               &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((farea(if0,igrid),              &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((varea(iv0,igrid),              &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((ldist(ie0,igrid),              &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((ddist(ie0,igrid),              &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)

! Dimensionalize
farea = farea*rearth*rearth
varea = varea*rearth*rearth
ldist = ldist*rearth
ddist = ddist*rearth

! Determine smallest face area on each grid
DO igrid = 1, ngrids
  fareamin(igrid) = MINVAL(farea(1:nface(igrid),igrid))
ENDDO


! Allocate arrays for size of operator stencils
ALLOCATE(nlsten(nfacex,ngrids), nmsten(nedgex,ngrids), &
         njsten(nvertx,ngrids), nhsten(nedgex,ngrids), &
         nrsten(nfacex,ngrids), nrxsten(nvertx,ngrids), &
         nwsten(nedgex,ngrids), ntsten(nfacex,ngrids))

! Read the sizes of the operator stencils on each grid
nlsten = 0
nmsten = 0
njsten = 0
nhsten = 0
nrsten = 0
nrxsten = 0
nwsten = 0
ntsten = 0
READ(changrid) ((nlsten(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((nmsten(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((njsten(iv0,igrid),             &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((nhsten(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((nrsten(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((nrxsten(iv0,igrid),            &
                     iv0 = 1, nvert(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((nwsten(ie0,igrid),             &
                     ie0 = 1, nedge(igrid)),    &
                     igrid = 1, ngrids)
READ(changrid) ((ntsten(if0,igrid),             &
                     if0 = 1, nface(igrid)),    &
                     igrid = 1, ngrids)

! Find maximum values in order to allocate subsequent arrays
nlsmx = MAXVAL(nlsten)
nmsmx = MAXVAL(nmsten)
njsmx = MAXVAL(njsten)
nhsmx = MAXVAL(nhsten)
nrsmx = MAXVAL(nrsten)
nrxsmx = MAXVAL(nrxsten)
nwsmx = MAXVAL(nwsten)
ntsmx = MAXVAL(ntsten)

PRINT *,'Maximum stencil sizes:'
PRINT *,'massL ...  ',nlsmx
PRINT *,'massM ...  ',nmsmx
PRINT *,'HodgeJ ... ',njsmx
PRINT *,'HodgeH ... ',nhsmx
PRINT *,'operR ...  ',nrxsmx
PRINT *,'operW ...  ',nwsmx
PRINT *,'operT ...  ',ntsmx
PRINT *,' '


! Allocate arrays for operator stencils and coefficients
ALLOCATE(lsten(nfacex,nlsmx,ngrids), msten(nedgex,nmsmx,ngrids), &
         jsten(nvertx,njsmx,ngrids), hsten(nedgex,nhsmx,ngrids), &
         rsten(nfacex,nrsmx,ngrids), rxsten(nvertx,nrxsmx,ngrids), &
         wsten(nedgex,nwsmx,ngrids), tsten(nfacex,ntsmx,ngrids))
ALLOCATE(lmass(nfacex,nlsmx,ngrids), mmass(nedgex,nmsmx,ngrids), &
         jstar(nvertx,njsmx,ngrids), hstar(nedgex,nhsmx,ngrids), &
         rcoeff(nfacex,nrsmx,ngrids), rxcoeff(nvertx,nrxsmx,ngrids), &
         wcoeff(nedgex,nwsmx,ngrids), tcoeff(nfacex,ntsmx,ntsmx,ngrids), &
         jlump(nvertx,ngrids), mlump(nedgex,ngrids), hlump(nedgex,ngrids))

! Read the operator stencils and coefficients
READ(changrid) (((lsten(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nlsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((msten(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nmsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((jsten(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((hsten(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((rsten(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((rxsten(iv0,ix,igrid),         &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nrxsmx),           &
                     igrid = 1, ngrids)
READ(changrid) (((wsten(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nwsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((tsten(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ntsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((lmass(if0,ix,igrid),          &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nlsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((mmass(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nmsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((jstar(iv0,ix,igrid),          &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, njsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((hstar(ie0,ix,igrid),          &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nhsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((rcoeff(if0,ix,igrid),         &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, nrsmx),            &
                     igrid = 1, ngrids)
READ(changrid) (((rxcoeff(iv0,ix,igrid),        &
                     iv0 = 1, nvert(igrid)),    &
                     ix = 1, nrxsmx),           &
                     igrid = 1, ngrids)
READ(changrid) (((wcoeff(ie0,ix,igrid),         &
                     ie0 = 1, nedge(igrid)),    &
                     ix = 1, nwsmx),            &
                     igrid = 1, ngrids)
READ(changrid) ((((tcoeff(if0,ix,ixx,igrid),    &
                     if0 = 1, nface(igrid)),    &
                     ix = 1, ntsmx),            &
                     ixx = 1, ntsmx),           &
                     igrid = 1, ngrids)

! Dimensionalize
lmass = lmass/(rearth*rearth)

! Construct the tables eoffin and eofvin
ALLOCATE(eoffin(nfacex,nefmx,ngrids), eofvin(nvertx,nevmx,ngrids))
DO igrid = 1, ngrids

  DO if1 = 1, nface(igrid)
    DO ix = 1, neoff(if1,igrid)
      ie1 = eoff(if1,ix,igrid)
      if2 = fnxte(ie1,1,igrid)
      IF (if1 == if2) THEN
        ! This edge points out of face if1
	eoffin(if1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points into face if1
	eoffin(if1,ix,igrid) = 1.0d0
      ENDIF
    ENDDO
  ENDDO

  DO iv1 = 1, nvert(igrid)
    DO ix = 1, neofv(iv1,igrid)
      ie1 = eofv(iv1,ix,igrid)
      iv2 = vofe(ie1,1,igrid)
      IF (iv1 == iv2) THEN
        ! This edge points away from vertex iv1
	eofvin(iv1,ix,igrid) = -1.0d0
      ELSE
        ! This edge points towards vertex iv1
	eofvin(iv1,ix,igrid) = 1.0d0
      ENDIF
    ENDDO
  ENDDO

ENDDO


! Allocate array for size of restriction stencil
ALLOCATE(ninj(nfacex,ngrids-1))

! Read the size of the restriction stencil on each grid
ninj = 0
READ(changrid) ((ninj(if0,igrid),              &
                    if0 = 1, nface(igrid)),    &
                    igrid = 1, ngrids-1)

! Find maximum value in order to allocate subsequent arrays
ninjmx = MAXVAL(ninj)

! Allocate arrays for restriction stencils and weights
ALLOCATE(injsten(nfacex,ninjmx,ngrids-1))
ALLOCATE(injwgt(nfacex,ninjmx,ngrids-1))

! Read the restriction stencil and weights
READ(changrid) (((injsten(if0,ix,igrid),       &
                    if0 = 1, nface(igrid)),    &
                    ix = 1, ninjmx),           &
                    igrid = 1, ngrids-1)
READ(changrid) (((injwgt(if0,ix,igrid),        &
                    if0 = 1, nface(igrid)),    &
                    ix = 1, ninjmx),           &
                    igrid = 1, ngrids-1)


! -------------------------------------------------------------------

END SUBROUTINE readgrid

!     ===============================================================
!
      SUBROUTINE LL2XYZ(LONG,LAT,X,Y,Z)
!
!     To convert longitude and latitude to cartesian coordinates
!     on the unit sphere
!
      IMPLICIT NONE
!
      REAL*8 LONG,LAT,X,Y,Z,CLN,SLN,CLT,SLT
!
!     ------------------------------------------------------------------
!
      SLN=SIN(LONG)
      CLN=COS(LONG)
      SLT=SIN(LAT)
      CLT=COS(LAT)
!
      X=CLN*CLT
      Y=SLN*CLT
      Z=SLT
!
!     ------------------------------------------------------------------
!
      RETURN
      END
!
!     ==================================================================
!
      SUBROUTINE XYZ2LL(X,Y,Z,LONG,LAT)
!
!     To convert cartesian coordinates to longitude and latitude
!
      IMPLICIT NONE
!   
      REAL*8 X,Y,Z,LONG,LAT,PI,TLN,TLT,R
!
!     -------------------------------------------------------------------
!
      PI=4.0D0*ATAN(1.0D0)
!
      IF (X.EQ.0.0D0) THEN
        IF (Y.GE.0.0D0) THEN
          LONG=0.5D0*PI
        ELSE
          LONG=1.5D0*PI
        ENDIF
      ELSE
        TLN=Y/X
        LONG=ATAN(TLN)
        IF (X.LT.0.0D0) THEN
          LONG=LONG+PI
        ENDIF
        IF (LONG.LT.0.0D0) THEN
          LONG=LONG+2.0D0*PI
        ENDIF
      ENDIF
!
      R=SQRT(X*X+Y*Y)
      IF (R.EQ.0.0D0) THEN
        IF (Z.GT.0.0D0) THEN
          LAT=0.5D0*PI
        ELSE
          LAT=-0.5D0*PI
        ENDIF
      ELSE
        TLT=Z/R
        LAT=ATAN(TLT)
      ENDIF
!
!     --------------------------------------------------------------------
!
      RETURN
      END
!
!     ====================================================================
!
      SUBROUTINE STAREA2(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AREA)
!
!     Calculate the area of the spherical triangle whose corners
!     have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
!     The formula below is more robust to roundoff error than the
!     better known sum of angle - PI formula
!
      IMPLICIT NONE
!
      REAL*8 X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2, &
             D0,D1,D2,S,T0,T1,T2,T3,AREA
!
!
!     Distances between pairs of points
      CALL SPDIST(X0,Y0,Z0,X1,Y1,Z1,D2)
      CALL SPDIST(X1,Y1,Z1,X2,Y2,Z2,D0)
      CALL SPDIST(X2,Y2,Z2,X0,Y0,Z0,D1)
!
!     Half perimeter
      S=0.5D0*(D0+D1+D2)
!
!     Tangents
      T0 = TAN(0.5D0*(S-D0))
      T1 = TAN(0.5D0*(S-D1))
      T2 = TAN(0.5D0*(S-D2))
      T3 = TAN(0.5D0*S)
!
!     Area
      AREA = 4.0D0*ATAN(SQRT(T0*T1*T2*T3))
!
      RETURN
      END
!
!     ===================================================================
!
      SUBROUTINE SPDIST(X1,Y1,Z1,X2,Y2,Z2,S)
!
!     Calculate the spherical distance S between two points with Cartesian
!     coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

      IMPLICIT NONE

      REAL*8 X1, Y1, Z1, X2, Y2, Z2, S, DX, DY, DZ, AD


      DX = X2 - X1
      DY = Y2 - Y1
      DZ = Z2 - Z1
      AD = SQRT(DX*DX + DY*DY + DZ*DZ)
      S = 2.0D0*ASIN(0.5D0*AD)


      RETURN
      END
!
!     ===================================================================

SUBROUTINE centroid(if0,long,lat,igrid)

! Find the centroid of cell if0 on grid igrid

USE grid
IMPLICIT NONE
INTEGER, INTENT(IN) :: if0, igrid
REAL*8, INTENT(OUT) :: long, lat
INTEGER :: ixe, ie1, iv1, iv2
REAL*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          xc, yc, zc, a, aby3, mag

! -----------------------------------------------------------------------

! Coordinates of `centre' of face (i.e. dual vertex)
long1 = flong(if0,igrid)
lat1 = flat(if0,igrid)
CALL ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the face
! Hence find area of face and centroid
xc = 0.0d0
yc = 0.0d0
zc = 0.0d0
DO ixe = 1, neoff(if0,igrid)
  ie1 = eoff(if0,ixe,igrid)
  iv1 = vofe(ie1,1,igrid)
  iv2 = vofe(ie1,2,igrid)
  long1 = vlong(iv1,igrid)
  lat1 = vlat(iv1,igrid)
  CALL ll2xyz(long1,lat1,x1,y1,z1)
  long1 = vlong(iv2,igrid)
  lat1 = vlat(iv2,igrid)
  CALL ll2xyz(long1,lat1,x2,y2,z2)
  CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
  aby3 = a/3.0d0
  xc = xc + (x0 + x1 + x2)*aby3
  yc = yc + (y0 + y1 + y2)*aby3
  zc = zc + (z0 + z1 + z2)*aby3
ENDDO
mag = SQRT(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag
CALL xyz2ll(xc,yc,zc,long,lat)

! -----------------------------------------------------------------------

END SUBROUTINE centroid

!     ===================================================================

SUBROUTINE dual_centroid(iv0,long,lat,igrid)

! Find the centroid of dual cell iv0 on grid igrid

USE grid
IMPLICIT NONE
INTEGER, INTENT(IN) :: iv0, igrid
REAL*8, INTENT(OUT) :: long, lat
INTEGER :: ixe, ie1, iv1, iv2
REAL*8 :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
          xc, yc, zc, a, aby3, mag

! -----------------------------------------------------------------------

! Coordinates of `centre' of dual cell (i.e. vertex)
long1 = vlong(iv0,igrid)
lat1 = vlat(iv0,igrid)
CALL ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the dual cell
! Hence find area of dual cell and centroid
xc = 0.0d0
yc = 0.0d0
zc = 0.0d0
DO ixe = 1, neofv(iv0,igrid)
  ie1 = eofv(iv0,ixe,igrid)
  iv1 = fnxte(ie1,1,igrid)
  iv2 = fnxte(ie1,2,igrid)
  long1 = flong(iv1,igrid)
  lat1 = flat(iv1,igrid)
  CALL ll2xyz(long1,lat1,x1,y1,z1)
  long1 = flong(iv2,igrid)
  lat1 = flat(iv2,igrid)
  CALL ll2xyz(long1,lat1,x2,y2,z2)
  CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
  aby3 = a/3.0d0
  xc = xc + (x0 + x1 + x2)*aby3
  yc = yc + (y0 + y1 + y2)*aby3
  zc = zc + (z0 + z1 + z2)*aby3
ENDDO
mag = SQRT(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag
CALL xyz2ll(xc,yc,zc,long,lat)

! -----------------------------------------------------------------------

END SUBROUTINE dual_centroid

!     ===================================================================

SUBROUTINE testpoisson

! To test the solution of thePoisson problem

USE grid
USE laplacecoeff
IMPLICIT NONE

! Number of passes
INTEGER :: npass = 10

INTEGER :: nf, ne, nv, if1, iv1, ipass, nprt
REAL*8, ALLOCATABLE :: psi0(:), zeta(:), psi(:), ff1(:), ff2(:), ff3(:), ff4(:), temp1(:), temp2(:)
REAL*8 :: long, lat, psibar

! -------------------------------------------------------------------

nf = nface(ngrids)
ne = nedge(ngrids)
nv = nvert(ngrids)

print *,' '
print *,'--------------------------'
print *,' '
print *,'Testing mgsolve '
print *,' '

! Number of values to print for testing
nprt = 5

ALLOCATE(psi0(nf),zeta(nf),psi(nf),ff1(nf),ff2(nf),ff3(nf),ff4(nf),temp1(ne),temp2(ne))


! Build coefficients used in Laplacian operator on all grids
CALL buildlap


! Set up test data
! Large-scale part
DO if1 = 1, nf
  long = flong(if1,ngrids)
  lat = flat(if1,ngrids)
  ! psi0(if1) = SIN(lat)
  psi0(if1) = COS(lat)*SIN(long)
ENDDO
! Plus small-scale part
psi0(10) = 10.0d0*psi0(10)
! Remove global mean (to ensure unambiguous result)
psibar = SUM(psi0*farea(:,ngrids))/SUM(farea(:,ngrids))
psi0 = psi0 - psibar
! Convert to area integrals
ff1 = psi0*farea(:,ngrids)
print *,'Original field ff1 =     ',ff1(1:nprt)
print *,' '

! Calculate laplacian
CALL laplace(ff1,zeta,ngrids,nf,ne)


! Now invert laplacian to check we get back to where we started

! Initialize result to zero
! Note psi will be stream function times grid cell area
psi = 0.0d0
temp2 = 0.0d0

! Iterate several passes
DO ipass = 1, npass

  print *,'Pass ',ipass

  ! Compute residual based on latest estimate
  CALL massL(psi,ff2,ngrids,nf)

  CALL Ddual1(ff2,temp1,ngrids,nf,ne)

  ! Improve the estimate temp2 that we obtained in the previous pass
  CALL massMinv(temp1,temp2,ngrids,ne,4)

  CALL Dprimal2(temp2,ff3,ngrids,ne,nf)

  ! Residual
  ff4 = zeta - ff3

  ! Now solve the approximate Poisson equation
  ! D2 xMinv D1bar L psi' = residual
  CALL mgsolve(ff3,ff4,ngrids)

  ! And increment best estimate
  ! *** We could think about adding beta*ff3 here and tuning beta for optimal convergence ***
  psi = psi + ff3

  ! Remove global mean (to ensure unambiguous result)
  psibar = SUM(psi)/SUM(farea(:,ngrids))
  psi = psi - psibar*farea(:,ngrids)

  print *,'Original field ff1      = ',ff1(1:nprt)
  print *,'Soln of Poisson eqn psi = ', psi(1:nprt)
  ff4 = (ff1-psi)/farea(:,ngrids)
  print *,'RMS err in global problem = ',sqrt(sum(ff4*ff4)/nf)
  print *,' '

ENDDO




END SUBROUTINE testpoisson

! ====================================================================
