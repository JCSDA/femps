module femps_testgrid_mod

use femps_kinds_mod
use femps_grid_mod
use femps_utils_mod

implicit none
private
public cstestgrid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine cstestgrid(grid,flavour_in,nsmooth_in)

implicit none
type(fempsgrid),   intent(inout) :: grid
integer, optional, intent(in)    :: flavour_in
integer, optional, intent(in)    :: nsmooth_in

! flavour = cube sphere grid flavour
! 1 - equiangular cube
! 2 - barycentric as in TCD14
! 3 - centroidal

! nsmooth = number of smoothing iterations.
! For flavour 1, no iterations are taken. nsmooth is ignored.
! For flavour 2, nsmooth = 1 is recommended. It must be at least 1 for a consistent FV H operator.
! For flavour 3, nsmooth = 3 is recommended.

integer :: igrid, i, j, ixv, p1, p2, pp, jr, iv, iv0, n, n2, &
           ie0, ie1, ie2, if0, if1, if2, if3, iv1, iv2, ix1, ix2, &
           ixmin, ifmin, if21, if22, iv11, iv12, iv21, iv22, &
           ismooth, flavour, nsmooth

real(kind=kind_real) :: pi, dlambda, lambda1, lambda2, t1, t2, long, lat, &
                        x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                        x5, y5, z5, x6, y6, z6, piby4, xc, yc, zc, x0, y0, z0, &
                        rmag, atri, aface, lmn, lmx, dmn, dmx, dav, cs, s, &
                        theta, thetamin, sn, d1x, d1y, d1z, d2x, d2y, d2z, &
                        lambdaf1, lambdaf2, tf1, tf2, xf1, yf1, zf1, xf2, yf2, zf2, &
                        xf3, yf3, zf3, xf4, yf4, zf4, xf5, yf5, zf5, xf6, yf6, zf6, &
                        amin, amax

logical :: lfound

! Default setup
! -------------
flavour = 1
if (present(flavour_in)) flavour = flavour_in

nsmooth = 3
if (present(nsmooth_in)) nsmooth = nsmooth_in

! Constants
! ---------
piby4 = ATAN(1.0d0)
pi = 4.0d0*piby4


! Generate the grid
! -----------------
DO igrid = 1, grid%ngrids

! Size of panels on this grid
n = grid%n0*(2**(igrid-1))
n2 = n*n
dlambda = 0.5d0*pi/n

grid%nface(igrid) = 6*n2
grid%nedge(igrid) = 2*grid%nface(igrid)
grid%nvert(igrid) = grid%nface(igrid) + 2

!
! Loop over vertices/faces of one panel
DO j = 1, n
  lambda2 = (j-1)*dlambda - piby4
  lambdaf2 = (j-0.5d0)*dlambda - piby4
  t2 = TAN(lambda2)
  tf2 = TAN(lambdaf2)
  DO i = 1, n
    lambda1 = (i-1)*dlambda - piby4
    lambdaf1 = (i-0.5d0)*dlambda - piby4
    t1 = TAN(lambda1)
    tf1 = TAN(lambdaf1)

!   Set up coordinates of vertices and faces
!   Panel 1
!   Index of vertex
    ixv = (j-1)*n + i
!   Cartesian coordinates of vertex
    x1 = 1.0d0/SQRT(1.0d0 + t1*t1 + t2*t2)
    y1 = x1*t1
    z1 = x1*t2
!   Lat long coordinates of vertex
    CALL xyz2ll(x1,y1,z1,long,lat)
    grid%vlong(ixv,igrid) = long
    grid%vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf1 = 1.0d0/SQRT(1.0d0 + tf1*tf1 + tf2*tf2)
    yf1 = xf1*tf1
    zf1 = xf1*tf2
!   Lat long coordinates of face
    CALL xyz2ll(xf1,yf1,zf1,long,lat)
    grid%flong(ixv,igrid) = long
    grid%flat(ixv,igrid) = lat

!   Panel 2
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x2 = -y1
    y2 = x1
    z2 = z1
!   Lat long coordinates of vertex
    CALL xyz2ll(x2,y2,z2,long,lat)
    grid%vlong(ixv,igrid) = long
    grid%vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf2 = -yf1
    yf2 = xf1
    zf2 = zf1
!   Lat long coordinates of face
    CALL xyz2ll(xf2,yf2,zf2,long,lat)
    grid%flong(ixv,igrid) = long
    grid%flat(ixv,igrid) = lat

!   Panel 3
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x3 = x2
    y3 = -z2
    z3 = y2
!   Lat long coordinates of vertex
    CALL xyz2ll(x3,y3,z3,long,lat)
    grid%vlong(ixv,igrid) = long
    grid%vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf3 = xf2
    yf3 = -zf2
    zf3 = yf2
!   Lat long coordinates of face
    CALL xyz2ll(xf3,yf3,zf3,long,lat)
    grid%flong(ixv,igrid) = long
    grid%flat(ixv,igrid) = lat

!   Panel 4
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x4 = -z3
    y4 = y3
    z4 = x3
!   Lat long coordinates of vertex
    CALL xyz2ll(x4,y4,z4,long,lat)
    grid%vlong(ixv,igrid) = long
    grid%vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf4 = -zf3
    yf4 = yf3
    zf4 = xf3
!   Lat long coordinates of face
    CALL xyz2ll(xf4,yf4,zf4,long,lat)
    grid%flong(ixv,igrid) = long
    grid%flat(ixv,igrid) = lat

!   Panel 5
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x5 = -y4
    y5 = x4
    z5 = z4
!   Lat long coordinates of vertex
    CALL xyz2ll(x5,y5,z5,long,lat)
    grid%vlong(ixv,igrid) = long
    grid%vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf5 = -yf4
    yf5 = xf4
    zf5 = zf4
!   Lat long coordinates of face
    CALL xyz2ll(xf5,yf5,zf5,long,lat)
    grid%flong(ixv,igrid) = long
    grid%flat(ixv,igrid) = lat

!   Panel 6
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x6 = x5
    y6 = -z5
    z6 = y5
!   Lat long coordinates of vertex
    CALL xyz2ll(x6,y6,z6,long,lat)
    grid%vlong(ixv,igrid) = long
    grid%vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf6 = xf5
    yf6 = -zf5
    zf6 = yf5
!   Lat long coordinates of face
    CALL xyz2ll(xf6,yf6,zf6,long,lat)
    grid%flong(ixv,igrid) = long
    grid%flat(ixv,igrid) = lat

!   Set up incidence tables ignoring complications at
!   panel edges
    DO p1 = 1, 6
      ixv = (p1 - 1)*n2 + (j-1)*n + i

!     Edges of the face
      grid%eoff(ixv,1,igrid) = 2*ixv - 1
      grid%eoff(ixv,2,igrid) = 2*ixv
      grid%eoff(ixv,3,igrid) = 2*ixv + 1
      grid%eoff(ixv,4,igrid) = 2*ixv + 2*n
!     Vertices of the face
      grid%voff(ixv,1,igrid) = ixv
      grid%voff(ixv,2,igrid) = ixv + 1
      grid%voff(ixv,3,igrid) = ixv + n + 1
      grid%voff(ixv,4,igrid) = ixv + n
!     Faces neighboring this face
      grid%fnxtf(ixv,1,igrid) = ixv - 1
      grid%fnxtf(ixv,2,igrid) = ixv - n
      grid%fnxtf(ixv,3,igrid) = ixv + 1
      grid%fnxtf(ixv,4,igrid) = ixv + n
!     Edges incident on the vertex
      grid%eofv(ixv,1,igrid) = 2*ixv - 2
      grid%eofv(ixv,2,igrid) = 2*(ixv-n) - 1
      grid%eofv(ixv,3,igrid) = 2*ixv
      grid%eofv(ixv,4,igrid) = 2*ixv - 1

    ENDDO

  ENDDO
ENDDO


! Now sort out complications at panel edges
DO j = 1, n
  jr = n + 1 - j

  DO pp = 1, 3

    ! Odd numbered panels
    p1 = 2*pp - 1

    ! Left edge of panel p1 joins to top edge of panel p1 - 2
    ! Reverse order
    p2 = MODULO(p1 + 3, 6) + 1
    ixv = (p1 - 1)*n2 + n*(j - 1) + 1
    grid%fnxtf(ixv,1,igrid) = p2*n2 - n + jr
    grid%eofv(ixv,1,igrid) = 2*p2*n2 - 2*n - 1 + 2*(jr + 1)

    ! Bottom edge of panel p1 joins to top edge of panel p1 - 1
    p2 = MODULO(p1 + 4, 6) + 1
    ixv = (p1 - 1)*n2 + j
    grid%fnxtf(ixv,2,igrid) = p2*n2 - n + j
    grid%eofv(ixv,2,igrid) = 2*p2*n2 - 2*n - 1 + 2*j

    ! Right edge of panel p1 joins to left edge of panel p1 + 1
    p2 = MODULO(p1, 6) + 1
    ixv = (p1 - 1)*n2 + n*j
    grid%eoff(ixv,3,igrid) = 2*(p2 - 1)*n2 + 2*(j - 1)*n + 1
    grid%voff(ixv,2,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1
    grid%voff(ixv,3,igrid) = (p2 - 1)*n2 + j*n + 1
    grid%fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1

    ! Top edge of panel p1 joins to left edge of panel p1 + 2
    ! Reverse order
    p2 = MODULO(p1 + 1, 6) + 1
    ixv = p1*n2 - n + j
    grid%eoff(ixv,4,igrid) = 2*(p2 - 1)*n2 + 2*(jr - 1)*n + 1
    grid%voff(ixv,3,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1
    grid%voff(ixv,4,igrid) = (p2 - 1)*n2 + jr*n + 1
    grid%fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1

    ! Even numbered panels
    p1 = 2*pp

    ! Left edge of panel p1 joins to right edge of panel p1 - 1
    p2 = MODULO(p1 + 4, 6) + 1
    ixv = (p1 - 1)*n2 + n*(j - 1) + 1
    grid%fnxtf(ixv,1,igrid) = (p2 - 1)*n2 + n*j
    grid%eofv(ixv,1,igrid) = 2*(p2 - 1)*n2 + 2*n*j

    ! Bottom edge of panel p1 joins to right edge of panel p1 - 2
    ! Reverse order
    p2 = MODULO(p1 + 3, 6) + 1
    ixv = (p1 - 1)*n2 + j
    grid%fnxtf(ixv,2,igrid) = (p2 - 1)*n2 + n*jr
    grid%eofv(ixv,2,igrid) = 2*(p2 - 1)*n2 + 2*n*(jr + 1)

    ! Top edge of panel p1 joins to bottom edge of panel p1 + 1
    p2 = MODULO(p1, 6) + 1
    ixv = p1*n2 - n + j
    grid%eoff(ixv,4,igrid) = 2*(p2 - 1)*n2 + 2*j
    grid%voff(ixv,3,igrid) = (p2 - 1)*n2 + j + 1
    grid%voff(ixv,4,igrid) = (p2 - 1)*n2 + j
    grid%fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + j

    ! Right edge of panel p1 joins to bottom edge of panel p1 + 2
    ! Reverse order
    p2 = MODULO(p1 + 1, 6) + 1
    ixv = (p1-1)*n2 + n*j
    grid%eoff(ixv,3,igrid) = 2*(p2 - 1)*n2 + 2*jr
    grid%voff(ixv,2,igrid) = (p2 - 1)*n2 + jr + 1
    grid%voff(ixv,3,igrid) = (p2 - 1)*n2 + jr
    grid%fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + jr

  ENDDO

ENDDO

! All faces have 4 edges and vertices
grid%neoff(1:grid%nface(igrid),igrid) = 4

! Almost all vertices have 4 edges (exceptions dealt with below)
grid%neofv(1:grid%nvert(igrid),igrid) = 4

! Vertices not correctly captured by the above
! Corner vertices only have 3 edges
DO pp = 1, 3

  ! Bottom left of odd numbered panels
  p1 = 2*pp - 1
  ixv = (p1 - 1)*n2 + 1
  ! First edge needs to be deleted
  grid%eofv(ixv,1,igrid) = grid%eofv(ixv,2,igrid)
  grid%eofv(ixv,2,igrid) = grid%eofv(ixv,3,igrid)
  grid%eofv(ixv,3,igrid) = grid%eofv(ixv,4,igrid)
  grid%eofv(ixv,4,igrid) = 0
  grid%neofv(ixv,igrid) = 3

  ! Bottom left of even numbered panels
  p1 = 2*pp
  ixv = (p1 - 1)*n2 + 1
  ! Second edge needs to be deleted
  grid%eofv(ixv,2,igrid) = grid%eofv(ixv,3,igrid)
  grid%eofv(ixv,3,igrid) = grid%eofv(ixv,4,igrid)
  grid%eofv(ixv,4,igrid) = 0
  grid%neofv(ixv,igrid) = 3

ENDDO

! Vertex 6*n2 + 1 is at top left of panels 1, 3, and 5
iv = 6*n2 + 1
lambda2 = piby4
t2 = TAN(lambda2)
lambda1 = - piby4
t1 = TAN(lambda1)
! Cartesian coordinates of vertex
x1 = 1.0d0/SQRT(1.0d0 + t1*t1 + t2*t2)
y1 = x1*t1
z1 = x1*t2
! Lat long coordinates of vertex
CALL xyz2ll(x1,y1,z1,long,lat)
grid%vlong(iv,igrid) = long
grid%vlat(iv,igrid) = lat
DO pp = 1, 3
  p1 = 2*pp - 1
  ixv = p1*n2 - n + 1
  grid%voff(ixv,4,igrid) = iv
  grid%eofv(iv,pp,igrid) = 2*p1*n2 - 2*n + 1
ENDDO
grid%neofv(iv,igrid) = 3

! Vertex 6*n2 + 2 is at bottom right of panels 2, 4, and 6
iv = 6*n2 + 2
x1 = -x1
y1 = -y1
z1 = -z1
! Lat long coordinates of vertex
CALL xyz2ll(x1,y1,z1,long,lat)
grid%vlong(iv,igrid) = long
grid%vlat(iv,igrid) = lat
DO pp = 1, 3
  p1 = 2*pp
  ixv = (p1 - 1)*n2 + n
  grid%voff(ixv,2,igrid) = iv
  grid%eofv(iv,pp,igrid) = 2*(p1 - 1)*n2 + 2*n
ENDDO
grid%neofv(iv,igrid) = 3


ENDDO ! End of main loop over grids


! Now construct inverse tables
! First initialize entries to zero
grid%fnxte = 0
grid%vofe = 0
grid%fofv = 0

DO igrid = 1, grid%ngrids

  DO j = 1, grid%nface(igrid)
    DO i = 1, 4
      ie1 = grid%eoff(j,i,igrid)
      CALL addtab(grid%nedgex,2,grid%fnxte(1,1,igrid),ie1,j)
      ixv = grid%voff(j,i,igrid)
      CALL addtab(grid%nvertx,4,grid%fofv(1,1,igrid),ixv,j)
    ENDDO
  ENDDO

  DO j = 1, grid%nvert(igrid)
    DO i = 1, 4
      ixv = grid%eofv(j,i,igrid)
      IF (ixv > 0) THEN
        CALL addtab(grid%nedgex,2,grid%vofe(1,1,igrid),ixv,j)
      ENDIF
    ENDDO
  ENDDO

ENDDO

! Calculate geometrical quantities
DO igrid = 1, grid%ngrids

  ! Smoothing iterations
  IF (flavour == 1) THEN

    ! No smoothing for equiangular cube

  ELSEIF (flavour == 2) THEN

    ! TCD barycentric cube
    DO ismooth = 1, nsmooth

      ! First locate face centres at barycentres of
      ! surrounding vertices
      DO if1 = 1, grid%nface(igrid)
        xc = 0.0d0
        yc = 0.0d0
        zc = 0.0d0
        DO i = 1, 4
          ixv = grid%voff(if1,i,igrid)
          long = grid%vlong(ixv,igrid)
          lat = grid%vlat(ixv,igrid)
          CALL ll2xyz(long,lat,x1,y1,z1)
          xc = xc + x1
          yc = yc + y1
          zc = zc + z1
        ENDDO
        rmag = 1.0d0/SQRT(xc*xc + yc*yc + zc*zc)
        xc = xc*rmag
        yc = yc*rmag
        zc = zc*rmag
        CALL xyz2ll(xc,yc,zc,long,lat)
        grid%flong(if1,igrid) = long
        grid%flat(if1,igrid) = lat
      ENDDO

      ! Next relocate vertices at barycentres of
      ! surrounding face centres - needed for H operator
      DO iv1 = 1, grid%nvert(igrid)
        xc = 0.0d0
        yc = 0.0d0
        zc = 0.0d0
        DO i = 1, grid%neofv(iv1,igrid)
          if1 = grid%fofv(iv1,i,igrid)
          long = grid%flong(if1,igrid)
          lat = grid%flat(if1,igrid)
          CALL ll2xyz(long,lat,x1,y1,z1)
          xc = xc + x1
          yc = yc + y1
          zc = zc + z1
        ENDDO
        rmag = 1.0d0/SQRT(xc*xc + yc*yc + zc*zc)
        xc = xc*rmag
        yc = yc*rmag
        zc = zc*rmag
        CALL xyz2ll(xc,yc,zc,long,lat)
        grid%vlong(iv1,igrid) = long
        grid%vlat(iv1,igrid) = lat
      ENDDO

    ENDDO

  ELSEIF (flavour == 3) THEN

    ! Centroidal cube
    DO ismooth = 1, nsmooth
      ! Move faces to centroids
      DO if0 = 1, grid%nface(igrid)
        CALL centroid(grid,if0,long,lat,igrid)
        grid%flong(if0,igrid) = long
        grid%flat(if0,igrid) = lat
      ENDDO
      ! Move vertices to centroids
      DO iv0 = 1, grid%nvert(igrid)
        CALL dual_centroid(grid,iv0,long,lat,igrid)
        grid%vlong(iv0,igrid) = long
        grid%vlat(iv0,igrid) = lat
      ENDDO
    ENDDO

  ELSE

    PRINT *,'flavour = ',flavour,'  not recognized. Please pick another.'
    STOP

  ENDIF

! Tabulate areas
  DO if1 = 1, grid%nface(igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
!   Compute face area
    aface = 0.0d0
    DO i = 1, 4
      ie1 = grid%eoff(if1,i,igrid)
      ixv = grid%vofe(ie1,1,igrid)
      long = grid%vlong(ixv,igrid)
      lat = grid%vlat(ixv,igrid)
      CALL ll2xyz(long,lat,x1,y1,z1)
      ixv = grid%vofe(ie1,2,igrid)
      long = grid%vlong(ixv,igrid)
      lat = grid%vlat(ixv,igrid)
      CALL ll2xyz(long,lat,x2,y2,z2)
      CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,atri)
      aface = aface + atri
    ENDDO
    grid%farea(if1,igrid) = aface
  ENDDO
  amin = MINVAL(grid%farea(1:grid%nface(igrid),igrid))
  amax = MAXVAL(grid%farea(1:grid%nface(igrid),igrid))
  PRINT *,'Grid ',igrid,' min and max face area ',amin,amax,' ratio ',amin/amax

! Tabulate lengths of edges and distances between face centres
! across each edge
  lmn=5.0D0
  lmx=0.0D0
  dmn=5.0D0
  dmx=0.0D0
  dav=0.0D0
  DO ie0 = 1, grid%nedge(igrid)
!   Vertices at ends of this edge
    iv1 = grid%vofe(ie0,1,igrid)
    iv2 = grid%vofe(ie0,2,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    long = grid%vlong(iv2,igrid)
    lat = grid%vlat(iv2,igrid)
    CALL ll2xyz(long,lat,x2,y2,z2)
    CALL spdist(x1,y1,z1,x2,y2,z2,s)
    grid%ldist(ie0,igrid) = s
    lmn = MIN(lmn,grid%ldist(ie0,igrid))
    lmx = MIN(lmx,grid%ldist(ie0,igrid))
!   Faces either side of this edge
    if1 = grid%fnxte(ie0,1,igrid)
    if2 = grid%fnxte(ie0,2,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    long = grid%flong(if2,igrid)
    lat = grid%flat(if2,igrid)
    CALL ll2xyz(long,lat,x2,y2,z2)
    CALL spdist(x1,y1,z1,x2,y2,z2,s)
    grid%ddist(ie0,igrid) = s
    dmn = MIN(dmn,grid%ddist(ie0,igrid))
    dmx = MIN(dmx,grid%ddist(ie0,igrid))
    dav = dav + grid%ddist(ie0,igrid)/grid%nedge(igrid)
  ENDDO

ENDDO


! Sort FNXTF into anticlockwise order on each grid
! and sort EOFF to correspond to FNXTF
! Also sort fofv into anticlockwise order
DO igrid = 1, grid%ngrids

  DO if0 = 1, grid%nface(igrid)

!   Coordinates of face if0
    long = grid%flong(if0,igrid)
    lat = grid%flat(if0,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    DO ix1 = 1, 2
!     Coordinates of IX1'th neighbour
      if1 = grid%fnxtf(if0,ix1,igrid)
      long = grid%flong(if1,igrid)
      lat = grid%flat(if1,igrid)
      CALL ll2xyz(long,lat,x1,y1,z1)
      d1x = x1 - x0
      d1y = y1 - y0
      d1z = z1 - z0
!     Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0
      DO ix2 = ix1 + 1, 4
!       Coordinates of IX2'th neighbour
        if2 = grid%fnxtf(if0,ix2,igrid)
        long = grid%flong(if2,igrid)
        lat = grid%flat(if2,igrid)
        CALL ll2xyz(long,lat,x2,y2,z2)
        d2x=x2 - x0
        d2y=y2 - y0
        d2z=z2 - z0
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0*(d1y*d2z - d1z*d2y) &
           + y0*(d1z*d2x - d1x*d2z) &
           + z0*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        IF ((theta < thetamin) .AND. (theta > 0.0d0)) THEN
          ixmin = ix2
          ifmin = if2
          thetamin = theta
        ENDIF
      ENDDO
!     The face in position IXMIN belongs in position IX1+1 so swap them
      if3 = grid%fnxtf(if0,ix1+1,igrid)
      grid%fnxtf(if0,ix1+1,igrid) = ifmin
      grid%fnxtf(if0,ixmin,igrid) = if3
    ENDDO

    DO ix1 = 1, 4
      if1 = grid%fnxtf(if0,ix1,igrid)
      ix2 = ix1 - 1
      lfound = .FALSE.
      DO WHILE (.NOT. lfound)
         ix2 = ix2 + 1
         ie1 = grid%eoff(if0,ix2,igrid)
         if21 = grid%fnxte(ie1,1,igrid)
         if22 = grid%fnxte(ie1,2,igrid)
         IF ((if21 + if22) == (if0 + if1)) lfound = .TRUE.
      ENDDO
!     Edge IE2 corresponds to face IF1
      grid%eoff(if0,ix2,igrid) = grid%eoff(if0,ix1,igrid)
      grid%eoff(if0,ix1,igrid) = ie1
    ENDDO

  ENDDO

  DO iv0 = 1, grid%nvert(igrid)
!   Coordinates of vertex iv0
    long = grid%vlong(iv0,igrid)
    lat = grid%vlat(iv0,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    DO ix1 = 1, grid%neofv(iv0,igrid) - 2
!     Coordinates of IX1'th face
      if1 = grid%fofv(iv0,ix1,igrid)
      long = grid%flong(if1,igrid)
      lat = grid%flat(if1,igrid)
      CALL ll2xyz(long,lat,x1,y1,z1)
      d1x = x1 - x0
      d1y = y1 - y0
      d1z = z1 - z0
!     Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0
      DO ix2 = ix1 + 1, grid%neofv(iv0,igrid)
!       Coordinates of IX2'th neighbour
        if2 = grid%fofv(iv0,ix2,igrid)
        long = grid%flong(if2,igrid)
        lat = grid%flat(if2,igrid)
        CALL ll2xyz(long,lat,x2,y2,z2)
        d2x=x2 - x0
        d2y=y2 - y0
        d2z=z2 - z0
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0*(d1y*d2z - d1z*d2y) &
           + y0*(d1z*d2x - d1x*d2z) &
           + z0*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        IF ((theta < thetamin) .AND. (theta > 0.0d0)) THEN
          ixmin = ix2
          ifmin = if2
          thetamin = theta
        ENDIF
      ENDDO
!     The face in position IXMIN belongs in position IX1+1 so swap them
      if3 = grid%fofv(iv0,ix1+1,igrid)
      grid%fofv(iv0,ix1+1,igrid) = ifmin
      grid%fofv(iv0,ixmin,igrid) = if3
    ENDDO
  ENDDO

ENDDO

!
! Order VOFF so that the k'th vertex is between the
! k'th and (k+1)'th edges in EOFF
DO igrid = 1, grid%ngrids
  DO if0 = 1, grid%nface(igrid)
    DO ix1 = 1, grid%neoff(if0,igrid)
      ix2 = ix1 + 1
      IF (ix2 > grid%neoff(if0,igrid)) ix2 = 1
      ie1 = grid%eoff(if0,ix1,igrid)
      ie2 = grid%eoff(if0,ix2,igrid)
      ! Find the common vertex of IE1 and IE2
      iv11 = grid%vofe(ie1,1,igrid)
      iv12 = grid%vofe(ie1,2,igrid)
      iv21 = grid%vofe(ie2,1,igrid)
      iv22 = grid%vofe(ie2,2,igrid)
      IF ((iv11 == iv21) .OR. (iv11 == iv22)) THEN
        iv0 = iv11
      ELSEIF ((iv12 == iv21) .OR. (iv12 == iv22)) THEN
        iv0 = iv12
      ELSE
        PRINT *,'Common vertex not found'
        STOP
      ENDIF
      grid%voff(if0,ix1,igrid) = iv0
    ENDDO
  ENDDO
ENDDO


! Sort VOFE so that VOFE(1) -> VOFE(2) (tangent vector)
! is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2) (normal vector)
DO igrid = 1, grid%ngrids
  DO ie0 = 1, grid%nedge(igrid)
    if1 = grid%fnxte(ie0,1,igrid)
    if2 = grid%fnxte(ie0,2,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    long = grid%flong(if2,igrid)
    lat = grid%flat(if2,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    d1x = x1 - x0
    d1y = y1 - y0
    d1z = z1 - z0
    iv1 = grid%vofe(ie0,1,igrid)
    iv2 = grid%vofe(ie0,2,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    long = grid%vlong(iv2,igrid)
    lat = grid%vlat(iv2,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    d2x = x1 - x0
    d2y = y1 - y0
    d2z = z1 - z0
    sn = x0*(d1y*d2z - d1z*d2y) &
       + y0*(d1z*d2x - d1x*d2z) &
       + z0*(d1x*d2y - d1y*d2x)
    IF (sn < 0.0d0) THEN
      ! Swap the two vertices
      grid%vofe(ie0,1,igrid) = iv2
      grid%vofe(ie0,2,igrid) = iv1
    ENDIF
  ENDDO
ENDDO

end subroutine cstestgrid

! --------------------------------------------------------------------------------------------------

end module femps_testgrid_mod
