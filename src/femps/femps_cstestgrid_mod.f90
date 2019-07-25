module femps_cstestgrid_mod

use femps_types_mod, only: fempsgrid

implicit none
private
public cstestgrid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine cstestgrid(grid)

implicit none
type(fempsgrid), intent(inout) :: grid

integer :: igrid, i, j, ixv, p1, p2, pp, jr, iv, iv0, n, n2, &
           ie0, ie1, ie2, if0, if1, if2, if3, iv1, iv2, ix1, ix2, &
           ixmin, ifmin, if21, if22, iv11, iv12, iv21, iv22, &
           ismooth

real*8 :: pi, dlambda, lambda1, lambda2, t1, t2, long, lat, &
          x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
          x5, y5, z5, x6, y6, z6, piby4, xc, yc, zc, x0, y0, z0, &
          rmag, atri, aface, lmn, lmx, dmn, dmx, dav, cs, s, &
          theta, thetamin, sn, d1x, d1y, d1z, d2x, d2y, d2z, &
          lambdaf1, lambdaf2, tf1, tf2, xf1, yf1, zf1, xf2, yf2, zf2, &
          xf3, yf3, zf3, xf4, yf4, zf4, xf5, yf5, zf5, xf6, yf6, zf6, &
          amin, amax

logical :: lfound

!     ---------------------------------------------------------------

! Constants

piby4 = ATAN(1.0d0)
pi = 4.0d0*piby4

!     ---------------------------------------------------------------

DO igrid = 1, grid%ngrids

! Size of panels on this grid
n = n0*(2**(igrid-1))
n2 = n*n
dlambda = 0.5d0*pi/n

nface(igrid) = 6*n2
nedge(igrid) = 2*nface(igrid)
nvert(igrid) = nface(igrid) + 2

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
    vlong(ixv,igrid) = long
    vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf1 = 1.0d0/SQRT(1.0d0 + tf1*tf1 + tf2*tf2)
    yf1 = xf1*tf1
    zf1 = xf1*tf2
!   Lat long coordinates of face
    CALL xyz2ll(xf1,yf1,zf1,long,lat)
    flong(ixv,igrid) = long
    flat(ixv,igrid) = lat

!   Panel 2
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x2 = -y1
    y2 = x1
    z2 = z1
!   Lat long coordinates of vertex
    CALL xyz2ll(x2,y2,z2,long,lat)
    vlong(ixv,igrid) = long
    vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf2 = -yf1
    yf2 = xf1
    zf2 = zf1
!   Lat long coordinates of face
    CALL xyz2ll(xf2,yf2,zf2,long,lat)
    flong(ixv,igrid) = long
    flat(ixv,igrid) = lat

!   Panel 3
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x3 = x2
    y3 = -z2
    z3 = y2
!   Lat long coordinates of vertex
    CALL xyz2ll(x3,y3,z3,long,lat)
    vlong(ixv,igrid) = long
    vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf3 = xf2
    yf3 = -zf2
    zf3 = yf2
!   Lat long coordinates of face
    CALL xyz2ll(xf3,yf3,zf3,long,lat)
    flong(ixv,igrid) = long
    flat(ixv,igrid) = lat

!   Panel 4
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x4 = -z3
    y4 = y3
    z4 = x3
!   Lat long coordinates of vertex
    CALL xyz2ll(x4,y4,z4,long,lat)
    vlong(ixv,igrid) = long
    vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf4 = -zf3
    yf4 = yf3
    zf4 = xf3
!   Lat long coordinates of face
    CALL xyz2ll(xf4,yf4,zf4,long,lat)
    flong(ixv,igrid) = long
    flat(ixv,igrid) = lat

!   Panel 5
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x5 = -y4
    y5 = x4
    z5 = z4
!   Lat long coordinates of vertex
    CALL xyz2ll(x5,y5,z5,long,lat)
    vlong(ixv,igrid) = long
    vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf5 = -yf4
    yf5 = xf4
    zf5 = zf4
!   Lat long coordinates of face
    CALL xyz2ll(xf5,yf5,zf5,long,lat)
    flong(ixv,igrid) = long
    flat(ixv,igrid) = lat

!   Panel 6
!   Index of vertex
    ixv = ixv + n2
!   Cartesian coordinates of vertex
    x6 = x5
    y6 = -z5
    z6 = y5
!   Lat long coordinates of vertex
    CALL xyz2ll(x6,y6,z6,long,lat)
    vlong(ixv,igrid) = long
    vlat(ixv,igrid) = lat
!   Cartesian coordinates of face
    xf6 = xf5
    yf6 = -zf5
    zf6 = yf5
!   Lat long coordinates of face
    CALL xyz2ll(xf6,yf6,zf6,long,lat)
    flong(ixv,igrid) = long
    flat(ixv,igrid) = lat

!   Set up incidence tables ignoring complications at
!   panel edges
    DO p1 = 1, 6
      ixv = (p1 - 1)*n2 + (j-1)*n + i

!     Edges of the face
      eoff(ixv,1,igrid) = 2*ixv - 1
      eoff(ixv,2,igrid) = 2*ixv
      eoff(ixv,3,igrid) = 2*ixv + 1
      eoff(ixv,4,igrid) = 2*ixv + 2*n
!     Vertices of the face
      voff(ixv,1,igrid) = ixv
      voff(ixv,2,igrid) = ixv + 1
      voff(ixv,3,igrid) = ixv + n + 1
      voff(ixv,4,igrid) = ixv + n
!     Faces neighboring this face
      fnxtf(ixv,1,igrid) = ixv - 1
      fnxtf(ixv,2,igrid) = ixv - n
      fnxtf(ixv,3,igrid) = ixv + 1
      fnxtf(ixv,4,igrid) = ixv + n
!     Edges incident on the vertex
      eofv(ixv,1,igrid) = 2*ixv - 2
      eofv(ixv,2,igrid) = 2*(ixv-n) - 1
      eofv(ixv,3,igrid) = 2*ixv
      eofv(ixv,4,igrid) = 2*ixv - 1

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
    fnxtf(ixv,1,igrid) = p2*n2 - n + jr
    eofv(ixv,1,igrid) = 2*p2*n2 - 2*n - 1 + 2*(jr + 1)

    ! Bottom edge of panel p1 joins to top edge of panel p1 - 1
    p2 = MODULO(p1 + 4, 6) + 1
    ixv = (p1 - 1)*n2 + j
    fnxtf(ixv,2,igrid) = p2*n2 - n + j
    eofv(ixv,2,igrid) = 2*p2*n2 - 2*n - 1 + 2*j

    ! Right edge of panel p1 joins to left edge of panel p1 + 1
    p2 = MODULO(p1, 6) + 1
    ixv = (p1 - 1)*n2 + n*j
    eoff(ixv,3,igrid) = 2*(p2 - 1)*n2 + 2*(j - 1)*n + 1
    voff(ixv,2,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1
    voff(ixv,3,igrid) = (p2 - 1)*n2 + j*n + 1
    fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1

    ! Top edge of panel p1 joins to left edge of panel p1 + 2
    ! Reverse order
    p2 = MODULO(p1 + 1, 6) + 1
    ixv = p1*n2 - n + j
    eoff(ixv,4,igrid) = 2*(p2 - 1)*n2 + 2*(jr - 1)*n + 1
    voff(ixv,3,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1
    voff(ixv,4,igrid) = (p2 - 1)*n2 + jr*n + 1
    fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1

    ! Even numbered panels
    p1 = 2*pp

    ! Left edge of panel p1 joins to right edge of panel p1 - 1
    p2 = MODULO(p1 + 4, 6) + 1
    ixv = (p1 - 1)*n2 + n*(j - 1) + 1
    fnxtf(ixv,1,igrid) = (p2 - 1)*n2 + n*j
    eofv(ixv,1,igrid) = 2*(p2 - 1)*n2 + 2*n*j

    ! Bottom edge of panel p1 joins to right edge of panel p1 - 2
    ! Reverse order
    p2 = MODULO(p1 + 3, 6) + 1
    ixv = (p1 - 1)*n2 + j
    fnxtf(ixv,2,igrid) = (p2 - 1)*n2 + n*jr
    eofv(ixv,2,igrid) = 2*(p2 - 1)*n2 + 2*n*(jr + 1)

    ! Top edge of panel p1 joins to bottom edge of panel p1 + 1
    p2 = MODULO(p1, 6) + 1
    ixv = p1*n2 - n + j
    eoff(ixv,4,igrid) = 2*(p2 - 1)*n2 + 2*j
    voff(ixv,3,igrid) = (p2 - 1)*n2 + j + 1
    voff(ixv,4,igrid) = (p2 - 1)*n2 + j
    fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + j

    ! Right edge of panel p1 joins to bottom edge of panel p1 + 2
    ! Reverse order
    p2 = MODULO(p1 + 1, 6) + 1
    ixv = (p1-1)*n2 + n*j
    eoff(ixv,3,igrid) = 2*(p2 - 1)*n2 + 2*jr
    voff(ixv,2,igrid) = (p2 - 1)*n2 + jr + 1
    voff(ixv,3,igrid) = (p2 - 1)*n2 + jr
    fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + jr

  ENDDO

ENDDO

! All faces have 4 edges and vertices
neoff(1:nface(igrid),igrid) = 4

! Almost all vertices have 4 edges (exceptions dealt with below)
neofv(1:nvert(igrid),igrid) = 4

! Vertices not correctly captured by the above
! Corner vertices only have 3 edges
DO pp = 1, 3

  ! Bottom left of odd numbered panels
  p1 = 2*pp - 1
  ixv = (p1 - 1)*n2 + 1
  ! First edge needs to be deleted
  eofv(ixv,1,igrid) = eofv(ixv,2,igrid)
  eofv(ixv,2,igrid) = eofv(ixv,3,igrid)
  eofv(ixv,3,igrid) = eofv(ixv,4,igrid)
  eofv(ixv,4,igrid) = 0
  neofv(ixv,igrid) = 3

  ! Bottom left of even numbered panels
  p1 = 2*pp
  ixv = (p1 - 1)*n2 + 1
  ! Second edge needs to be deleted
  eofv(ixv,2,igrid) = eofv(ixv,3,igrid)
  eofv(ixv,3,igrid) = eofv(ixv,4,igrid)
  eofv(ixv,4,igrid) = 0
  neofv(ixv,igrid) = 3

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
vlong(iv,igrid) = long
vlat(iv,igrid) = lat
DO pp = 1, 3
  p1 = 2*pp - 1
  ixv = p1*n2 - n + 1
  voff(ixv,4,igrid) = iv
  eofv(iv,pp,igrid) = 2*p1*n2 - 2*n + 1
ENDDO
neofv(iv,igrid) = 3

! Vertex 6*n2 + 2 is at bottom right of panels 2, 4, and 6
iv = 6*n2 + 2
x1 = -x1
y1 = -y1
z1 = -z1
! Lat long coordinates of vertex
CALL xyz2ll(x1,y1,z1,long,lat)
vlong(iv,igrid) = long
vlat(iv,igrid) = lat
DO pp = 1, 3
  p1 = 2*pp
  ixv = (p1 - 1)*n2 + n
  voff(ixv,2,igrid) = iv
  eofv(iv,pp,igrid) = 2*(p1 - 1)*n2 + 2*n
ENDDO
neofv(iv,igrid) = 3


ENDDO ! End of main loop over grids


! ---------------------------------------------------------------

! Now construct inverse tables
! First initialize entries to zero
fnxte = 0
vofe = 0
fofv = 0

DO igrid = 1, ngrids

  DO j = 1, nface(igrid)
    DO i = 1, 4
      ie1 = eoff(j,i,igrid)
      CALL addtab(fnxte(1,1,igrid),ie1,j,nedgex,2)
      ixv = voff(j,i,igrid)
      CALL addtab(fofv(1,1,igrid),ixv,j,nvertx,4)
    ENDDO
  ENDDO

  DO j = 1, nvert(igrid)
    DO i = 1, 4
      ixv = eofv(j,i,igrid)
      IF (ixv > 0) THEN
        CALL addtab(vofe(1,1,igrid),ixv,j,nedgex,2)
      ENDIF
    ENDDO
  ENDDO

ENDDO


! Calculate geometrical quantities
DO igrid = 1, ngrids

  ! Smoothing iterations
  IF (flavour == 1) THEN

    ! No smoothing for equiangular cube

  ELSEIF (flavour == 2) THEN

    ! TCD barycentric cube
    DO ismooth = 1, nsmooth

      ! First locate face centres at barycentres of
      ! surrounding vertices
      DO if1 = 1, nface(igrid)
        xc = 0.0d0
        yc = 0.0d0
        zc = 0.0d0
        DO i = 1, 4
          ixv = voff(if1,i,igrid)
          long = vlong(ixv,igrid)
          lat = vlat(ixv,igrid)
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
        flong(if1,igrid) = long
        flat(if1,igrid) = lat
      ENDDO

      ! Next relocate vertices at barycentres of
      ! surrounding face centres - needed for H operator
      DO iv1 = 1, nvert(igrid)
        xc = 0.0d0
        yc = 0.0d0
        zc = 0.0d0
        DO i = 1, neofv(iv1,igrid)
          if1 = fofv(iv1,i,igrid)
          long = flong(if1,igrid)
          lat = flat(if1,igrid)
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
        vlong(iv1,igrid) = long
        vlat(iv1,igrid) = lat
      ENDDO

    ENDDO

  ELSEIF (flavour == 3) THEN

    ! Centroidal cube
    DO ismooth = 1, nsmooth
      ! Move faces to centroids
      DO if0 = 1, nface(igrid)
        CALL centroid(if0,long,lat,igrid)
        flong(if0,igrid) = long
        flat(if0,igrid) = lat
      ENDDO
      ! Move vertices to centroids
      DO iv0 = 1, nvert(igrid)
        CALL dual_centroid(iv0,long,lat,igrid)
        vlong(iv0,igrid) = long
        vlat(iv0,igrid) = lat
      ENDDO
    ENDDO

  ELSE

    PRINT *,'flavour = ',flavour,'  not recognized. Please pick another.'
    STOP

  ENDIF

! Tabulate areas
  DO if1 = 1, nface(igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
!   Compute face area
    aface = 0.0d0
    DO i = 1, 4
      ie1 = eoff(if1,i,igrid)
      ixv = vofe(ie1,1,igrid)
      long = vlong(ixv,igrid)
      lat = vlat(ixv,igrid)
      CALL ll2xyz(long,lat,x1,y1,z1)
      ixv = vofe(ie1,2,igrid)
      long = vlong(ixv,igrid)
      lat = vlat(ixv,igrid)
      CALL ll2xyz(long,lat,x2,y2,z2)
      CALL starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,atri)
      aface = aface + atri
    ENDDO
    farea(if1,igrid) = aface
  ENDDO
  amin = MINVAL(farea(1:nface(igrid),igrid))
  amax = MAXVAL(farea(1:nface(igrid),igrid))
  PRINT *,'Grid ',igrid,' min and max face area ',amin,amax,' ratio ',amin/amax

! Tabulate lengths of edges and distances between face centres
! across each edge
  lmn=5.0D0
  lmx=0.0D0
  dmn=5.0D0
  dmx=0.0D0
  dav=0.0D0
  DO ie0 = 1, nedge(igrid)
!   Vertices at ends of this edge
    iv1 = vofe(ie0,1,igrid)
    iv2 = vofe(ie0,2,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    long = vlong(iv2,igrid)
    lat = vlat(iv2,igrid)
    CALL ll2xyz(long,lat,x2,y2,z2)
    CALL spdist(x1,y1,z1,x2,y2,z2,s)
    ldist(ie0,igrid) = s
    lmn = MIN(lmn,ldist(ie0,igrid))
    lmx = MIN(lmx,ldist(ie0,igrid))
!   Faces either side of this edge
    if1 = fnxte(ie0,1,igrid)
    if2 = fnxte(ie0,2,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    long = flong(if2,igrid)
    lat = flat(if2,igrid)
    CALL ll2xyz(long,lat,x2,y2,z2)
    CALL spdist(x1,y1,z1,x2,y2,z2,s)
    ddist(ie0,igrid) = s
    dmn = MIN(dmn,ddist(ie0,igrid))
    dmx = MIN(dmx,ddist(ie0,igrid))
    dav = dav + ddist(ie0,igrid)/nedge(igrid)
  ENDDO

ENDDO


! Sort FNXTF into anticlockwise order on each grid
! and sort EOFF to correspond to FNXTF
! Also sort fofv into anticlockwise order
DO igrid = 1, ngrids

  DO if0 = 1, nface(igrid)

!   Coordinates of face if0
    long = flong(if0,igrid)
    lat = flat(if0,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    DO ix1 = 1, 2
!     Coordinates of IX1'th neighbour
      if1 = fnxtf(if0,ix1,igrid)
      long = flong(if1,igrid)
      lat = flat(if1,igrid)
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
        if2 = fnxtf(if0,ix2,igrid)
        long = flong(if2,igrid)
        lat = flat(if2,igrid)
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
      if3 = fnxtf(if0,ix1+1,igrid)
      fnxtf(if0,ix1+1,igrid) = ifmin
      fnxtf(if0,ixmin,igrid) = if3
    ENDDO

    DO ix1 = 1, 4
      if1 = fnxtf(if0,ix1,igrid)
      ix2 = ix1 - 1
      lfound = .FALSE.
      DO WHILE (.NOT. lfound)
         ix2 = ix2 + 1
         ie1 = eoff(if0,ix2,igrid)
         if21 = fnxte(ie1,1,igrid)
         if22 = fnxte(ie1,2,igrid)
         IF ((if21 + if22) == (if0 + if1)) lfound = .TRUE.
      ENDDO
!     Edge IE2 corresponds to face IF1
      eoff(if0,ix2,igrid) = eoff(if0,ix1,igrid)
      eoff(if0,ix1,igrid) = ie1
    ENDDO

  ENDDO

  DO iv0 = 1, nvert(igrid)
!   Coordinates of vertex iv0
    long = vlong(iv0,igrid)
    lat = vlat(iv0,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    DO ix1 = 1, neofv(iv0,igrid) - 2
!     Coordinates of IX1'th face
      if1 = fofv(iv0,ix1,igrid)
      long = flong(if1,igrid)
      lat = flat(if1,igrid)
      CALL ll2xyz(long,lat,x1,y1,z1)
      d1x = x1 - x0
      d1y = y1 - y0
      d1z = z1 - z0
!     Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0
      DO ix2 = ix1 + 1, neofv(iv0,igrid)
!       Coordinates of IX2'th neighbour
        if2 = fofv(iv0,ix2,igrid)
        long = flong(if2,igrid)
        lat = flat(if2,igrid)
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
      if3 = fofv(iv0,ix1+1,igrid)
      fofv(iv0,ix1+1,igrid) = ifmin
      fofv(iv0,ixmin,igrid) = if3
    ENDDO
  ENDDO

ENDDO

!
! Order VOFF so that the k'th vertex is between the
! k'th and (k+1)'th edges in EOFF
DO igrid = 1, ngrids
  DO if0 = 1, nface(igrid)
    DO ix1 = 1, neoff(if0,igrid)
      ix2 = ix1 + 1
      IF (ix2 > neoff(if0,igrid)) ix2 = 1
      ie1 = eoff(if0,ix1,igrid)
      ie2 = eoff(if0,ix2,igrid)
      ! Find the common vertex of IE1 and IE2
      iv11 = vofe(ie1,1,igrid)
      iv12 = vofe(ie1,2,igrid)
      iv21 = vofe(ie2,1,igrid)
      iv22 = vofe(ie2,2,igrid)
      IF ((iv11 == iv21) .OR. (iv11 == iv22)) THEN
        iv0 = iv11
      ELSEIF ((iv12 == iv21) .OR. (iv12 == iv22)) THEN
        iv0 = iv12
      ELSE
        PRINT *,'Common vertex not found'
        STOP
      ENDIF
      voff(if0,ix1,igrid) = iv0
    ENDDO
  ENDDO
ENDDO


! Sort VOFE so that VOFE(1) -> VOFE(2) (tangent vector)
! is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2) (normal vector)
DO igrid = 1, ngrids
  DO ie0 = 1, nedge(igrid)
    if1 = fnxte(ie0,1,igrid)
    if2 = fnxte(ie0,2,igrid)
    long = flong(if1,igrid)
    lat = flat(if1,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    long = flong(if2,igrid)
    lat = flat(if2,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    d1x = x1 - x0
    d1y = y1 - y0
    d1z = z1 - z0
    iv1 = vofe(ie0,1,igrid)
    iv2 = vofe(ie0,2,igrid)
    long = vlong(iv1,igrid)
    lat = vlat(iv1,igrid)
    CALL ll2xyz(long,lat,x0,y0,z0)
    long = vlong(iv2,igrid)
    lat = vlat(iv2,igrid)
    CALL ll2xyz(long,lat,x1,y1,z1)
    d2x = x1 - x0
    d2y = y1 - y0
    d2z = z1 - z0
    sn = x0*(d1y*d2z - d1z*d2y) &
       + y0*(d1z*d2x - d1x*d2z) &
       + z0*(d1x*d2y - d1y*d2x)
    IF (sn < 0.0d0) THEN
      ! Swap the two vertices
      vofe(ie0,1,igrid) = iv2
      vofe(ie0,2,igrid) = iv1
    ENDIF
  ENDDO
ENDDO

end subroutine cstestgrid

! --------------------------------------------------------------------------------------------------

end module femps_cstestgrid_mod
