module femps_utils_mod

use femps_kinds_mod
use femps_grid_mod

implicit none

public

interface ll2xyz
   module procedure ll2xyz_sca
   module procedure ll2xyz_vec
end interface

interface xyz2ll
   module procedure xyz2ll_sca
   module procedure xyz2ll_vec
end interface

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine ll2xyz_sca(long,lat,x,y,z)

! To convert longitude and latitude to cartesian coordinates
! on the unit sphere

implicit none
real(kind=kind_real), intent(in)  :: long,lat
real(kind=kind_real), intent(out) :: x, y, z

real(kind=kind_real) :: cln,sln,clt,slt

sln=sin(long)
cln=cos(long)
slt=sin(lat)
clt=cos(lat)

x=cln*clt
y=sln*clt
z=slt

end subroutine ll2xyz_sca

! --------------------------------------------------------------------------------------------------

subroutine ll2xyz_vec(long,lat,x)

implicit none
real(kind=kind_real), intent(in)  :: long,lat
real(kind=kind_real), intent(out) :: x(3)

call ll2xyz_sca(long,lat,x(1),x(2),x(3))

end subroutine ll2xyz_vec

! --------------------------------------------------------------------------------------------------

subroutine xyz2ll_sca(x,y,z,long,lat)

implicit none
real(kind=kind_real), intent(in)  :: x,y,z
real(kind=kind_real), intent(out) :: long,lat

real(kind=kind_real) :: pi,tln,tlt,r

! To convert cartesian coordinates to longitude and latitude
! ----------------------------------------------------------

pi=4.0d0*atan(1.0d0)

if (x.eq.0.0d0) then
  if (y.ge.0.0d0) then
    long=0.5d0*pi
  else
    long=1.5d0*pi
  endif
else
  tln=y/x
  long=atan(tln)
  if (x.lt.0.0d0) then
    long=long+pi
  endif
  if (long.lt.0.0d0) then
    long=long+2.0d0*pi
  endif
endif

r=sqrt(x*x+y*y)
if (r.eq.0.0d0) then
  if (z.gt.0.0d0) then
    lat=0.5d0*pi
  else
    lat=-0.5d0*pi
  endif
else
  tlt=z/r
  lat=atan(tlt)
endif

end subroutine xyz2ll_sca

! --------------------------------------------------------------------------------------------------

subroutine xyz2ll_vec(x,long,lat)

implicit none
real(kind=kind_real), intent(in)  :: x(3)
real(kind=kind_real), intent(out) :: long,lat

call xyz2ll_sca(x(1),x(2),x(3),long,lat)

end subroutine xyz2ll_vec

! --------------------------------------------------------------------------------------------------

subroutine starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area)

! Calculate the area of the spherical triangle whose corners
! have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
! The formula below is more robust to roundoff error than the
! better known sum of angle - PI formula

implicit none

real(kind=kind_real) x0,y0,z0,x1,y1,z1,x2,y2,z2, &
       d0,d1,d2,s,t0,t1,t2,t3,area


! Distances between pairs of points
call spdist(x0,y0,z0,x1,y1,z1,d2)
call spdist(x1,y1,z1,x2,y2,z2,d0)
call spdist(x2,y2,z2,x0,y0,z0,d1)

! Half perimeter
s=0.5d0*(d0+d1+d2)

! Tangents
t0 = tan(0.5d0*(s-d0))
t1 = tan(0.5d0*(s-d1))
t2 = tan(0.5d0*(s-d2))
t3 = tan(0.5d0*s)

! Area
area = 4.0d0*atan(sqrt(t0*t1*t2*t3))

end subroutine starea2

! --------------------------------------------------------------------------------------------------

subroutine spdist(x1,y1,z1,x2,y2,z2,s)

! Calculate the spherical distance S between two points with Cartesian
! coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

implicit none

real(kind=kind_real) x1, y1, z1, x2, y2, z2, s, dx, dy, dz, ad

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1
ad = sqrt(dx*dx + dy*dy + dz*dz)
s = 2.0d0*asin(0.5d0*ad)

end subroutine spdist

! --------------------------------------------------------------------------------------------------

subroutine addtab(dim1,dim2,tab,index,entry)

implicit none
integer, intent(in)    :: dim1,dim2
integer, intent(inout) :: tab(dim1,dim2)
integer, intent(in)    :: index
integer, intent(in)    :: entry

integer :: i

! Add an entry to a table
! -----------------------

i=0
100 continue
i=i+1
if (i.gt.dim2) then
  print *,'**********'
  print *,'table full'
  print *,'**********'
  print *,index,entry,dim1,dim2
  print *,tab(index,:)
  stop
endif
if (tab(index,i).ne.0) goto 100
tab(index,i)=entry

end subroutine addtab

! --------------------------------------------------------------------------------------------------

subroutine centroid(grid,if0,long,lat,igrid)

implicit none

type(fempsgrid),      intent(in)  :: grid
integer,              intent(in)  :: if0
integer,              intent(in)  :: igrid
real(kind=kind_real), intent(out) :: long
real(kind=kind_real), intent(out) :: lat

integer :: ixe, ie1, iv1, iv2
real(kind=kind_real) :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
                        xc, yc, zc, a, aby3, mag

! Find the centroid of cell if0 on grid igrid
! -------------------------------------------

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
mag = sqrt(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag
call xyz2ll(xc,yc,zc,long,lat)

end subroutine centroid

! --------------------------------------------------------------------------------------------------

subroutine dual_centroid(grid,iv0,long,lat,igrid)

implicit none
type(fempsgrid),      intent(in)  :: grid
integer,              intent(in)  :: iv0
integer,              intent(in)  :: igrid
real(kind=kind_real), intent(out) :: long
real(kind=kind_real), intent(out) :: lat

integer :: ixe, ie1, iv1, iv2
real(kind=kind_real) :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
                        xc, yc, zc, a, aby3, mag

! Find the centroid of dual cell iv0 on grid igrid
! ------------------------------------------------

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

mag = sqrt(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag

call xyz2ll(xc,yc,zc,long,lat)

end subroutine dual_centroid

! --------------------------------------------------------------------------------------------------

subroutine findinlist(i,list,nl,ix)

implicit none
integer, intent(in)    :: i
integer, intent(in)    :: nl
integer, intent(inout) :: list(nl)
integer, intent(out)   :: ix

! Find the index ix to the entry i in the input list.
! if i is not already in the list then it is added into
! the first empty (i.e. zero) position.

ix = 1
do while ((list(ix) .ne. i) .and. (list(ix) .ne. 0))
  if (ix == nl) then
    print *,'search past end of list in subroutine findinlist.'
    print *,'i = ',i
    print *,'list = ',list
    stop
  endif
  ix = ix + 1
enddo
list(ix) = i

end subroutine findinlist

! --------------------------------------------------------------------------------------------------

subroutine triangle(x1,x2,x3,l1sq,l2sq,l3sq,area)

implicit none
real(kind=kind_real), intent(in)  :: x1(3), x2(3), x3(3)
real(kind=kind_real), intent(out) :: l1sq, l2sq, l3sq, area

real(kind=kind_real) :: a, b, c

! Compute the squares of the lengths of the sides and the area of
! a triangle defined by three points x1, x2 and x3 in 3d euclidean
! space. Side 1 is opposite vertex 1, etc.

! squares of side by pythagoras
l1sq = (x2(1) - x3(1))**2 + (x2(2) - x3(2))**2 + (x2(3) - x3(3))**2
l2sq = (x3(1) - x1(1))**2 + (x3(2) - x1(2))**2 + (x3(3) - x1(3))**2
l3sq = (x1(1) - x2(1))**2 + (x1(2) - x2(2))**2 + (x1(3) - x2(3))**2

! area by heron's formula
a = sqrt(l1sq)
b = sqrt(l2sq)
c = sqrt(l3sq)
area = 0.25d0*sqrt((a + b + c)*(b + c - a)*(c + a - b)*(a + b - c))

end subroutine triangle

! --------------------------------------------------------------------------------------------------

end module femps_utils_mod
