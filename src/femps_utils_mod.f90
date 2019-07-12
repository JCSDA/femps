module femps_utils_mod

implicit none
private

public ll2xyz, xyz2ll, starea2, spdist

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine ll2xyz(long,lat,x,y,z)

! To convert longitude and latitude to cartesian coordinates
! on the unit sphere

implicit none

real*8 long,lat,x,y,z,cln,sln,clt,slt

sln=sin(long)
cln=cos(long)
slt=sin(lat)
clt=cos(lat)

x=cln*clt
y=sln*clt
z=slt

end subroutine ll2xyz

! --------------------------------------------------------------------------------------------------

subroutine xyz2ll(x,y,z,long,lat)

! To convert cartesian coordinates to longitude and latitude

implicit none

real*8 x,y,z,long,lat,pi,tln,tlt,r

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

end subroutine xyz2ll

! --------------------------------------------------------------------------------------------------

subroutine starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area)

! Calculate the area of the spherical triangle whose corners
! have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
! The formula below is more robust to roundoff error than the
! better known sum of angle - PI formula

implicit none

real*8 x0,y0,z0,x1,y1,z1,x2,y2,z2, &
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

real*8 x1, y1, z1, x2, y2, z2, s, dx, dy, dz, ad

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1
ad = sqrt(dx*dx + dy*dy + dz*dz)
s = 2.0d0*asin(0.5d0*ad)

end subroutine spdist

! --------------------------------------------------------------------------------------------------

end module femps_utils_mod
