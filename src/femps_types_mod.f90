module femps_types_mod

implicit none
private

public grid, operators

! Type to hold all the grid information
! -------------------------------------
type fempsgrid

! 1 - equiangular cube
! 2 - barycentric as in TCD14
! 3 - centroidal
integer, parameter :: cs_flavour

! Number of grids in multigrid hierarchy
integer, parameter :: ngrids = 6

! n x n cells on each panel
! Smallest and largest n
integer, parameter :: n0 = 3, nx = n0*(2**(ngrids-1)), nx2 = nx*nx

! Number of smoothing iterations.
! For flavour 1, no iterations are taken. nsmooth is ignored.
! For flavour 2, nsmooth = 1 is recommended. It must be at least 1 for a consistent FV H operator.
! For flavour 3, nsmooth = 3 is recommended.
integer, parameter :: nsmooth = 3

integer, parameter :: nfacex = 6*nx*nx, &
                      nedgex = 2*nfacex, &
                      nvertx = nfacex + 2

integer :: neoff(nfacex,ngrids), neofv(nvertx,ngrids), &
           nface(ngrids), nedge(ngrids), nvert(ngrids)
integer :: fnxtf(nfacex,4,ngrids), eoff(nfacex,4,ngrids), &
           voff(nfacex,4,ngrids), fnxte(nedgex,2,ngrids), &
           vofe(nedgex,2,ngrids), fofv(nvertx,4,ngrids), &
           eofv(nvertx,4,ngrids)
real*8 :: flong(nfacex,ngrids), flat(nfacex,ngrids), &
          vlong(nvertx,ngrids), vlat(nvertx,ngrids), &
          farea(nfacex,ngrids), &
          ldist(nedgex,ngrids), ddist(nedgex,ngrids)


end type fempsgrid

! Type to hold the FEM operators
! ------------------------------
type fempsops



end type fempsops

end module femps_types_mod
