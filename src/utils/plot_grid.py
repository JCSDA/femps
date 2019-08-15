import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import math

# User input
# ----------
gridfile = '/discover/nobackup/drholdaw/JediDev/femps/build-int-release/test/grid.nc4'

igrid = 3 # Grid number to plot


fv3gridfile = '/gpfsm/dnb31/drholdaw/JediDev/fv3-bundle/build-intel-17.0.7.259-release-default/fv3-jedi/test/fv3grid_c0012.nc4'

# File handle
# -----------
fh_grid = Dataset(gridfile)
fh_fv3grid = Dataset(fv3gridfile)


# Dimensions
# ----------
ngrids   = len(fh_grid.dimensions['ngrids'])
nfacex   = len(fh_grid.dimensions['nfacex'])
nvertx   = len(fh_grid.dimensions['nvertx'])
nedgex   = len(fh_grid.dimensions['nedgex'])
dimfnxtf = len(fh_grid.dimensions['dimfnxtf'])
dimeoff  = len(fh_grid.dimensions['dimeoff'])
dimvoff  = len(fh_grid.dimensions['dimvoff'])
dimfnxte = len(fh_grid.dimensions['dimfnxte'])
dimvofe  = len(fh_grid.dimensions['dimvofe'])
dimfofv  = len(fh_grid.dimensions['dimfofv'])
dimeofv  = len(fh_grid.dimensions['dimeofv'])

nface = np.zeros([ngrids], dtype=int)
nedge = np.zeros([ngrids], dtype=int)
nvert = np.zeros([ngrids], dtype=int)
neoff = np.zeros([ngrids, nfacex], dtype=int)
neofv = np.zeros([ngrids, nvertx], dtype=int)
fnxtf = np.zeros([ngrids, dimfnxtf, nfacex], dtype=int)
eoff  = np.zeros([ngrids, dimeoff, nfacex], dtype=int)
voff  = np.zeros([ngrids, dimvoff, nfacex], dtype=int)
fnxte = np.zeros([ngrids, dimfnxte, nedgex], dtype=int)
vofe  = np.zeros([ngrids, dimvofe, nedgex], dtype=int)
fofv  = np.zeros([ngrids, dimfofv, nvertx], dtype=int)
eofv  = np.zeros([ngrids, dimeofv, nvertx], dtype=int)
flong = np.zeros([ngrids, nfacex])
flat  = np.zeros([ngrids, nfacex])
vlong = np.zeros([ngrids, nvertx])
vlat  = np.zeros([ngrids, nvertx])
farea = np.zeros([ngrids, nfacex])
ldist = np.zeros([ngrids, nedgex])
ddist = np.zeros([ngrids, nedgex])

nface[:]     = fh_grid.variables['nface'][:]
nedge[:]     = fh_grid.variables['nedge'][:]
nvert[:]     = fh_grid.variables['nvert'][:]
neoff[:,:]   = fh_grid.variables['neoff'][:,:]
neofv[:,:]   = fh_grid.variables['neofv'][:,:]
fnxtf[:,:,:] = fh_grid.variables['fnxtf'][:,:,:]
eoff [:,:,:] = fh_grid.variables['eoff' ][:,:,:]
voff [:,:,:] = fh_grid.variables['voff' ][:,:,:]
fnxte[:,:,:] = fh_grid.variables['fnxte'][:,:,:]
vofe [:,:,:] = fh_grid.variables['vofe' ][:,:,:]
fofv [:,:,:] = fh_grid.variables['fofv' ][:,:,:]
eofv [:,:,:] = fh_grid.variables['eofv' ][:,:,:]
flong[:,:]   = fh_grid.variables['flong'][:,:]
flat [:,:]   = fh_grid.variables['flat' ][:,:]
vlong[:,:]   = fh_grid.variables['vlong'][:,:]
vlat [:,:]   = fh_grid.variables['vlat' ][:,:]
farea[:,:]   = fh_grid.variables['farea'][:,:]
ldist[:,:]   = fh_grid.variables['ldist'][:,:]
ddist[:,:]   = fh_grid.variables['ddist'][:,:]

cube = int(np.sqrt(nface[igrid-1]/6))

fv3flon = np.zeros([6*cube,cube])
fv3flon[:,:] = fh_fv3grid.variables['flons'][:,:]
fv3flat = np.zeros([6*cube,cube])
fv3flat[:,:] = fh_fv3grid.variables['flats'][:,:]

fv3flon = np.reshape(fv3flon, (6*cube*cube))
fv3flat = np.reshape(fv3flat, (6*cube*cube))


fv3vlon = np.zeros([6*cube,cube])
fv3vlon[:,:] = fh_fv3grid.variables['vlons'][:,:]
fv3vlat = np.zeros([6*cube,cube])
fv3vlat[:,:] = fh_fv3grid.variables['vlats'][:,:]

fv3vlon = np.reshape(fv3vlon, (6*cube*cube))
fv3vlat = np.reshape(fv3vlat, (6*cube*cube))


print('Cube =', cube)

fig = plt.figure(figsize=(15, 7.5))

#plt.scatter(flong[igrid-1,0:nface[igrid-1]], flat[igrid-1,0:nface[igrid-1]], s=20, marker='x', c='red')
plt.scatter(vlong[igrid-1,0:nvert[igrid-1]], vlat[igrid-1,0:nvert[igrid-1]], s=20, marker='.', c='blue')

#plt.scatter(fv3flon+0.17, fv3flat, s=20, marker='.', c='red')
plt.scatter(fv3vlon-1.4, fv3vlat, s=20, marker='.', c='red')

plt.savefig('scatter.png')

fh_grid.close()
