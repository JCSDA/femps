import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# User input
# ----------
gridfile = '/discover/nobackup/drholdaw/JediDev/femps/build-int-release/test/grid.nc4'

igrid = 1 # Grid number to plot


# File handle
# -----------
fh_grid = Dataset(gridfile)


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

nface = np.zeros([ngrids])
nedge = np.zeros([ngrids])
nvert = np.zeros([ngrids])
neoff = np.zeros([ngrids, nfacex])
neofv = np.zeros([ngrids, nvertx])
fnxtf = np.zeros([ngrids, dimfnxtf, nfacex])
eoff  = np.zeros([ngrids, dimeoff, nfacex])
voff  = np.zeros([ngrids, dimvoff, nfacex])
fnxte = np.zeros([ngrids, dimfnxte, nedgex])
vofe  = np.zeros([ngrids, dimvofe, nedgex])
fofv  = np.zeros([ngrids, dimfofv, nvertx])
eofv  = np.zeros([ngrids, dimeofv, nvertx])
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

fig = plt.figure(figsize=(30, 15))

plt.scatter(flong[igrid-1,:], flat[igrid-1,:], s=6, marker='.', c='C0')
plt.scatter(vlong[igrid-1,:], vlat[igrid-1,:], s=6, marker='x', c='C4')

plt.savefig('scatter.pdf')

fh_grid.close()
