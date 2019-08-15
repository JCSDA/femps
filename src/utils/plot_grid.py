import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import math
import sys


# User input
# ----------
path = sys.argv[1]

gridfile = path+'/griddata/grid.nc4'

ig = 0 # 0-5  Grid number to plot


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

nface = np.zeros([ngrids], dtype=int)
nedge = np.zeros([ngrids], dtype=int)
nvert = np.zeros([ngrids], dtype=int)
neoff = np.zeros([ngrids, nfacex], dtype=int)
neofv = np.zeros([ngrids, nvertx], dtype=int)
fnxtf = np.zeros([ngrids, dimfnxtf, nfacex], dtype=int)
eoff  = np.zeros([ngrids, dimeoff,  nfacex], dtype=int)
voff  = np.zeros([ngrids, dimvoff,  nfacex], dtype=int)
fnxte = np.zeros([ngrids, dimfnxte, nedgex], dtype=int)
vofe  = np.zeros([ngrids, dimvofe,  nedgex], dtype=int)
fofv  = np.zeros([ngrids, dimfofv,  nvertx], dtype=int)
eofv  = np.zeros([ngrids, dimeofv,  nvertx], dtype=int)
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


# Plot the grid
# -------------

def plot_all(savepath):
    plt.scatter(flong[ig,0:nface[ig]], flat[ig,0:nface[ig]], s=3, marker='x', c='grey')
    plt.scatter(vlong[ig,0:nvert[ig]], vlat[ig,0:nvert[ig]], s=3, marker='.', c='black')
    plt.xlim(0.0-limoff, 2*math.pi+limoff)
    plt.ylim(-math.pi/2-limoff, math.pi/2+limoff)
    plt.savefig(savepath)
    plt.close()
    return

# Plotting choices
limoff = 0.1


# fnxtf
# -----

for f1 in range(nface[ig]):
    vstr = str(f1).zfill(3)
    plt.scatter(flong[ig,f1], flat[ig,f1], s=20, marker='x', c='blue')

    for f2 in range(dimfnxtf):
        plt.scatter(flong[ig,fnxtf[ig,f2,f1]-1], flat[ig,fnxtf[ig,f2,f1]-1], s=20, marker='*', c='red')

    plot_all(path+'/plots/fnxtf/scatter_'+vstr+'.png')

# eoff
# ----

# voff
# ----

for f1 in range(nface[ig]):
    vstr = str(f1).zfill(3)
    plt.scatter(flong[ig,f1], flat[ig,f1], s=20, marker='x', c='blue')

    for v1 in range(dimvoff):
        plt.scatter(vlong[ig,voff[ig,v1,f1]-1], vlat[ig,voff[ig,v1,f1]-1], s=20, marker='*', c='red')

    plot_all(path+'/plots/voff/scatter_'+vstr+'.png')

# fnxte
# -----

# vofe
# ----

# fofv
# ----
for v1 in range(nvert[ig]):
    vstr = str(v1).zfill(3)
    plt.scatter(vlong[ig,v1], vlat[ig,v1], s=20, marker='x', c='blue')

    for f1 in range(dimfofv):
        plt.scatter(flong[ig,fofv[ig,f1,v1]-1], flat[ig,fofv[ig,f1,v1]-1], s=20, marker='*', c='red')

    plot_all(path+'/plots/fofv/scatter_'+vstr+'.png')

# eofv
# ----


fh_grid.close()
