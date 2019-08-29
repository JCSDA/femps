import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import math
import sys
import os
import shutil as sh

# User input
# ----------
path = sys.argv[1]

gridfile = path+'/griddata/grid.nc4'
oprsfile = path+'/griddata/operators.nc4'

ig = 0 # 0-5  Grid number to plot


# File handle
# -----------
fh_grid = Dataset(gridfile)
fh_oprs = Dataset(oprsfile)


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

elong = np.zeros([ngrids, nedgex])
elat  = np.zeros([ngrids, nedgex])

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

elong[:,:]   = fh_oprs.variables['elong'][:,:]
elat [:,:]   = fh_oprs.variables['elat' ][:,:]

# Plot the grid
# -------------

# Plotting choices
limoff = 0.1
txtoff = 0.01
fontsize = 9

fmark = 'x'
vmark = '.'
emark = 'v'

sbg = 20
sfg = 2*sbg

#Clean up paths
sh.rmtree(path+'/plots/fnxtf', ignore_errors=True)
sh.rmtree(path+'/plots/eoff', ignore_errors=True)
sh.rmtree(path+'/plots/voff', ignore_errors=True)
sh.rmtree(path+'/plots/fnxte', ignore_errors=True)
sh.rmtree(path+'/plots/vofe', ignore_errors=True)
sh.rmtree(path+'/plots/fofv', ignore_errors=True)
sh.rmtree(path+'/plots/eofv', ignore_errors=True)

os.mkdir( path+'/plots/fnxtf');
os.mkdir( path+'/plots/eoff');
os.mkdir( path+'/plots/voff');
os.mkdir( path+'/plots/fnxte');
os.mkdir( path+'/plots/vofe');
os.mkdir( path+'/plots/fofv');
os.mkdir( path+'/plots/eofv');

def plot_all():
    plt.figure(figsize=(15,7.5))
    plt.scatter(flong[ig,0:nface[ig]], flat[ig,0:nface[ig]], s=sbg, marker=fmark, c='grey')
    plt.scatter(vlong[ig,0:nvert[ig]], vlat[ig,0:nvert[ig]], s=sbg, marker=vmark, c='black')
    plt.scatter(elong[ig,0:nedge[ig]], elat[ig,0:nedge[ig]], s=sbg, marker=emark, c='grey')
    return

def plot_save(savepath):
    plt.xlim(0.0-limoff, 2*math.pi+limoff)
    plt.ylim(-math.pi/2-limoff, math.pi/2+limoff)
    plt.savefig(savepath)
    plt.close()
    return

# fnxtf
# -----
print('fnxtf')
for i1 in range(nface[ig]):
    plot_all()
    i1str = str(i1).zfill(3)
    plt.scatter(flong[ig,i1], flat[ig,i1], s=sfg, marker=fmark, c='blue')

    for i2 in range(dimfnxtf):
        i2str = str(i2).zfill(1)
        ind = fnxtf[ig,i2,i1]-1
        x = flong[ig,ind]
        y =  flat[ig,ind]
        plt.scatter(x, y, s=sfg, marker=fmark, c='red')
        plt.text(x+txtoff, y+txtoff, i2str, fontsize=fontsize)

    plot_save(path+'/plots/fnxtf/scatter_'+i1str+'.png')

# eoff
# ----
print('eoff')
for i1 in range(nface[ig]):
    plot_all()
    i1str = str(i1).zfill(3)
    plt.scatter(flong[ig,i1], flat[ig,i1], s=sfg, marker=fmark, c='blue')

    for i2 in range(dimeoff):
        i2str = str(i2).zfill(1)
        ind = eoff[ig,i2,i1]-1
        x = elong[ig,ind]
        y =  elat[ig,ind]
        plt.scatter(x, y, s=sfg, marker=emark, c='red')
        plt.text(x+txtoff, y+txtoff, i2str, fontsize=fontsize)

    plot_save(path+'/plots/eoff/scatter_'+i1str+'.png')

# voff
# ----
print('voff')
for i1 in range(nface[ig]):
    plot_all()
    i1str = str(i1).zfill(3)
    plt.scatter(flong[ig,i1], flat[ig,i1], s=sfg, marker=fmark, c='blue')

    for i2 in range(dimvoff):
        i2str = str(i2).zfill(1)
        ind = voff[ig,i2,i1]-1
        x = vlong[ig,ind]
        y =  vlat[ig,ind]
        plt.scatter(x, y, s=sfg, marker=vmark, c='red')
        plt.text(x+txtoff, y+txtoff, i2str, fontsize=fontsize)

    plot_save(path+'/plots/voff/scatter_'+i1str+'.png')

# fnxte
# -----
print('fnxte')
for i1 in range(nedge[ig]):
    plot_all()
    i1str = str(i1).zfill(3)
    plt.scatter(elong[ig,i1], elat[ig,i1], s=sfg, marker=emark, c='blue')

    for i2 in range(dimfnxte):
        i2str = str(i2).zfill(1)
        ind = fnxte[ig,i2,i1]-1
        x = flong[ig,ind]
        y =  flat[ig,ind]
        plt.scatter(x, y, s=sfg, marker=fmark, c='red')
        plt.text(x+txtoff, y+txtoff, i2str, fontsize=fontsize)

    plot_save(path+'/plots/fnxte/scatter_'+i1str+'.png')

# vofe
# ----
print('vofe')
for i1 in range(nedge[ig]):
    plot_all()
    i1str = str(i1).zfill(3)
    plt.scatter(elong[ig,i1], elat[ig,i1], s=sfg, marker=emark, c='blue')

    for i2 in range(dimvofe):
        i2str = str(i2).zfill(1)
        ind = vofe[ig,i2,i1]-1
        x = vlong[ig,ind]
        y =  vlat[ig,ind]
        plt.scatter(x, y, s=sfg, marker=vmark, c='red')
        plt.text(x+txtoff, y+txtoff, i2str, fontsize=fontsize)

    plot_save(path+'/plots/vofe/scatter_'+i1str+'.png')

# fofv
# ----
print('fofv')
for i1 in range(nvert[ig]):
    plot_all()
    i1str = str(i1).zfill(3)
    plt.scatter(vlong[ig,i1], vlat[ig,i1], s=sfg, marker=vmark, c='blue')

    for i2 in range(dimfofv):
        i2str = str(i2).zfill(1)
        ind = fofv[ig,i2,i1]-1
        x = flong[ig,ind]
        y =  flat[ig,ind]
        plt.scatter(x, y, s=sfg, marker=fmark, c='red')
        plt.text(x+txtoff, y+txtoff, i2str, fontsize=fontsize)

    plot_save(path+'/plots/fofv/scatter_'+i1str+'.png')

# eofv
# ----
print('eofv')
for i1 in range(nvert[ig]):
    plot_all()
    i1str = str(i1).zfill(3)
    plt.scatter(vlong[ig,i1], vlat[ig,i1], s=sfg, marker=vmark, c='blue')

    for i2 in range(dimeofv):
        i2str = str(i2).zfill(1)
        ind = eofv[ig,i2,i1]-1
        x = elong[ig,ind]
        y =  elat[ig,ind]
        plt.scatter(x, y, s=sfg, marker=emark, c='red')
        plt.text(x+txtoff, y+txtoff, i2str, fontsize=fontsize)

    plot_save(path+'/plots/eofv/scatter_'+i1str+'.png')

fh_grid.close()
