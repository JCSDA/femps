import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import math


def midpoint(lat1,lon1_in,lat2,lon2_in):

    lon2 = lon1_in - math.pi
    lon1 = lon2_in - math.pi

    bx = math.cos(lat2) * math.cos(lon2-lon1);
    by = math.cos(lat2) * math.sin(lon2-lon1);
    latm = math.atan2( math.sin(lat1) + math.sin(lat2),
                  math.sqrt( (math.cos(lat1)+bx)*(math.cos(lat1)+bx) + by*by ) )
    lonm = lon1 + math.atan2(by, math.cos(lon1) + bx);

    lonm = lonm + math.pi

    return(latm,lonm)

fig = plt.figure(figsize=(15, 7.5))

c06_file = '/gpfsm/dnb31/drholdaw/JediDev/femps/test/fv3grid_c0006.nc4'
c12_file = '/gpfsm/dnb31/drholdaw/JediDev/femps/test/fv3grid_c0012.nc4'

fh_c06 = Dataset(c06_file)
fh_c12 = Dataset(c12_file)

cube12 = 13
cube06 = 7

fv3vlon_c12 = np.zeros([6*cube12,cube12])
fv3vlon_c12[:,:] = fh_c12.variables['vlons'][:,:]
fv3vlat_c12 = np.zeros([6*cube12,cube12])
fv3vlat_c12[:,:] = fh_c12.variables['vlats'][:,:]

for j in range(6*cube12):
    for i in range(cube12):
        if (fv3vlon_c12[j,i] > 10.0):
            fv3vlon_c12[j,i] = 0.0

for j in range(6*cube12):
    for i in range(cube12):
        if (fv3vlat_c12[j,i] > 10.0):
            fv3vlat_c12[j,i] = 0.0

fv3vlon_c12_plot = np.reshape(fv3vlon_c12, (6*cube12*cube12))
fv3vlat_c12_plot = np.reshape(fv3vlat_c12, (6*cube12*cube12))

plt.scatter(fv3vlon_c12_plot, fv3vlat_c12_plot, s=20, marker='x', c='blue')



fv3vlon_c06 = np.zeros([6*cube06,cube06])
fv3vlon_c06[:,:] = fh_c06.variables['vlons'][:,:]
fv3vlat_c06 = np.zeros([6*cube06,cube06])
fv3vlat_c06[:,:] = fh_c06.variables['vlats'][:,:]

for j in range(6*cube06):
    for i in range(cube06):
        if (fv3vlon_c06[j,i] > 10.0):
            fv3vlon_c06[j,i] = 0.0

for j in range(6*cube06):
    for i in range(cube06):
        if (fv3vlat_c06[j,i] > 10.0):
            fv3vlat_c06[j,i] = 0.0

fv3vlon_c06_plot = np.reshape(fv3vlon_c06, (6*cube06*cube06))
fv3vlat_c06_plot = np.reshape(fv3vlat_c06, (6*cube06*cube06))

plt.scatter(fv3vlon_c06_plot, fv3vlat_c06_plot, s=20, marker='.', c='red')

plt.savefig('scatter_orig.png')


# Bisect the C6 grid to produce the C12 grid
# ------------------------------------------

fv3vlon_c12_bis = np.zeros([6,cube12,cube12])
fv3vlat_c12_bis = np.zeros([6,cube12,cube12])


for t in range(3,4):

    fv3vlon_c06_tile = fv3vlon_c06[t*cube06:(t+1)*cube06]
    fv3vlat_c06_tile = fv3vlat_c06[t*cube06:(t+1)*cube06]

    fv3vlon_c12_tile = fv3vlon_c12[t*cube12:(t+1)*cube12]
    fv3vlat_c12_tile = fv3vlat_c12[t*cube12:(t+1)*cube12]

    for j in range(0,cube06):
        for i in range(0,cube06):
            print(fv3vlon_c06_tile[j,i],fv3vlon_c12_tile[2*j,2*i])

    # Coincidental points
    for j in range(0,cube12,2):
        for i in range(0,cube12,2):
            fv3vlon_c12_bis[t,j,i] = fv3vlon_c06_tile[int(j/2),int(i/2)]
            fv3vlat_c12_bis[t,j,i] = fv3vlat_c06_tile[int(j/2),int(i/2)]

    if (t==0 or t==3 or t==4):

        # UP/DOWN direction
        for j in range(0,cube12,2):
            for i in range(1,cube12-1,2):
                (fv3vlat_c12_bis[t,j,i], fv3vlon_c12_bis[t,j,i]) = midpoint( \
                  fv3vlat_c06_tile[int(j/2),int((i-1)/2)], fv3vlon_c06_tile[int(j/2),int((i-1)/2)], \
                  fv3vlat_c06_tile[int(j/2),int((i+1)/2)], fv3vlon_c06_tile[int(j/2),int((i+1)/2)] )

                #print(fv3vlat_c12_bis[t,j,i]-fv3vlat_c12[t*cube12+j,i])

        # LEFT/RIGHT direction
        for j in range(1,cube12-1,2):
            for i in range(0,cube12):
                (fv3vlat_c12_bis[t,j,i], fv3vlon_c12_bis[t,j,i]) = midpoint( \
                  fv3vlat_c12_bis[t,j+1,i], fv3vlon_c12_bis[t,j+1,i], \
                  fv3vlat_c12_bis[t,j-1,i], fv3vlon_c12_bis[t,j-1,i]  )#

#        for j in range(0,cube12):
#            for i in range(0,cube12):
#                fv3vlat_c12_bis[t,j,i] = fv3vlat_c12[t*cube12+j,i]
#                fv3vlon_c12_bis[t,j,i] = fv3vlon_c12[t*cube12+j,i]

#        for j in range(1,cube12-1,2):
#            for i in range(1,cube12-1,2):
#                fv3vlon_c12_bis[t,j,i] = 0.5*(fv3vlon_c12_bis[t,j,i+1] + fv3vlon_c12_bis[t,j,i])
#                fv3vlat_c12_bis[t,j,i] = 0.5*(fv3vlat_c12_bis[t,j,i+1] + fv3vlat_c12_bis[t,j,i])#

#        for j in range(1,cube12-1,2):
#            for i in range(1,cube12-1,2):
#                fv3vlon_c12_bis[t,j,i] = 0.5*(fv3vlon_c12_bis[t,j,i+1] + fv3vlon_c12_bis[t,j,i])
#                fv3vlat_c12_bis[t,j,i] = 0.5*(fv3vlat_c12_bis[t,j,i+1] + fv3vlat_c12_bis[t,j,i])

    #if (t==1 or t==2 or t==5):
    #    for j in range(0,cube12-1,2):
    #        for i in range(1,cube12,2):
    #            fv3vlon_c12_bis[t,j,i] = 0.5*(fv3vlon_c06_tile[int((j+1)/2),int((i)/2)] + fv3vlon_c06_tile[int((j-1)/2),int((i)/2)])
    #            fv3vlat_c12_bis[t,j,i] = 0.5*(fv3vlat_c06_tile[int((j+1)/2),int((i)/2)] + fv3vlat_c06_tile[int((j-1)/2),int((i)/2)])


fv3vlon_c12_bis_plot = np.reshape(fv3vlon_c12_bis, (6*cube12*cube12))
fv3vlat_c12_bis_plot = np.reshape(fv3vlat_c12_bis, (6*cube12*cube12))

fv3vlon_c06_tile_plot = np.reshape(fv3vlon_c06_tile, (cube06*cube06))
fv3vlat_c06_tile_plot = np.reshape(fv3vlat_c06_tile, (cube06*cube06))


fig = plt.figure(figsize=(15, 7.5))

c_start = 3
c_final = 3

pt06_start = c_start*cube06*cube06
pt06_final = c_start*cube06*cube06 + cube06*cube06

pt12_start = c_start*cube12*cube12
pt12_final = c_start*cube12*cube12 + cube12*cube12

#plt.scatter(fv3vlon_c06_plot[pt06_start:pt06_final], fv3vlat_c06_plot[pt06_start:pt06_final], s=20, marker='o', c='green')
plt.scatter(fv3vlon_c06_tile_plot, fv3vlat_c06_tile_plot, s=20, marker='o', c='green')
plt.scatter(fv3vlon_c12_bis_plot[pt12_start:pt12_final], fv3vlat_c12_bis_plot[pt12_start:pt12_final], s=20, marker='x', c='blue')
plt.scatter(fv3vlon_c12_plot[pt12_start:pt12_final], fv3vlat_c12_plot[pt12_start:pt12_final], s=20, marker='.', c='red')

plt.savefig('scatter_bis.png')

fh_c06.close()
fh_c12.close()
