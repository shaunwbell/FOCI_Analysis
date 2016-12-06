#System Stack
import csv


#Science Stack
import numpy as np
from netCDF4 import Dataset

# Visual Stack
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

def etopo5_data():
    """ read in etopo5 topography/bathymetry. """
    file = '/Users/bell/Data_Local/MapGrids/etopo5.nc'
    etopodata = Dataset(file)
    
    topoin = etopodata.variables['bath'][:]
    lons = etopodata.variables['X'][:]
    lats = etopodata.variables['Y'][:]
    etopodata.close()
    
    topoin,lons = shiftgrid(0.,topoin,lons,start=False) # -360 -> 0
    
    lons, lats = np.meshgrid(lons, lats)
    
    return(topoin, lats, lons)
    
d = {}
count = 0
with open('/Users/bell/Data_Local/Reanalysis_Files/ncepstorms_2000_2005.txt','rb') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    for row in tsvin:
        d[count] = row[0].strip().split()
        count = count + 1
        
lat = np.array([ d[k][13] for i,k in enumerate(d.keys())], float)
lon = -1 * np.array([ d[k][14] for i,k in enumerate(d.keys())], float) + 180.
year = np.array([ d[k][2] for i,k in enumerate(d.keys())], float)
boxnumlat = np.array([ d[k][-5] for i,k in enumerate(d.keys())], int)
boxnumlon = np.array([ d[k][-4] for i,k in enumerate(d.keys())], int)
idnum = np.array([ d[k][-1] for i,k in enumerate(d.keys())], int)

color = ['b','g','r','y','c']

freq_array = np.zeros(shape=(70,70))
lat_array = np.zeros(shape=(70,70))
lon_array = np.zeros(shape=(70,70))


    
(topoin, elats, elons) = etopo5_data()

fig = plt.figure()
ax = plt.subplot(111)
m = Basemap(resolution='c',projection='merc', llcrnrlat=50, \
    urcrnrlat=65,llcrnrlon=-175,urcrnrlon=-130, lat_ts=45)


#ETOPO 5 contour data 
ex, ey = m(elons, elats)
CS = m.contourf(ex,ey,topoin, levels=range(250,5000,250), cmap='gray_r', alpha=.75) #colors='black'
CS = m.contour(ex,ey,topoin, levels=range(250,5000,250), linewidths=0.2, colors='black', alpha=.75) #
CS = m.contour(ex,ey,topoin, levels=[-1000, -200, -100], linestyle='--', linewidths=0.2, colors='black', alpha=.75) #
plt.clabel(CS, inline=1, fontsize=8, fmt='%1.0f')


#plot points
for yy in (0,):
    for s in range(1,idnum.max(),1):
        s_ind = (idnum[year == yy] == s)
        lat_p = lat[s_ind]
        lon_p = lon[s_ind]
        if all(lon_p > -130):
            continue
        if all(lat_p < 45):
            continue
        x, y = m(lon_p, lat_p)
        m.plot(x, y,'r')

m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(46,66,4.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
m.drawmeridians(np.arange(-180,-140,5.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
#m.fillcontinents(color='black')

DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*1.5) )