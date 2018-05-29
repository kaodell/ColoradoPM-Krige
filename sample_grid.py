# sample_grid.py
#	create a figure of a sample grid with the surrounding kriging sites for Denver area
# Written by: Katelyn O'Dell
# Version history: vi - inital outline 02/01/18
###################################################################################################
#load important modules
import numpy as np
import pylab as pl
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import csv
import datetime as dt
import sys
sys.path.append('/home/kaodell/modules/')
import udf
###################################################################################################

###################################################################################################
# Import grid and surface site locations
###################################################################################################
f = np.load('/home/kaodell/nasa_fires/Martenies_Den3kmKrige/processed_datafiles/sfc_pm15.npz') 
sPM = f['sPM']  
lons = f['lons']
lats = f['lats']
ID = f['ID']
ASites = f['ASites']
f.close()

#get lat lon grid
file = '/home/kaodell/nasa_fires/Martenies_Den3kmKrige/processed_datafiles/geo_em.d03.nc'
nc_fid = Dataset(file, 'r')
glats = nc_fid.variables['XLAT_M'][:]
glons = nc_fid.variables['XLONG_M'][:]
nc_fid.close()

###################################################################################################
# Plot grid and surface sites
###################################################################################################
pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=35,urcrnrlat=44,\
    llcrnrlon=-112,urcrnrlon=-96,lat_ts=10,resolution='l')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
m.shadedrelief()
sx, sy = m(lons, lats)
kx, ky = m(glons[0,:,:], glats[0,:,:])
dx, dy = m(-104.9903,39.7392)
cs = m.pcolor(kx,ky,np.zeros([glons.shape[1],glats.shape[2]]),color = 'grey',alpha=0.1)
m.scatter(sx,sy, color = 'k',s=4)
m.scatter(dx, dy, color = 'r', marker = '*', s = 10)
pl.text(dx, dy, ' Denver', fontsize=11);
pl.title('Colorado kriging grid and surface sites')
pl.savefig('/home/kaodell/nasa_fires/Martenies_Den3kmKrige/figures/grid_w_sfcsites_CO.png')

pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=38.5,urcrnrlat=41,\
    llcrnrlon=-106,urcrnrlon=-103,lat_ts=10,resolution='i')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
m.drawcounties(color = 'r')
#m.shadedrelief()
sx, sy = m(lons, lats)
kx, ky = m(glons[0,:,:], glats[0,:,:])
dx, dy = m(-104.9903,39.7392)
cs = m.pcolor(kx,ky,np.zeros([glons.shape[1],glats.shape[2]]),color = 'grey',alpha=0.1)
m.scatter(sx,sy, color = 'k',s=4)
m.scatter(dx, dy, color = 'r', marker = '*', s = 10)
pl.title('Front range surface sites')
pl.savefig('/home/kaodell/nasa_fires/Martenies_Den3kmKrige/figures/grid_w_sfcsites_FR.png')

pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=39.5,urcrnrlat=40,\
    llcrnrlon=-105.3,urcrnrlon=-104.4,lat_ts=10,resolution='h')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries() 
m.drawcounties(color = 'r') 
#m.drawmapscale(39.7, -105, 0, 0, 50, barstyle='simple', units='km', fontsize=9, yoffset=None, labelstyle='simple', fontcolor='k', fillcolor1='w', fillcolor2='k', ax=None, format='%d', zorder=None)
m.shadedrelief()
sx, sy = m(lons, lats)
kx, ky = m(glons[0,:,:], glats[0,:,:])
dx, dy = m(-104.9903,39.7392)
cs = m.pcolor(kx,ky,np.zeros([glons.shape[1],glats.shape[2]]),color = 'grey',alpha=0.1)
m.scatter(sx,sy, color = 'k',s=10)
m.scatter(dx, dy, color = 'r', marker = '*', s = 20)
#m.arcgisimage(service='ESRI_Imagery_World_2D', verbose= True)
pl.title('Denver county with surrounding monitors')
pl.savefig('/home/kaodell/nasa_fires/Martenies_Den3kmKrige/figures/grid_w_sfcsites_DEN.png')

pl.show()






