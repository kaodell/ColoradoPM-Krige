#eval_kriging.py
#	A python script to make figures to evaluate kriged data, written to work with Denver kriging output
#written by: Kate O'Dell
#version history: 1.0 initial 04/10/2018 
###################################################################################################
# import modules
###################################################################################################
import numpy as np
import pylab as pl
from mpl_toolkits.basemap import Basemap, cm
import datetime as dt
import sys
sys.path.append('/home/kaodell/modules/')
import udf
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import linregress
###################################################################################################
# set year, load data and make datetime array
###################################################################################################
year = 2010

t0 = dt.datetime( year = year, month = 1, day = 1 )
tf = dt.datetime( year = year, month = 12, day = 31 )
NT = (tf - t0).days
alldays = np.empty( NT+1, dtype = object )
alldays[0] = t0
for i in range(1, NT+1):
	alldays[i] = alldays[i-1] + dt.timedelta(days = 1)

fname =  '/home/kaodell/nasa_fires/Sheena_Den3kmKrige/updated/kpmDEN'+ str(year) + '_3km_update.npz'
fid = np.load(fname,'r')
kPM = fid['kPM']
sPMin = fid['sPMin']
sPMout = fid['sPMout']
glons = fid['glons']
glats = fid['glats']
ilons = fid['ilon']
ilats = fid['ilat']
olons = fid['olon']
olats = fid['olat']
rsq = fid['rsq']
MB = fid['MB']
MAE = fid['MAE']
slope = fid['slope']
kPMloocv = fid['kPM_site']
fid.close()

###################################################################################################
# plot a time series of kriged values and site values for a site
###################################################################################################

# plot 1-1 of krige and sPM
inds = []
pl.figure()
for si in range(sPMin.shape[1]):
	slat = ilats[si]
	slon = ilons[si]
	dists = udf.haversine(slon,glons,slat,glats)
	ind = np.where(dists == dists.min())
	inds.append(ind)
	kPMsite = kPM[:,ind[0],ind[1]]
	pl.plot(sPMin[:,si],kPMsite,'.')
pl.plot(np.arange(0,60),np.arange(0,60),linewidth=4.0,color='black')
pl.xlabel('site PM')
pl.ylabel('kriged PM')
pl.title('Check on krige PM and site PM for Denver sites')

# plot 1-1 of LOOCV krige and sPM
inds = []
pl.figure()
for si in range(sPMin.shape[1]):
	pl.plot(sPMin[:,si],kPMloocv[:,si],'.')
pl.plot(np.arange(0,60),np.arange(0,60),linewidth=4.0,color='black')
pl.xlabel('site PM')
pl.ylabel('kriged PM, loocv')
pl.title('Check on loocv krige PM and site PM for Denver sites')
pl.text(2,59,'rsq: ' + str(rsq)[:4])
pl.text(2,56,'slope: ' + str(slope)[:4])
pl.text(2,53,'MB: ' + str(MB)[:4])
pl.text(2,50,'MAE: ' + str(MAE)[:4])

'''
###################################################################################################
# plot surface obs with kriging for a few days
###################################################################################################
di = 160

pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=35,urcrnrlat=44,\
    llcrnrlon=-110,urcrnrlon=-96,lat_ts=10,resolution='l')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
isx, isy = m(ilons, ilats)
osx, osy = m(olons,olats)
kx, ky = m(glons, glats)
cs = m.pcolor(kx,ky,kPM[di,:,:], vmin = -1, vmax = 15)
m.scatter(isx,isy, c = sPMin[di,:] ,s=20, edgecolor = 'black', vmin = -1, vmax = 15)
m.scatter(osx,osy, color = 'black', s = 20, vmin = -1, vmax = 15)
cbar = m.colorbar(cs,location='bottom',pad="10%", ticks = [-5,0,5,10,15,20,25])
cbar.set_label('PM$_{2.5}$ [$\mu$g/m$^3$]')
pl.title('Colorado kriging grid and surface sites')
pl.savefig('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/figures/' + str(year) + 'COmtncutoff.png')
#pl.show()

'''
###################################################################################################
# determine r2 at each site over time period and plot
###################################################################################################
rsq = np.zeros(sPMin.shape[1])
numdays = np.zeros(sPMin.shape[1])
for si in range(sPMin.shape[1]):
	trueinds = np.where(np.isnan(sPMin[:,si])==False)
	dsx = sPMin[trueinds,si]
	dsy = kPMloocv[trueinds,si]

	numdays[si] = len(dsx[0])	

	trueinds2 = np.where(np.isnan(dsy) == False)
	dsx = dsx[trueinds2]
	dsy = dsy[trueinds2]


	slope, intercept, rvalue, pvalue, stderr = linregress( dsx,dsy )
	rsq[si] = rvalue**2.
'''
pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=36,urcrnrlat=41.5,\
    llcrnrlon=-110,urcrnrlon=-100,lat_ts=10,resolution='l')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
isx, isy = m(ilons, ilats)
cs = m.scatter(isx,isy, c = rsq ,s=40, edgecolor = 'black', vmin = 0, vmax = 1.)
cbar = m.colorbar(cs,location='bottom',pad="10%", ticks = [0,0.25,0.5,0.75,1.0])
cbar.set_label('R$^2$')
pl.title('Colorado kriging loocv R$^2$ at sites ' + str(year))
pl.savefig('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/figures/' + str(year) + 'COmtncutoff_site_r2.png')
#pl.show()

	
###################################################################################################
# determine # days availabe at each site over time period and plot
###################################################################################################

pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=36,urcrnrlat=41.5,\
    llcrnrlon=-110,urcrnrlon=-100,lat_ts=10,resolution='l')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
isx, isy = m(ilons, ilats)
cs = m.scatter(isx,isy, c = numdays ,s=30, edgecolor = 'black', vmin = 0, vmax = 270.,cmap = 'inferno')
cbar = m.colorbar(cs,location='bottom',pad="10%", ticks = [0,45,90,135,180,225,270])
cbar.set_label('# days')
pl.title('Colorado # avail days at sites ' + str(year))
pl.savefig('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/figures/' + str(year) + 'COmtncutoff_site_ndays.png')
#pl.show()
'''

###################################################################################################
# plot scatter plot of R2 vs # of days for each site, colored by site location
###################################################################################################

# regions to color sites by
DENinds = np.where((ilats > 39.5) & (ilats < 40))
NOCO = np.where((ilats > 40) & (ilats < 41))
WY = np.where(ilats > 41)
CS_P =  np.where((ilats > 38) & (ilats < 39))
SOCO = np.where((ilats > 37) & (ilats < 38))

# make figure
fig,ax = pl.subplots()
#pl.plot(ilats, numdays,'.')
pl.scatter(numdays[DENinds], rsq[DENinds],marker = 'o', facecolor = 'none',edgecolor = 'cornflowerblue', s = 50, label = 'DEN')
pl.scatter(numdays[NOCO], rsq[NOCO],marker = 'o', facecolor = 'none',edgecolor = 'purple', s = 50, label = 'NOCO')
pl.scatter(numdays[CS_P], rsq[CS_P] ,marker = 'o',facecolor = 'none', edgecolor = 'indianred', s = 50, label = 'CS and Pueblo')
pl.scatter(numdays[WY], rsq[WY] ,marker = 'o', facecolor='none',edgecolor = 'forestgreen', s = 50, label = 'WY')
pl.legend(loc='best')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
pl.xlabel('# days with data',fontsize=12)
pl.ylabel('R$^2$',fontsize=12)
pl.title('R$^2$ and # avail days at sites, ' + str(year))
pl.savefig('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/updated/' + str(year) + 'COmtncutoff_rsq_ndays.png')
#pl.show()

# map of site locations
pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=36,urcrnrlat=41.5,\
    llcrnrlon=-110,urcrnrlon=-100,lat_ts=10,resolution='l')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
DENx, DENy = m(ilons[DENinds], ilats[DENinds])
NOCOx, NOCOy = m(ilons[NOCO],ilats[NOCO])
WYx, WYy = m(ilons[WY],ilats[WY])
CS_Px, CS_Py = m(ilons[CS_P],ilats[CS_P])
cs = m.scatter(DENx,DENy, marker = 'o', facecolor = 'none',edgecolor = 'cornflowerblue' ,s=55)
cs = m.scatter(NOCOx,NOCOy, marker = 'o', facecolor = 'none',edgecolor = 'purple' ,s=55)
cs = m.scatter(WYx,WYy, marker = 'o', facecolor = 'none',edgecolor = 'forestgreen' ,s=55)
cs = m.scatter(CS_Px,CS_Py, marker = 'o', facecolor = 'none',edgecolor = 'indianred' ,s=55)
pl.title('Colorado sites ' + str(year))
pl.savefig('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/updated/' + str(year) + 'COmtncutoff_sites.png')
pl.show()



