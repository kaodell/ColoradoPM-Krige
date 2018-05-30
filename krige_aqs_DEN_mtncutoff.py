#Krige_aqs_DEN
#	Reading aqs data from multiple years and kriging to trimmed Colorado 3km wrf-chem grid (truncating monitors used at the mountains) 2010 - 2014
#written by: Kate O'Dell 
#version history: 1.0 initial 09/30/16 
#		 2.0 re-written for 2015 data 10/25/16
#		 3.0 re-written for 2013 Oregon (renamed from krige_aqs_2015 to krige_aqs_OR2013) 1/16/17
#		 4.0 re-written for Denver on a 3km grid (renamed from krige_aqs_OR2013 to krige_aqs_DEN) 1/29/18
#		 4.5 edited to cut grid at the mountains and throw out any sites to the west (ie in the mountains)
########################################################################
#load important modules
import numpy as np
import pylab as pl
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import csv
import datetime as dt
from pykrige.ok import OrdinaryKriging
from scipy.stats import linregress
import sys
sys.path.append('/home/kaodell/modules/')
import udf
########################################################################

if len(sys.argv) > 1:
	print sys.argv[1]
	#sys.exit()
	year = int(sys.argv[1])
else:
	print 'no'

#year = 2016

#Time period of interest
t0 = dt.datetime( year = year, month = 1, day = 1 )
tf = dt.datetime( year = year, month = 12, day = 31 )
NT = (tf - t0).days
alldays = np.empty( NT+1, dtype = object )
alldays[0] = t0
for i in range(1, NT+1):
	alldays[i] = alldays[i-1] + dt.timedelta(days = 1)

############################################################################################################
#Part 1-2: makeing sPM arrays from EPA aqs data, already done
#	see earlier versions of code for the calculations
# this section pulls in daily PM2.5 concentrations from the US EPA AQS datamart (available on the EPA AQS site)
# and organizes the input into an array called sPM (for site PM2.5) that is number of sites x number of days
# with a separate array for each year
############################################################################################################

############################################################################################################
#Part 3: load in data and process
############################################################################################################

#f = np.load('/home/kaodell/nasa_fires/processed_datafiles/US_multiyear_krige/sfc_pm'+str(year)[2:]+'.npz') 
f = np.load('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/updated/sfc_pm'+str(year)[2:]+'_wCDPHE.npz') 
sPM = f['sPM']  
lons = f['lons']
lats = f['lats']
#ID = f['ID']
#ASites = f['ASites']
f.close()

# get lat lon grid - this is the grid we are kriging to, it is the grid used by the WRF-Chem model
file = '/home/kaodell/nasa_fires/Sheena_Den3kmKrige/processed_datafiles/geo_em.d03.nc'
nc_fid = Dataset(file, 'r')
glats = nc_fid.variables['XLAT_M'][:][0,:,:]
glons = nc_fid.variables['XLONG_M'][:][0,:,:]
nc_fid.close()	

NT, NS = sPM.shape
NX, NY = glons.shape


# trim grid to Colorado
# Boundaries
LatMax = 47.
LatMin = 30.

LonMax = -60.
LonMin = -105.1


#Trim WRF Grid/ WRF 
tAlat, tBlat = np.where( glats < LatMax )
tClat, tDlat = np.where( glats > LatMin )

tAlon, tBlon = np.where( glons < LonMax )
tClon, tDlon = np.where( glons > LonMin )

wlats = glats[tClat.min():tAlat.max(), tDlon.min():tBlon.max()]
wlons = glons[tClat.min():tAlat.max(), tDlon.min():tBlon.max()]

glats = wlats
glons = wlons

NT, NS = sPM.shape
NX, NY = glons.shape


############################################################################################################
#Part 4: Trim and sort data
############################################################################################################

# sPMin is sPM for sites inside glon and glat grid  
# sPMout sPM for sites close to glon and glat (but not inside) we will use these for the kriging but not evaluate
# the sPM array at these sites

NT, NS = sPM.shape

dmax_in = 15000. #m
dmax_out = 150000. #m

iIND = []
oIND = []
# now loop through sites and find those inside the grid and those close enough to use for kriging
for siteIND in range(0,NS):
	# use the haversine formula (https://en.wikipedia.org/wiki/Haversine_formula) to determine distances
	# between sites and the grid
	dists = udf.haversine(lons[siteIND],glons, lats[siteIND], glats)
	dmin = dists.min()
	if dmin <= dmax_in:
		iIND.append(siteIND)
	elif dmin <= dmax_out:
		oIND.append(siteIND)

sPMin = sPM[:,iIND]
sPMout = sPM[:,oIND]
ilat = lats[iIND]
ilon = lons[iIND]
olat = lats[oIND]
olon = lons[oIND]

dinds = []
# take out osites to the West, these are in the mountains
for i in range(len(olon)):
	if olon[i] < -105:
		dinds.append(i)
dinds = np.array(dinds)

olon = np.delete(olon, dinds)
olat = np.delete(olat, dinds)
sPMout = np.delete(sPMout,dinds,axis = 1)

NWS = len(iIND)
NNWS = len(oIND)
'''
# plot kriging sites in and out
pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=35,urcrnrlat=44,\
    llcrnrlon=-112,urcrnrlon=-96,lat_ts=10,resolution='l')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
m.shadedrelief()
isx, isy = m(ilon, ilat)
osx, osy = m(olon,olat)
kx, ky = m(glons, glats)
cs = m.pcolor(kx,ky,np.zeros([glons.shape[0],glats.shape[1]]),color = 'grey',alpha=0.1)
m.scatter(isx,isy, color = 'r',s=4)
m.scatter(osx,osy, color = 'k',s=4)
pl.title('Colorado kriging grid and surface sites')
pl.savefig('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/figures/3kmgrid_w_krigesites_COmtncutoff.png')
'''

############################################################################################################
#Part 5: LOOCV - Leave One Out Cross Validation
############################################################################################################

# Loop though each site, remove that site (and co-located sites), and krige the remaining sites. 
# Then evaluate the estimated PM2.5 at the removed monitor location against the montior-observed value for each day.

#define variogram parameters
#vgp = [3.0, 8.0, 0.1]		#Wash 2012 parameters
#vgp = [3.0,2.5,0.1]		#Oregon 2013 parameters, tested. 
vgp = [2.8,2.5,0.2]		#variogram parameters tested for Eastern Colorado monitor set-up

# for testing variogram parameters
#if len(sys.argv) > 1:
#	print sys.argv[1], sys.argv[2], sys.argv[3]
	#sys.exit()
#	s = float(sys.argv[1])
#	r = float(sys.argv[2])
#	n = float(sys.argv[3])
#	vgp = [s, r, n]
#else:
#	print 'no'

#preallocate matricies
dsx = np.array([-999])
dsy = np.array([-999])

kRSQ = np.empty(NWS)

kPM_site = np.zeros([NT,NWS])
INVALIDS = []				#where sPMin is over 60% nans
yout_all = np.zeros([NT,NWS])
for si in range(0,NWS):
	#remove site being tested
	sPMinTEMP = np.delete(sPMin, si, axis = 1)
	ilatTEMP = np.delete(ilat, si)
	ilonTEMP = np.delete(ilon, si)

	#save data from site being tested
	dsx = np.hstack([dsx, sPMin[:, si]])
	clat = ilat[si]
	clon = ilon[si]

	#recombine sites
	sPMtemp = np.hstack([sPMinTEMP, sPMout])
	latTEMP = np.hstack([ilatTEMP, olat])
	lonTEMP = np.hstack([ilonTEMP, olon])
	
	#remove sites with same lat/lon as test site
	dups = np.where( latTEMP == ilat[si])[0]
	sPMtemp2 = np.delete(sPMtemp, dups, axis = 1)
	latTEMP2 = np.delete(latTEMP, dups)
	lonTEMP2 = np.delete(lonTEMP,dups)

	#loop through time
	yout = np.empty(NT)
	
	for ti in range(NT):		
		sPMtemp3 = sPMtemp2[ti,:]
		
		#find and remove all nans				
		trueINDS = np.where(np.isnan(sPMtemp3) == False )[0]
		sPMtemp4 = sPMtemp3[trueINDS]
		latTEMP3 = latTEMP2[trueINDS]
		lonTEMP3 = lonTEMP2[trueINDS]

		if len(lonTEMP3) <= 1:
			yout[ti] = np.nan
		else:

			#krige using pykrige function
			#print lonTEMP3.shape, latTEMP3.shape, sPMtemp4.shape
			ok = OrdinaryKriging( lonTEMP3, latTEMP3, sPMtemp4, variogram_model='spherical', \
				variogram_parameters = vgp, verbose=False, enable_plotting=False)
			z, ss = ok.execute('points', clon, clat)
			yout[ti] = z
	
	dsy = np.hstack([dsy, yout])
	
	xout = sPMin[:,si] 
	kPM_site[:,si] = yout
	nanINDS = np.where( np.isnan(xout) == True )[0]		
	xout = np.delete(xout, nanINDS)
	yout = np.delete(yout, nanINDS)

#remove invalids (-999s and nans)	
dsx = dsx[1:]
dsy = dsy[1:]

trueINDS1 = np.where(np.isnan(dsx) == False)[0]
dsx = dsx[trueINDS1]
dsy = dsy[trueINDS1]

trueINDS2 = np.where(np.isnan(dsy) == False)[0]
dsx = dsx[trueINDS2]
dsy = dsy[trueINDS2]

# calculate stats
slope, intercept, rvalue, pvalue, stderr = linregress( dsx,dsy )
rsq = rvalue**2.

MB = np.mean(dsy - dsx)
MAE = np.mean(np.abs(dsy - dsx))

#print LOOCV stats to log file when testing svp
#print rsq, slope, MB, MAE

############################################################################################################
#Part 6: Krige to grid
############################################################################################################
lonK = np.hstack( [ilon, olon] )
latK = np.hstack( [ilat, olat] )

# preallocate for kriging output
kPM = np.empty( [NT, NX, NY] )
NS = NWS + NNWS
NED = 0
EDind = []
for ti in range(NT):
	datatemp = np.hstack([sPMin[ti,:], sPMout[ti,:]])
	INDS = np.where( np.isnan( datatemp ) == False )
	if len(INDS[0]) <= 0.3 * NS:
		#print "Only 30% sites available. Skipping" % len(INDS[0])
		#print NED
		NED += 1
		EDind.append(ti)
		#continue
	datatemp = datatemp[INDS]
	templon = lonK[INDS]
	templat = latK[INDS]
	if len(templat) <= 1:
		kPM[ti,:,:] = np.nan
	else:
		ok = OrdinaryKriging( templon, templat, datatemp, variogram_model='spherical', \
			variogram_parameters = vgp, verbose=False, enable_plotting=False)
		z, ss = ok.execute( 'points', glons.reshape((1,-1)), glats.reshape((1,-1)) )
		kPM[ti, :, :] = z.reshape( (NX, NY) )
# print NED, EDind
# check to make sure there are a decent number of sites with data for each day, NED = 0

# save output
fname = '/home/kaodell/nasa_fires/Sheena_Den3kmKrige/updated/kpmDEN'+ str(year) + '_3km_update.npz'
f = open(fname, 'w')
np.savez(f, alldays=alldays, kPM=kPM,glats=glats,glons=glons,kRSQ=kRSQ,sPMin=sPMin,\
	sPMout=sPMout,ilon=ilon,ilat=ilat,olat=olat,olon=olon,slope=slope,kPM_site=kPM_site,\
	intercept=intercept,rvalue=rvalue,pvalue=pvalue,stderr=stderr,rsq=rsq,MB=MB,MAE=MAE)
f.close()


'''
############################################################################################################
# Figures
############################################################################################################

f = np.load('/home/kaodell/nasa_fires/Martenies_Den3kmKrige/processed_datafiles/kpmDEN2012_3km.npz') 
sPMin = f['sPMin']  
glons = f['glons']
glats = f['glats']
kPM = f['kPM']
ASites = f['ASites']
ilon = f['ilon']
ilat = f['ilat']
f.close()

############################################################################################################
# FIGURE 1: Map of WRF grid and krige sites 

pl.figure()
m = Basemap(projection='merc',lon_0=-129,llcrnrlat=35,urcrnrlat=44,\
    llcrnrlon=-112,urcrnrlon=-96,lat_ts=10,resolution='l')
# draw coastlines, state and country boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()  
m.shadedrelief()
isx, isy = m(ilon, ilat)
osx, osy = m(olon,olat)
kx, ky = m(glons, glats)
cs = m.pcolor(kx,ky,np.zeros([glons.shape[0],glats.shape[1]]),color = 'grey',alpha=0.1)
m.scatter(isx,isy, color = 'r',s=4)
m.scatter(osx,osy, color = 'k',s=4)
pl.title('Colorado kriging grid and surface sites')
pl.savefig('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/figures/3kmgrid_w_krigesites_COmtncutoff.png')
'''



