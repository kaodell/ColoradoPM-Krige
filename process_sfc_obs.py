#process_sfc_obs.py
#	pull in csv files of PM2.5 data from the AQS and CDPHE and save as array of daily PM for each site
#		this is the first part of the krige_aqs_US code, with processing of CDPHE data added
#written by: Kate O'Dell 
#version history: 1.0 initial 09/30/16 (version in krige code)
#		  2.0 this version 05/11/18
###################################################################################################
# load modules
import datetime as dt
import numpy as np
import csv
import sys
###################################################################################################
# set year
#year = 2012

#indicate include CDPHE data
CDPHE = True

# to run multiple years at once

if len(sys.argv) > 1:
	print sys.argv[1]
	#sys.exit()
	year = int(sys.argv[1])
else:
	print 'no'

###################################################################################################
#Time period of interest
###################################################################################################
t0 = dt.datetime( year = year, month = 1, day = 1 )
tf = dt.datetime( year = year, month = 12, day = 31 )
NT = (tf - t0).days
alldays = np.empty( NT+1, dtype = object )
alldays[0] = t0
for i in range(1, NT+1):
	alldays[i] = alldays[i-1] + dt.timedelta(days = 1)

###################################################################################################
#Part 1a: Get AQS Data
###################################################################################################

#First EPA AQS files
fn1 = '/fischer-scratch/kaodell/surface_obs/01-24-18dwnload/daily_88101_' + str(year) +'.csv'
fn2 = '/fischer-scratch/kaodell/surface_obs/01-24-18dwnload/daily_88502_' + str(year) + '.csv'
fns = [fn1,fn2]
#fid = np.loadtxt( fn, dtype=str, delimiter='","' )

lats = []
lons = []
dates = []
pmraw = []
ID = []
data_type = []

ulons = []
ulats = []
Asites = []

for i in range(2):
	fid = open(fns[i], 'r')

	ro = csv.reader( fid, delimiter = ',' )
	count = 0
	N = ro.__sizeof__()

	for r in ro:
		#print count
		if count==0:
			header = r
			count += 1
			continue
		#only use lines where there is no flag	
		if r[13] == 'None':	
			IDtemp = 0
			lats.append(r[5])
			lons.append(r[6])
			dstr = r[11]
			dtemp = dt.datetime( year = int(dstr.split('-')[0]), month = int(dstr.split('-')[1]), \
				day = int( dstr.split('-')[2]) )
			dates.append(dtemp)
			pmraw.append(r[16])
			if r[0] == 'CC':
				#these are canada sites, just change state code to number so works with id
				r[0] = '99'
			IDtemp = (i+1)*10**18 + int(r[0])*10**16 + int(r[1])*10**11 + int(r[2])*10**7 + int(r[4])*10**2
			#print int(r[0])
			#i indicates which file it came from, site numbers recycle between files			
			#print IDtemp
			if (int(r[14])) > 1:	
				dtype = 24
				data_type.append(dtype)
			elif (int(r[14])) == 1:
				dtype = 1
				data_type.append(dtype)
			IDtemp = IDtemp + dtype
			ID.append(IDtemp)
		count += 1

	fid.close()

dates = np.array(dates, dtype = object)
pmraw = np.array(pmraw)
lats = np.array(lats)
lons = np.array(lons)
#ID = np.array( ID )
ASites = np.unique(ID)
data_type = np.array(data_type)


#Get unique lats/lons for each site
ulons = []
ulats = []
for cid in ASites:
	cIND = np.where( ID == cid )[0][0]
	ulons.append( lons[cIND] )
	ulats.append( lats[cIND] )

lons = np.array( ulons ).astype(float)
lats = np.array( ulats ).astype(float)

###################################################################################################
#Part 1b: Get CDPHE Data (avail 2009 - 2016) (if CDPHE=True)
###################################################################################################
if CDPHE:

	fn3 = '/fischer-scratch/kaodell/surface_obs/cdphe_sites/pm2.5_colorado2009-2016_lyfix.csv'
	fid = open(fn3, 'rU')
	ro = csv.reader( fid )
	count = 0
	COpm25 = []
	for r in ro:
		if count==0:
			header = r
			count += 1
			continue
		#pull only this years data	
		if r[0][-2:] == str(year)[-2:]:
			COpm25.append(r[2:18])

	COpm25 = np.array(COpm25)

	naninds1 = np.where(COpm25 == '-')
	COpm25[naninds1] = np.nan

	naninds2 = np.where(COpm25 == '')
	COpm25[naninds2] = np.nan

	COpm25 = np.array(COpm25,dtype = 'float')

	dailyCOpm25 = np.empty([len(alldays),16])
	dailyCOpm25[:] = np.nan
	dind = 0
	for i in range(0,24*len(alldays),24):
		dayPM = COpm25[i:i+24,:]
		for j in range(16):
			nmeas = len(np.where(np.isnan(dayPM[:,j]) == False )[0])
			if nmeas >= 19:
				dailyCOpm25[dind,j] = np.nanmean(dayPM[:,j])
		dind += 1
	fid.close()

	# load site lat lons
	clons = []
	clats = []
	fn4 = '/fischer-scratch/kaodell/surface_obs/cdphe_sites/pm2.5_colorado2010-2014_siteloc.csv'
	fid = open(fn4, 'rU')
	ro = csv.reader(fid)
	for r in ro:
		clons.append(r[1])
		clats.append(r[2])
	
	# make numpy arrays and remove header
	clons = np.array(clons[1:],dtype = 'float')
	clats = np.array(clats[1:],dtype = 'float')

####################################################################################################
#Part 2: Create sPM matrix
#####################################################################################################
NT = len(alldays)
NS = len(ASites)

sPM = np.empty( [NT, NS] )
errs = []
errs2 = []
check = 0
j = 0
nnan = 0
for csID in ASites:
	sINDS = np.where(ID == csID)[0]
	#if len(sINDS[0]) == NT:
	#	sPM[:,j] = pmraw[sINDS]
	#else:
	i = 0		
	for i in range(0,len(alldays)):
		dIND = np.where( dates[sINDS] == alldays[i])
		if len(dIND[0]) == 1:
			sPM[i,j] = pmraw[sINDS[dIND[0]][0]]
			check += 1
		elif len(dIND[0]) == 2:
			#dt1 = data_type[sINDS[dIND[0][0]]]
			#dt2 = data_type[sINDS[dIND[0][1]]]
			#if dt1 != dt2:
			#	sPM[i,j] = pmraw[sINDS[dIND[0][0]]]				
			#else:
			errs.append([i,j])
			
		else:
			sPM[i,j] = np.nan
			nnan += 1
	j+=1


sPM = np.array(sPM)
lons = np.array(lons)
lats = np.array(lats)
ID = np.array(ID)
ASites = np.array(ASites)


if CDPHE:
	# append CDPHE estimates and site locations
	sPM = np.hstack([sPM,dailyCOpm25])
	lons = np.hstack([lons,clons])
	lats = np.hstack([lats,clats])

#find and delete sites that have no data during the specified time frame
#print nnan
naninds = np.where(np.isnan(sPM)==True)

count = []
for i in range(sPM.shape[1]):
	count.append(len(np.where(naninds[1] == i)[0]))
count = np.array(count)
inds = np.where(count == NT)	# these are sites that have no data during the specified time frame

# delete these sites from the data
sPM = np.delete(sPM,inds,axis=1)
lons = np.delete(lons,inds)
lats = np.delete(lats, inds)
#ASites = np.delete(ASites, inds)

print len(naninds[0])
print inds
print sPM.shape[0]*sPM.shape[1] - len(pmraw)

#f = open('/home/kaodell/nasa_fires/processed_datafiles/US_multiyear_krige/sfc_pm' + str(year)[2:] + '.npz', 'w')
#np.savez(f, sPM=sPM, lons=lons, lats=lats, ID=ID, ASites=ASites, errs=errs,nnan=nnan, pmraw=pmraw,NT=NT,NS=NS)
#f.close()

f = open('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/processed_datafiles/update/sfc_pm' + str(year)[2:] + '_wCDPHE.npz', 'w')
np.savez(f, sPM=sPM, lons=lons, lats=lats, errs=errs,nnan=nnan, pmraw=pmraw)
f.close()

####################################################################################################
#optional: double check we don't have duplicate sites from CDPHE
#####################################################################################################

#for i in range(sPM.shape[0])



		






