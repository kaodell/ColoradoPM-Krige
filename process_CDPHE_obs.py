#process_CDPHE_obs.py
#	pull in csv files of PM2.5 data from the CDPHE and save as array of daily PM for each site
#		this is the first part of the krige_aqs_US code, with processing of CDPHE data added
#written by: Kate O'Dell 
#version history: 1.0 initial 09/30/16 (version in krige code)
#		  2.0 cdphe + AQS version 05/11/18
#		  3.0 only CDPHE 06/05/18
###################################################################################################
# load modules
import datetime as dt
import numpy as np
import csv
import sys
###################################################################################################
# set year
year = 2014

#indicate include CDPHE data
CDPHE = True
'''
# to run multiple years at once

if len(sys.argv) > 1:
	print sys.argv[1]
	#sys.exit()
	year = int(sys.argv[1])
else:
	print 'no'
'''
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
#Part 1b: Get CDPHE Data (avail 2009 - 2016) 
###################################################################################################


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
		COpm25.append(r[2:19])

COpm25 = np.array(COpm25)

naninds1 = np.where(COpm25 == '-')
COpm25[naninds1] = np.nan

naninds2 = np.where(COpm25 == '')
COpm25[naninds2] = np.nan

COpm25 = np.array(COpm25,dtype = 'float')

dailyCOpm25 = np.empty([len(alldays),17])
dailyCOpm25[:] = np.nan
dind = 0
for i in range(0,24*len(alldays),24):
	dayPM = COpm25[i:i+24,:]
	for j in range(17):
		nmeas = len(np.where(np.isnan(dayPM[:,j]) == False )[0])
		if nmeas >= 19:
			dailyCOpm25[dind,j] = np.nanmean(dayPM[:,j])
	dind += 1
fid.close()

# load site lat lons
clons = []
clats = []
fn4 = '/fischer-scratch/kaodell/surface_obs/cdphe_sites/pm2.5_colorado2009-2016_lyfix_locs.csv'
fid = open(fn4, 'rU')
ro = csv.reader(fid)
for r in ro:
	clons.append(r[1])
	clats.append(r[2])
clons = np.array(clons)
rinds = np.where(clons =='')	
clons = np.delete(clons,rinds)
clats = np.delete(clats,rinds)
# make numpy arrays and remove header
clons = np.array(clons[1:],dtype = 'float')
clats = np.array(clats[1:],dtype = 'float')

sPM = dailyCOpm25
lons = clons
lats = clats

####################################################################################################
# save data
#####################################################################################################


f = open('/home/kaodell/nasa_fires/Sheena_Den3kmKrige/processed_datafiles/update/sfc_pm' + str(year)[2:] + '_onlyCDPHE.npz', 'w')
np.savez(f, sPM=sPM, lons=lons, lats=lats)
f.close()
