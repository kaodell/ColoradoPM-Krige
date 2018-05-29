#npz2netCDF.py
#	A python script to load npz data and convert to a netCDF data format
#written by: Kate O'Dell
#version history: 1.0 initial 04/13/18
###################################################################################################
# import modules
###################################################################################################
import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime as dt
###################################################################################################
# load data and make datetime array
###################################################################################################
year = 2013

t0 = dt.datetime( year = year, month = 1, day = 1 )
tf = dt.datetime( year = year, month = 12, day = 31 )
NT = (tf - t0).days
alldays = np.empty( NT+1, dtype = object )
alldays[0] = t0
date = ['na']*( NT+1)
for i in range(1, NT+1):
	alldays[i] = alldays[i-1] + dt.timedelta(days = 1)
	date[i] = str(alldays[i])[:10]
fname =  '/home/kaodell/nasa_fires/Sheena_Den3kmKrige/processed_datafiles/kpmDEN'+ str(year) + '_3km.npz'
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
# save as netcdf
###################################################################################################

lonx = [0]*glons.shape[0]
lony = [0]*glons.shape[1]

outfn = '/home/kaodell/nasa_fires/Sheena_Den3kmKrige/processed_datafiles/DENkrigedPM25_' + str(year) + '.nc'
nc_w_fid = netCDF4.Dataset(outfn, 'w', format='NETCDF4')
nc_w_fid.description = 'Krigged (interpolated) PM2.5 concentrations from EPA AQS surface monitors for Jan-Dec' + str(year)

# Need to define dimensions that will be used in the file
nc_w_fid.createDimension('date', None) # If "None", then time can be appended
nc_w_fid.createDimension('lonx')
nc_w_fid.createDimension('lony')
nc_w_fid.createDimension('kfold')
#nc_w_fid.createDimension('glats', None)
# define variables for netcdf file
                          #      name  format  dimensions
date_nc = nc_w_fid.createVariable('date', 'S1', ('date',)) # making the name the same as the dimention
                                                           # makes time a dependent variable
lon_nc = nc_w_fid.createVariable('lon', 'f8', ('lonx','lony'))

lat_nc = nc_w_fid.createVariable('lat', 'f8', ('lonx','lony'))

PM25_nc = nc_w_fid.createVariable('PM25', 'f8', ('date', 'lonx','lony'))

rsq_nc = nc_w_fid.createVariable('r-squared', 'f8', ('kfold'))

MB_nc = nc_w_fid.createVariable('mean bias', 'f8', ('kfold'))

MAE_nc = nc_w_fid.createVariable('mean absolute error', 'f8', ('kfold'))


date_nc.setncatts({'units':'date','long_name':'Date time object for day of 24-hr avg PM2.5conc',\
               'var_desc':'Datetime Object [date]'})
lat_nc.setncatts({'units':'degrees','long_name':'degrees latitude for data grid',\
               'var_desc':'Latitude [degrees]'})
lon_nc.setncatts({'units':'degrees','long_name':'degrees longitude for data grid',\
               'var_desc':'Longitude [degrees]'})
PM25_nc.setncatts({'units':'ug/m3','long_name':'24 hour average PM2.5 concentration',\
               'var_desc':'PM2.5 concentration [ug/m3]'})
rsq_nc.setncatts({'units':'','long_name':'r-squared',\
               'var_desc':'r-squared for each k-fold cross validation against sfc monitors'})
MB_nc.setncatts({'units':'binary','long_name':'mean bias',\
               'var_desc':'mean bias for each k-fold cross validation against sfc monitors'})
MAE_nc.setncatts({'units':'binary','long_name':'mean absolute error',\
               'var_desc':'mean absolute error for each k-fold cross validation against sfc monitors'})

# data (use date read from the csv file)
date_nc[:] = date
lon_nc[:] = glons
lat_nc[:] = glats
PM25_nc[:] = kPM
rsq_nc[:] = rsq
MB_nc[:] = MB
MAE_nc[:] = MAE
nc_w_fid.close()




