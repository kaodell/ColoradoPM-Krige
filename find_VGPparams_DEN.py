#find_VGPparams_DEN.py
#	read in output file from batch testing of different vgp parameters and find the ones that
#	maximize r2, minimize mean error, and mean bias
#Written by : Katelyn O'Dell 3/14/17 (Pi Day!)

# updates 03/30/18 to work with differnt output file, renamed find_VGPparams_Den.py
###################################################################################################
# Import Modules
###################################################################################################
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
###################################################################################################
# Read in File
###################################################################################################

# open file
fid = open('vgpparams14.log','r') 

s = []
r = []
n = []
rsq = []
slope = []
MB = []
MAE = []

# When reading files this way, Python reads files line-by-line
ind = -1
count = 0
for line in fid.readlines(): # this reads one line at a time
    if count == 0:
        spl_line = line.split() # split into a list of strings by whitespace
        s.append(round(float(spl_line[0]),1)) 
        r.append(round(float(spl_line[1]),1)) 
	n.append(round(float(spl_line[2]),1))
	count += 1
	ind += 1
    elif count == 1:				# this is rsq
	spl_line = line.split()
	rsq.append( float(spl_line[0]))
	slope.append( float(spl_line[1]))
	MB.append( float(spl_line[2]))
	MAE.append( float(spl_line[3]))
	count = 0

s = np.array(s) 
r = np.array(r) 
n = np.array(n) 
rsq = np.array(rsq) 
slope = np.array(slope) 
MB = np.array(MB) 
MAE = np.array(MAE) 
'''
#find mean and std deviation of each set
rsq_m = np.mean(rsq)
slope_m = np.mean(slope)
MB_m = np.mean(MB)
MAE_m = np.mean(MAE)
'''
rsq_sd = np.std(rsq)
slope_sd = np.std(slope)
MB_sd = np.std(MB)
MAE_sd = np.std(MAE)


#to narrow down search for parameters,find top 10% of all stats, and look at parameters that fit all stats
su = np.unique(s)
ru = np.unique(r)
nu = np.unique(n)


rsq_mat = np.empty([len(su), len(ru), len(nu)])
MB_mat = np.empty([len(su), len(ru), len(nu)])
MAE_mat = np.empty([len(su), len(ru), len(nu)])
slope_mat = np.empty([len(su), len(ru), len(nu)])

for i in range(len(s)):
	sind = np.where(s[i] == su)
	rind =  np.where(r[i] == ru)
	nind = np.where(n[i] == nu)
	rsq_mat[sind,rind,nind] = rsq[i]
	MB_mat[sind,rind,nind] = MB[i]
	MAE_mat[sind,rind,nind] = MAE[i]
	slope_mat[sind,rind,nind] = slope[i]

rsq_sort = np.sort(-rsq)
rsq_sort = rsq_sort*-1.
rsq_top10 = rsq_sort[0:int(.2*len(rsq))]

MB_sort = np.sort(abs(MB))
MB_top10 = MB_sort[0:int(.2*len(rsq))]

MAE_sort = np.sort(MAE)
MAE_top10 = MAE_sort[0:int(.2*len(rsq))]

slope_sort = np.sort(abs(1.-slope))
slope_top10 = slope_sort[0:int(.2*len(rsq))]

MB_ind = []
MAE_ind = []
rsq_ind = []
slope_ind = []

for i in range(len(rsq_top10)):
	ind = np.where(abs(MB) == MB_top10[i])[0]
	for j in range(len(ind)):
		MB_ind.append(ind[j])
	ind = np.where(MAE == MAE_top10[i])[0]
	for j in range(len(ind)):
		MAE_ind.append(ind[j])
	ind = np.where(rsq == rsq_top10[i])[0]
	for j in range(len(ind)):
		rsq_ind.append(ind[j])
	ind = np.where(abs(1.-slope) == slope_top10[i])[0]
	for j in range(len(ind)):
		slope_ind.append(ind[j])

MB_ind = np.array(MB_ind)
MAE_ind = np.array(MAE_ind)
rsq_ind = np.array(rsq_ind)
slope_ind = np.array(slope_ind)

top10params = []
for i in range(len(MB_ind)):
	if MB_ind[i] in MAE_ind and MB_ind[i] in rsq_ind and MB_ind[i] in slope_ind:
		top10params.append(MB_ind[i])


s_top10 = s[top10params]
r_top10 = r[top10params]
n_top10 = n[top10params]
vgptop10 = [s_top10,r_top10,n_top10]





inds = np.where(r_top10 == 2.5)




