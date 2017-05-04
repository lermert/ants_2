# determine kernel values for ray-theoretical kernels, bin measurements, save the maps, and plot them
import numpy as np
from math import pi
import pandas as pd
import os
from ants_2.tools.geo import get_midpoint, get_antipode, area_of_sqdeg
from geographiclib import geodesic, geodesicline
from ants_2.tools.plot import plot_grid
from warnings import warn

class sourcemap(object):

	def __init__(self,csvfile,v,f,q,seg_km,min_snr=0.):

		self.v = v
		self.f = f
		self.q = q
		self.w = 2 * pi * f
		self.seg_km = seg_km
		self.min_snr = min_snr

		self.data = pd.read_csv(csvfile)


	def plot_sourcemap(self):

		if os.path.exists('sourcemap.npy'):

			smap = np.load('sourcemap.npy')

			plot_grid(smap[0],smap[1],smap[2],outfile='sourcemap.png')
			plot_grid(smap[0],smap[1],smap[3],outfile='raycntmap.png',
				cmap='seq',vmin=0,vmax=1.0)


	def _bin_kernels(self,ddeg_lon,ddeg_lat,lonmin=-180,lonmax=180,latmin=-90,latmax=90):

		lats = np.arange(latmin,latmax+ddeg_lat,ddeg_lat)
		lons = np.arange(lonmin,lonmax+ddeg_lon,ddeg_lon)
		vals = np.zeros((len(lons),len(lats)))
		hits = np.zeros((len(lons),len(lats)))


		data = open('tempfile.txt','r').read().split('\n')

		for entry in data:

			if entry.split() == []: continue

			lat = float(entry.split()[0])
			lon = float(entry.split()[1])
			val = float(entry.split()[2])

			if lon > lonmax: continue
			if lat > latmax: continue
			if lon < lonmin: continue
			if lat < latmin: continue

			# Index 1 - longitude index
			i1 = int(round((lon-lonmin)/ddeg_lon))
			if i1 > len(lons)-1: continue

			# Index 2 - latitude index
			i2 = int(round((lat-latmin)/ddeg_lat))
			if i2 > len(lats)-1: continue

			vals[i1,i2] += val
			hits[i1,i2] += 1

		# Write into an numpy file in x,y,z,h format:
		weight0 = area_of_sqdeg(0.)
		x = []
		y = []
		z = []
		h = [] 
		for i in range(len(lons)):
			for j in range(len(lats)):


				# weight bins by their surface area and hit count
				if hits[i,j] == 0:
					weight = 0
				else:
					weight = area_of_sqdeg(lats[j]) / weight0 / hits[i,j]

				x.append(lons[i])
				y.append(lats[j])
				z.append(vals[i,j] * weight)
				h.append(hits[i,j])

		smap = np.array(zip(x,y,z,h)).transpose()
		np.save('sourcemap.npy',smap)
		#os.system('rm tempfile.txt')
		

	def _temp_kernels(self):

		# Save all the ray kernel coordinates and values to a temporary file
		fh = open('tempfile.txt','w')

		n = len(self.data)

		gf_exp = self.w / (self.v * self.q)
		
		for i in range(n):
			lat1 = self.data.at[i,'lat1']
			lon1 = self.data.at[i,'lon1']
			lat2 = self.data.at[i,'lat2']
			lon2 = self.data.at[i,'lon2']
			snr  = self.data.at[i,'snr']

			if snr < self.min_snr:
				continue

			# determine antipode of midpoint
			# find midpoint
			mp = get_midpoint(lat1,lon1,lat2,lon2)
			# find antipode of midpoint
			ap = get_antipode(mp[0],mp[1])

			# distance to antipode of midpoint and number of segments on discrete ray:
			dist = geodesic.Geodesic.WGS84.Inverse(ap[0],ap[1],lat1,lon1)['s12']
			num_seg = int( dist / 1000. / self.seg_km )

			# kernel array: lat, lon, value
			kernel = np.zeros((2*num_seg,3))

        	# Fill two segments: station1 to antipode of midpoint, station 2 to antipode of midpoint
			for nr in [1,2]:

				# Discrete ray
				lat = self.data.at[i,'lat'+str(nr)]
				lon = self.data.at[i,'lon'+str(nr)]

				
				i0 = (nr-1) * num_seg
				
				
				kernel[i0:i0+num_seg,:] = self.discrete_ray(lat,lon,ap[0],ap[1],num_seg)

				#ToDo ask Andreas about distance measure
				
				# determine kernel
				if nr == 1:
					kernel[i0:i0+num_seg,2] = np.exp(-gf_exp * kernel[i0:i0+num_seg,2]) 

				elif nr == 2:
					kernel[i0:i0+num_seg,2] = - np.exp(-gf_exp * kernel[i0:i0+num_seg,2]) 
				
			msr = self.data.at[i,'obs']
			if np.isnan(msr):
				warn("NaN measurements found.")
				continue 

			# Multiply by measurement
			kernel[:,2] *= msr * -1. * -1. #The first -1 is because 
			# the definition of the kernel is (A-A_obs) * K_dataless, and here A=0, so we
			# need just -A_obs.
			# The second -1 is because we want to look at the negative of the misfit gradient
			# (we want to decrease misfit by our update; so that is the update direction)
			
			# append ray coordinates and kernel to the temporary file
			for k in range(2*num_seg):
				fh.write("%7.2f %7.2f %8.5f\n" %tuple(kernel[k,:]))


		fh.close()







	def discrete_ray(self,lat1,lon1,lat2,lon2,num):
	    """
	    Obtain segment coordinates that separate the great circle b/w lat1,lon1 
	    and lat2,lon2 into num equal segments.
	    index: one row of the data frame containing coordinates and measurement.
	    
	    """


	    # Get a geodesic line object
	    p=geodesic.Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2) 
	    l=geodesic.Geodesic.WGS84.Line(p['lat1'],p['lon1'],p['azi1'])
	    seg_len = p['s12']/(num)

	    ray = np.zeros((num,3))

	    
	    
	    
	    for i in range(num):

	        b = l.Position( i * seg_len )
	        
	        
	        ray[i,0] = b['lat2']
	        ray[i,1] = b['lon2']
	        ray[i,2] = i * seg_len
	    return ray

	            
	    #else: 
	    #   #determine omega
	    #   w = 2 * pi * freq
	    #   for i in range(num):
	    #       #determine x; the path obtained with geodesic is split in equal length segments.
	    #       x = sta_dist / 2. + i * p['s12']/num
	    #       #x = i * p['s12'] / (num)
	#
	    #       # determine K on the station-station line:
	    #       if line_kern == True:
	    #           k = exp(-w * x / (v * Q))/sqrt(x**2+sta_dist**2/4)
	    #       # The following kernel is integrated along the y-direction
	    #       else:
	    #           k = exp(-w * x / (v * Q))
	    #       b=l.Position(i*p['s12']/(num))
	   #        newcoords.append((b['lat2'],b['lon2'],k)) 
        
    	

