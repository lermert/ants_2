# determine kernel values for ray-theoretical kernels, bin measurements, save the maps, and plot them
import numpy as np
from math import pi, isnan
import pandas as pd
import os
from obspy import UTCDateTime
from ants_2.tools.geo import get_midpoint, get_antipode, area_of_sqdeg
from geographiclib import geodesic, geodesicline
from ants_2.tools.plot import plot_grid
from glob import glob

class sourcemap_2(object):

	def __init__(self,csv_file,kernel_dir,min_snr=0.,
		t0=None,prefix=None,msr='obs',syn='syn',cha='MXZ'):

		self.data = pd.read_csv(csv_file)
		self.kernel_dir = kernel_dir
		self.min_snr = min_snr
		
		self.t0 = t0
		self.prefix = prefix if prefix is not None else './'
		self.msr = msr
		self.syn = syn
		self.cha = cha

	def assemble_descent(self):
		# loop over stationpairs
		cnt_success = 0
		cnt_lowsnr = 0
		cnt_lown = 0
		cnt_overlap = 0
		cnt_unavail = 0

		if self.t0 is not None:
			self.data = self.data[data.t0 == self.t0]
		n = len(self.data)


		for i in range(n):
			print i
			if self.data.at[i,'snr_c'] < self.min_snr\
			 and self.data.at[i,'snr_a'] < self.min_snr:
				cnt_lowsnr += 1
				continue
			# ToDo: deal with station pairs with several measurements (with different instruments)
			# (At the moment, just all added. Probably fine on this large scale)
			# find kernel file
			sta1 = self.data.at[i,'sta1']
			sta2 = self.data.at[i,'sta2']
			print sta1,sta2
			if sta1.split('.')[-1][-1] in ['E','N','T','R']:
				msg = "Cannot yet handle horizontal components"
				raise NotImplementedError(msg)
			if sta2.split('.')[-1][-1] in ['E','N','T','R']:
				msg = "Cannot yet handle horizontal components"
				raise NotImplementedError(msg)
		
		
			# ToDo !!! Replace this by a decent formulation, where the channel is properly set !!! No error for E, R, T, N
			sta1 = "*.{}..{}".format(sta1.split('.')[1],self.cha) # ignoring network: IRIS has sometimes several network codes at same station
			sta2 = "*.{}..{}".format(sta2.split('.')[1],self.cha) # ignoring network: IRIS has sometimes several network codes at same station
		
			kernelfile1 = os.path.join(self.kernel_dir,"{}--{}.npy".\
				format(sta1,sta2))
			kernelfile2 = os.path.join(self.kernel_dir,"{}--{}.npy".\
				format(sta2,sta1))
			# Same problem with different network codes.
			# Due to station pairs being in alphabetic order of network.station.loc.cha, different network
			# codes also lead to different ordering.
			try:
				kernelfile = glob(kernelfile1)[0]
			except IndexError:
				try: 
					kernelfile = glob(kernelfile2)[0]
				except IndexError:
					kernelfile = kernelfile1 # this file does not actually exist, but there might be a good reason.
					# Check that first, and then complain.
			print kernelfile
			# Skip if entry is nan: This is most likely due to no measurement taken because station distance too short	
			if isnan(self.data.at[i,self.msr]):
				print("No measurement in dataset for:")
				print(sta1)
				print(sta2)
				cnt_overlap += 1
				continue

			# ...unless somehow the kernel went missing (undesirable case!)

			if not os.path.exists(kernelfile):
				print("File does not exist:")
				print(os.path.basename(kernelfile))
				cnt_unavail += 1
				continue


			# load kernel
			kernel = np.load(kernelfile)
			# ToDo Ugly!
			if 'gradient' not in locals():
				gradient = np.zeros(kernel.shape)
			if True in np.isnan(kernel):
				print("kernel contains nan, skipping")
				print(os.path.basename(kernelfile))
				continue


			# if everythin worked: 
			

			# Find synthetic measurem, multiply kernel and measurement, add to descent dir. 
			
			if self.msr == 'obs':
				kernel *= (self.data.at[i,'syn'] - self.data.at[i,self.msr])
			elif msr == 'enr_a':
				kernel *= 2. * (self.data.at[i,'syn_a'] - self.data.at[i,self.msr])
			elif msr == 'enr_c':
				kernel *= 2. * (self.data.at[i,'syn_c'] - self.data.at[i,self.msr])
			cnt_success += 1 # yuhu

			gradient += kernel
			
			del kernel


		# Save the positive kernel
		kernelfile = self.prefix+'grad_all.npy'
		np.save(kernelfile,gradient)




class sourcemap(object):

	def __init__(self,csvfile,v,f,q,seg_km,min_snr=0.,min_win=1,
		t0=None,prefix=None):

		self.v = v
		self.f = f
		self.q = q
		self.w = 2 * pi * f
		self.seg_km = seg_km
		self.min_snr = min_snr
		self.min_win = min_win
		self.t0 = t0
		self.prefix = prefix if prefix is not None else './'
		self.data = pd.read_csv(csvfile)
		print self.data.keys()


	def plot_sourcemap(self):

		if os.path.exists(self.prefix+'.sourcemap.npy'):

			smap = np.load(self.prefix+'.sourcemap.npy')

			plot_grid(smap[0],smap[1],smap[2],outfile=self.prefix+'.sourcemap.png',
				normalize=False)
			plot_grid(smap[0],smap[1],smap[3],outfile=self.prefix+'.raycntmap.png',
				cmap='seq',vmin=0,vmax=1.0,normalize=True)


	def _bin_kernels(self,ddeg_lon,ddeg_lat,lonmin=110,lonmax=160,latmin=15,latmax=60):

		print("Attention, Japan setting hardcoded...")

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

		np.save(self.prefix+'.sourcemap.npy',smap)
		os.system('rm tempfile.txt')
		

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
			snr_a  = self.data.at[i,'snr_a']
			snr_c  = self.data.at[i,'snr_c']
			n_win = self.data.at[i,'nstk']
			starttime = self.data.at[i,'t0']

			if snr_a < self.min_snr and snr_c < self.min_snr:
				continue

			if n_win < self.min_win:
				continue

			if (self.t0 is not None 
			and UTCDateTime(self.t0) != UTCDateTime(starttime)):
				continue

			if self.data.at[i,'obs'] == np.nan:
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
				
			
			# Multiply by measurement
			kernel[:,2] *= self.data.at[i,'obs'] * -1. * -1. #The first -1 is because 
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
        
    	

