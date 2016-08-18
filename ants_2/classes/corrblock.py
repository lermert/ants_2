# Correlation block object:
from obspy import Stream, read_inventory, UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth

import numpy as np
import os
import re


from ants_2.tools.preprocess import ram_norm, whiten
from ants_2.tools.correlations import cross_covar
# list of possible channels combinations indicating that the data needs to be rotated.
horizontals = ['RR','RT','TR','TT','TZ','ZT','RZ','ZR']


def get_geoinf(id1,id2):


		inv1 = '{}.{}.xml'.format(*id1.split('.')[0:2])
		inv2 = '{}.{}.xml'.format(*id2.split('.')[0:2])

		inv1 = read_inventory(os.path.join('meta','stationxml',inv1))
		inv2 = read_inventory(os.path.join('meta','stationxml',inv2))

		# Replace 'radial' and 'transverse' by 'N' and 'E'
		id1 = re.sub('\.??R$','N',id1)
		id2 = re.sub('\.??R$','N',id2)
		id1 = re.sub('\.??T$','E',id1)
		id2 = re.sub('\.??T$','E',id2)
		

		c1 = inv1.get_coordinates(id1)
		c2 = inv2.get_coordinates(id2)

		lat1, lon1, lat2, lon2 = (
			c1['latitude'],
			c1['longitude'],
			c2['latitude'],
			c2['longitude'])

		dist, az, baz = gps2dist_azimuth(lat1,lon1,lat2,lon2)

		return lat1, lon1, lat2, lon2, dist, az, baz



class CorrTrace(object):

	"""
	Object holds correlation data along with metainformation (station id, geographic location).
	"""

	def __init__(self,cha1,cha2,corr_type,nlag,t0=None,t1=None):

		self.stack = np.zeros(nlag)
		self.pstak = np.zeros(nlag)

		self.begin = t0
		self.end   = t1
		self.id1   = cha1
		self.id2   = cha2
		self.id    = cha1 + '--' + cha2
		self.ctype = corr_type


		self.cnt_tot = 0
		self.cnt_int = 0

		geo_inf = get_geoinf(cha1,cha2)
		
		self.lat1 = geo_inf[0]
		self.lat2 = geo_inf[2]
		self.lon1 = geo_inf[1]
		self.lon2 = geo_inf[3]
		self.az   = geo_inf[5]
		self.baz  = geo_inf[6]
		self.dist = geo_inf[4]


		# open the file to dump intermediate stack results
		int_file = '{}.{}.windows.bin'.format(self.id,self.ctype)
		int_file = os.path.join('data','correlations',int_file)
		self.int_file = open(int_file,'wb')

	


	def _add_corr(self):
		pass

	def write_tot(self):
		pass

	def write_int(self):
		pass

	

class CorrBlock(object):


# - initialize with station pairs
	def __init__(self,block,cfg):

		
		self.correlations = []
		self.inv = block.inventory
		self.channels = block.channels
		self.station_pairs = block.station_pairs
		self.channel_pairs = block.channel_pairs


		self.initialize_data()
		self.sampling_rate = self.data[0].stats.sampling_rate
		self.delta = self.data[0].stats.delta

		n_lag = 13

		for cp in block.channel_pairs:
			for pair in cp:
				try:
					self.correlations.append(CorrTrace(pair[0],pair[1],
					cfg.corr_type,n_lag))
				except:
					print('** Could not initialize correlation for %s,%s: check metadata'
						%(pair[0],pair[1]))


		
		if any(i in cfg.corr_tensorcomponents for i in horizontals):

			self.azms = []
			self.bazs = []

			for pair in self.station_pairs:
				try:
					geoinf = get_geoinf(pair[0],pair[1])
					# - find azimuth, backazimuth
					self.azms.append(geoinf[5])
					self.bazs.append(geoinf[6])
				except:
					self.azms.append(0)
					self.bazs.append(0)


		

	def run(self,cfg):

		t_0 = UTCDateTime(cfg.time_begin)
		t_end = UTCDateTime(cfg.time_end)
		win_len_seconds = cfg.time_window_length
		win_len_samples = int(round(win_len_seconds*self.sampling_rate))
		min_len_samples = int(round(cfg.time_min_window*self.sampling_rate))
		max_lag_samples = int(round(cfg.corr_maxlag * self.sampling_rate))



		
		t = t_0
		

		while t <= t_end - (win_len_seconds - self.delta):
			print(t)
			

			
			# - check endtime, if necessary, add data from 'later' file
			self.update_data(t, win_len_seconds)

			# - slice the traces
			windows = self.data.slice(t, t + win_len_seconds - self.delta)
			

			# self.preprocess?
			# - Apply preprocessing
			for w in windows:
				# - check minimum length requirement

				if cfg.whiten:
					w = whiten(w,cfg.white_freqmin,cfg.white_freqmax,cfg.white_taper)

				if cfg.onebit:
					w.data = np.sign(w.data)

				if cfg.ram_norm:
					w = ram_norm(w,cfg.ram_window,cfg.ram_prefilt)
				

			# - station pair loop
			for sp_i in range(len(self.station_pairs)):

				pair = self.station_pairs[sp_i]
				print(pair)
				# - select traces
				[net1, sta1] = pair[0].split('.')
				[net2, sta2] = pair[1].split('.')
				
				str1 = windows.select(network=net1, station=sta1)
				str2 = windows.select(network=net2, station=sta2)
				

				# - if horizontal components are involved, copy and rotate
				if any([i in cfg.corr_tensorcomponents for i in horizontals]):
					
					try:
						str1.rotate('NE->RT',back_azimuth=bazs[sp_i])
						
					except:
						print('** data not rotated for station: ')
						print(str1)
						pass
					try:
						str2.rotate('NE->RT',back_azimuth=azms[sp_i])
					except:
						print('** data not rotated for station: ')
						print(str2)
						pass

				# - channel loop
				
				for cpair in self.channel_pairs[sp_i]:
					
					print(cpair)
					loc1, cha1 = cpair[0].split('.')[2:4]
					loc2, cha2 = cpair[1].split('.')[2:4]
					try:
						tr1 = str1.select(location=loc1,channel=cha1)[0]
						tr2 = str2.select(location=loc2,channel=cha2)[0]
					except IndexError:
						continue
					# - check minimum length requirement
					if tr1.stats.npts < min_len_samples:
						continue
					if tr2.stats.npts < min_len_samples:
						continue
					# - correlate
					correlation = cross_covar(tr1.data,tr2.data,
						max_lag_samples,False)[0]
					
					# - add to stack

			# - update time
			t += cfg.time_window_length - cfg.time_overlap


	def update_data(self,t,win_len):

		for trace in self.data:

			# is the trace long enough?
			if trace.stats.endtime < t + win_len:


				# which trace to read next?
				try:
					f = self.inv[trace.id].pop(0)
				except:
					f = None
				# Get the last bit of the trace that we still need
				trace.trim(starttime=t)

				# read new trace
				try:
					newtrace = read(f) 
				except:
					print('** Could not read trace: %s' %f)

				# add new trace to stream
				self.data += newtrace

			else:

				# Only trim -- use no more memory than necessary
				trace.trim(starttime=t)

	def initialize_data(self):

		self.data = Stream()
		
		for channel in self.channels:
			
			f = self.inv[channel].pop(0)
			try:
				self.data += read(f)
			except IOError:
				print('** problems reading file %s' 
				%self.inv.data[channel])

