# Correlation block object:
from obspy import Stream, Trace, read_inventory, UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth

import numpy as np
import os
import re
from glob import glob

from ants_2.tools.preprocess2 import ram_norm, whiten



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

		geo_inf = self.get_geoinf()
		
			
		

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

	def get_geoinf(self):


		inv1 = '{}.{}.xml'.format(*self.id1.split('.')[0:2])
		inv2 = '{}.{}.xml'.format(*self.id2.split('.')[0:2])

		inv1 = read_inventory(os.path.join('meta','stationxml',inv1))
		inv2 = read_inventory(os.path.join('meta','stationxml',inv2))

		# Replace 'radial' and 'transverse' by 'N' and 'E'
		id1 = re.sub('\.??R$','N',self.id1)
		id2 = re.sub('\.??R$','N',self.id2)
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


	def _add_corr(self):
		pass

	def write_tot(self):
		pass

	def write_int(self):
		pass

	

class CorrBlock(object):


# - initialize with station pairs
	def __init__(self,block,inv,cfg):

		self.inv = inv
		self.correlations = []
		self.channels = []

		n_lag = 13

		for pair in block:
			try:
				self.correlations.append(CorrTrace(pair[0],pair[1],
				cfg.corr_type,n_lag))
			except:
				print('** Could not initialize correlation for %s,%s: check metadata'
					%(pair[0],pair[1]))
			
			self.channels.append(pair[0])
			self.channels.append(pair[1])


		self.channels = list(set(self.channels))
		
		
		self.initialize_data()
		self.sampling_rate = self.data[0].stats.sampling_rate
		self.delta = self.data[0].stats.delta

		self.correlate(cfg)

	def correlate(self,cfg):

		t_0 = UTCDateTime(cfg.time_begin)
		t_end = UTCDateTime(cfg.time_end)
		win_len_seconds = cfg.time_window_length
		win_len_samples = round(win_len_seconds*self.sampling_rate)
		min_len_samples = round(cfg.time_min_len*self.sampling_rate)
		


		

		t = t_0
		#t = min(t0,self.data)

		while t <= t_end - (win_len_seconds - self.delta):
			print(t)
			

			
			# - check endtime, if necessary, add data from 'later' file
			self.update_data(t, win_len_seconds)

			# - slice the traces
			windows = self.data.slice(t, t + win_len_seconds - self.delta)
			
			# - Apply preprocessing
			for w in windows:
				if cfg.whiten:
					w = whiten(w,cfg.white_freqmin,cfg.white_freqmax,cfg.white_taper)

				if cfg.onebit:
					w.data = np.sign(w.data)

				if cfg.ram_norm:
					w = ram_norm(w,cfg.ram_window,cfg.ram_prefilt)
				

			# - correlate each relevant pair
			while block:

				pair = block.pop()

				# - select traces
				net1, sta1, loc1, cha1 = pair[0].split('.')
				net2, sta2, loc2, cha2 = pair[1].split('.')

				if 'R' in cha1 or 'R' in cha2 or 'T' in cha1 or 'T' in cha2:

					#str1, str2 = 

				tr1 = windows.select(network = net1,
									 station = sta1,
									 location = loc1,
									 channel = cha1)
				tr2 = windows.select(network = net2,
									 station = sta2,
									 location = loc2,
									 channel = cha2)

				

				# - check minimum length requirement

				# - if horizontal components are involved, copy and rotate
				

			# - add to stack
			# - if window counter reaches n_intermediate_stack: save intermediate


			

			t += cfg.time_window_length - cfg.time_overlap


	def update_data(self,t,win_len):

		for trace in self.data:

			# is the trace long enough?
			if trace.stats.endtime < t + win_len:


				# which trace to read next?
				try:
					f = self.inv.data[trace.id].pop(0)
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
			
			f = self.inv.data[channel].pop(0)
			try:
				self.data += read(f)
			except IOError:
				print('** problems reading file %s' 
				%self.inv.data[channel])

