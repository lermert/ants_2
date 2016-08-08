from __future__ import print_function
import os
import ants_2.tools.prepare as pp

import numpy as np
from obspy import Stream, read, read_inventory, Inventory
from scipy.signal import cheb2ord, cheby2, zpk2sos

class PrepStream(object):

	def __init__(self,filename,ofid=None):

		
		self.stream = read(filename) # This is a stream now
		# Obspy will throw an error if there are any IO problems.
		self.ids = list(set([tr.id for tr in self.stream]))

		
		# Check if the ids are all compatible? Actually, have to sort out different channels here.
		self.ofid = ofid 



	def write(self):
		pass

	def prepare(self,cfg):

# Tasks:
# - merge 
# - if asked, slice 
# - if asked, trim 
# - if asked, add an antialias filter
# - if asked, add a prefilt
# - if asked, add an instr. response

		self.stream = pp.merge_traces(self.stream,
			cfg.Fs_old,5,maxgap=cfg.quality_maxgapsec)

		if cfg.quality_trimfullsec:
			self.stream = pp.trim_next_sec(self.stream,
				cfg.verbose,self.ofid)


		if cfg.wins:
			self.stream = pp.slice_traces(self.stream,
				cfg.wins_len_sec,cfg.quality_minlensec,
				cfg.verbose,self.ofid)


		Fs = self.stream[0].stats.sampling_rate

		if Fs > cfg.Fs_new[-1]:
			self.add_antialias(Fs,cfg.Fs_new[-1]*
				cfg.Fs_antialias_factor)

		if cfg.instr_correction:

			#self.pre_filt = cfg.instr_correction_prefilt
			self.add_inv(cfg.instr_correction_input)

	def process(self,cfg):
		
		self.check_nan_inf(cfg.verbose)

		if cfg.wins_detrend:
			self.detrend(cfg.verbose)

		if cfg.wins_demean:
			self.demean(cfg.verbose)

		# ToDo: Prettier event excluder
		if cfg.event_exclude:
		    self.event_exclude(cfg)

		if cfg.wins_taper is not None:
			self.taper(
				cfg.wins_taper_type,
				cfg.wins_taper,
				cfg.verbose
				)

		Fs = self.stream[0].stats.sampling_rate
		if Fs > cfg.Fs_new[-1]:
			self.downsampling(cfg.verbose)

		if cfg.instr_correction:
			self.remove_response(
				cfg.instr_correction_prefilt,
				cfg.instr_correction_waterlevel,
				cfg.instr_correction_unit,
				cfg.verbose)



	

	def add_inv(self,input):

		
		
		if input == 'staxml':

			inf = self.ids[0].split('.')[0:2]
			file = '{}.{}.xml'.format(*inf)
			file = os.path.join('meta','stationxml',file)

			self.inv = read_inventory(file)



		elif input == 'resp':

			self.inv = {}

			for id in self.ids:
				inf = id.split('.')
				file = 'RESP.{}.{}.{}.{}'.format(*inf)
				file = os.path.join('meta','resp',file)
				self.inv[id] = {'filename': file}
		
		else:
			msg = 'input must be \'resp\' or \'staxml\''
			raise ValueError(msg)

	def add_antialias(self,Fs,freq):
		# From obspy
		nyquist = Fs * 0.5
		# rp - maximum ripple of passband, rs - attenuation of stopband
		rp, rs, order = 1, 96, 1e99
		ws = freq / nyquist  # stop band frequency
		wp = ws  # pass band frequency
		# raise for some bad scenarios
		if ws > 1:
		    ws = 1.0
		    msg = "Selected corner frequency is above Nyquist. " + \
		          "Setting Nyquist as high corner."
		    warnings.warn(msg)
		while True:
		    if order <= maxorder:
		        break
		    wp = wp * 0.99
		    order, wn = cheb2ord(wp, ws, rp, rs, analog=0)
		antialias = cheby2(order, rs, wn, 
			btype='low', analog=0, output='zpk')
		self.antialias = zpk2sos(antialias)


	def check_nan_inf(self,verbose):
		"""
		Check if trace contains nan, inf and takes them out of the stream
		"""
		for i in range(len(self.stream)):

			trace = self.stream[i]

			#- check NaN
			if True in np.isnan(trace.data):
			    if verbose: 
			        print('** trace contains NaN, discarded',\
			    file=self.ofid)
			    
			    del_trace = self.stream.pop(i)
			    print(del_trace,file=self.ofid)
			    continue
			            
			#- check infinity
			if True in np.isinf(trace.data):
			    if verbose: print('** trace contains infinity, discarded',
			    file=self.ofid)
			    
			    del_trace = self.stream.pop(i)
			    print(del_trace,file=self.ofid)
			    continue 

	def cap_glitches(trace,cfg):
		pass

	def detrend(self,verbose):
		"""
		remove linear trend
		"""

		if verbose:
		    print('* detrend\n',file=self.ofid)

		self.stream.detrend('linear')



	def demean(self,verbose):
		"""
		remove the mean
		"""

		if verbose: 
		    print('* demean\n',file=self.ofid)

		self.stream.detrend('demean')



	def taper(self,ttype,perc,verbose):

		if verbose:
		    print('* tapering\n',file=self.ofid)
		    
		self.stream.taper(type=ttype,
			max_percentage=perc)

	
	def event_exclude(self,cfg):

		for trace in self.stream:
			pp.event_exclude(
			    	trace,
			    	windows=cfg.event_exclude_winsec,
			    	n_compare=cfg.event_exclude_n,
			    	min_freq=cfg.event_exclude_freq,
			    	factor_enrg=cfg.event_exclude_level,
			    	taper_perc=cfg.event_exclude_taper,
			    	thresh_stdv=cfg.event_exclude_std,
			    	ofid=self.ofid,
			    	verbose=cfg.verbose)
		if cfg.verbose:
			print('* excluded high energy events\n',file=self.ofid)




	def downsampling(self,verbose):


		# Apply antialias filter

		for trace in self.stream:

		    if zerophase_antialias:
		        firstpass = sosfilt(self.antialias,trace.data)
		        trace.data = sosfilt(self.antialias,firstpass[::-1])[::-1]
		    else:
		        trace.data = sosfilt(self.antialias,trace.data)

		# Decimate if possible, otherwise interpolate
		for Fs in cfg.Fs_new:
		    
		    Fs_old = self.stream[0].stats.sampling_rate
		    
		    
		    if Fs_old % Fs == 0:
		        dec = int( Fs_old / Fs)
		        self.stream.decimate(dec, no_filter=True, 
		        	strict_length=False)
		        if verbose:
		            print('* decimated traces to %g Hz' %Fs,
		            file=self.ofid)
		    else:
		        try:
		            self.stream.interpolate(sampling_rate = Fs,
		            method='lanczos')
		            print('* interpolated traces to %g Hz' %Fs,
		            file=self.ofid)
		        except:
		            self.stream.interpolate(sampling_rate = Fs)
		            print('* interpolated trace to %g Hz' %Fs,
		            file=self.ofid)



	def remove_response(self,pre_filt,waterlevel,unit,verbose):


		if isinstance(self.inv,dict):

			for trace in self.stream:
				inv = self.inv[trace.id]
				self.stream.simulate(
					paz_remove=None,
					pre_filt=pre_filt,
					seedresp=inv,
					sacsim=True,
					pitsasim=False,
					water_level=waterlevel,
					output=unit)
			if verbose:
			    print('* removed instrument response using seedresp',file=self.ofid)
		        
		elif isinstance(self.inv,Inventory):

			self.stream.remove_response(
				inventory=self.inv,
				pre_filt=pre_filt,
				water_level=waterlevel,
			    output=unit)
			if verbose:
				print('* removed instrument response using stationxml inv',file=self.ofid)
		        
		else:
			msg = 'No inventory or seedresp found.'
			raise ValueError(msg)
	    
	    

    

