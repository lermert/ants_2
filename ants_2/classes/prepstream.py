from __future__ import print_function
import os

import numpy as np
# Can comment plt out if you are not going to use testrun
import matplotlib.pyplot as plt

import ants_2.tools.prepare as pp
from ants_2.tools.bookkeep import name_processed_file
from obspy import Stream, read, read_inventory, Inventory
from scipy.signal import cheb2ord, cheby2, zpk2sos, sosfilt


class PrepStream(object):

	def __init__(self,filename,ofid=None):

		
		self.stream = read(filename) # This is a stream now
		# Obspy will throw an error if there are any IO problems.
		self.ids = list(set([tr.id for tr in self.stream]))

		self.ofid = ofid 



	def write(self,rdir,cfg):

		
		for trace in self.stream:

			fnew = name_processed_file(
				trace.stats,
				startonly=False)
			fnew = os.path.join(rdir,fnew) 
			
			trace.write(fnew,
				format=trace.stats._format)
			           
			if cfg.verbose:
				print('* renamed file: '+fnew,file = self.ofid)
		return()
		
	def prepare(self,cfg):

# Tasks:
# - merge 
# - if asked, slice 
# - if asked, trim 
# - if asked, add an antialias filter
# - if asked, add a prefilt
# - if asked, add an instr. response
		if cfg.verbose:
			print('* Merging stream', file = self.ofid)
		self.stream = pp.merge_traces(self.stream,
			cfg.Fs_old,5,maxgap=cfg.quality_maxgapsec)

		if len(self.stream) == 0:
			return()

		if cfg.wins_trim:
			if cfg.verbose:
				print('* Trimming to full second', file = self.ofid)
			self.stream = pp.trim_next_sec(self.stream,
				cfg.verbose,self.ofid)
		if len(self.stream) == 0:
			return()


		if cfg.wins:
			if cfg.verbose:
				print('* Slicing stream', file = self.ofid)
			
			self.stream = pp.slice_traces(self.stream,
				cfg.wins_len_sec,cfg.quality_minlensec,
				cfg.verbose,self.ofid)
		
		if len(self.stream) == 0:
			return()

		if cfg.testrun:
			# Retain only three, randomly selected parts of the stream
			sel_ind = np.random.randint(0,len(self.stream),1)[0]
			self.stream = Stream(self.stream[sel_ind:sel_ind+3])
			
		return()

		


		

	def process(self,cfg):
		
		# Preparatory steps
		if cfg.testrun: 
			teststream = Stream(self.stream[0].copy())
			testtitle = ['Raw data']
		
		Fs = self.stream[0].stats.sampling_rate

		if Fs > cfg.Fs_new[-1]:
			self.add_antialias(Fs,cfg.Fs_new[-1]*
				cfg.Fs_antialias_factor)
		
		if cfg.instr_correction:
			self.add_inv(cfg.instr_correction_input,
				cfg.instr_correction_unit)
		
		self.check_nan_inf(cfg.verbose)

		
		if len(self.stream) == 0: 
			return()


		if cfg.wins_detrend:
			if cfg.verbose:
				print('* Detrend', file = self.ofid)
			self.detrend()

		if cfg.wins_demean:
			if cfg.verbose:
				print('* Demean', file = self.ofid)
			self.demean()


		# ToDo: Prettier event excluder
		if cfg.event_exclude:
			# Re-merge for event exclusion
			if cfg.verbose:
				print('* Re-Merging for event exclusion',
					file = self.ofid)
			self.stream._cleanup()
			if cfg.verbose:
				print('* Excluding high energy windows', 
					file = self.ofid)
			# This is run twice
			self.event_exclude(cfg)
			self.event_exclude(cfg)

			if cfg.wins:
				if cfg.verbose:
					print('* Slicing stream', file = self.ofid)
				self.stream = pp.slice_traces(self.stream,
					cfg.wins_len_sec,cfg.quality_minlensec,
					cfg.verbose,self.ofid)
		
		if cfg.wins_taper is not None:
			if cfg.verbose:
				print('* Taper', file = self.ofid)
			self.taper(
				cfg.wins_taper_type,
				cfg.wins_taper
				)

		if cfg.testrun:
			teststream += self.stream[0].copy()
			testtitle.append('After detrend, event exclusion')

		if Fs > cfg.Fs_new[-1]:
			if cfg.verbose:
				print('* Downsample', file = self.ofid)
			self.downsampling(cfg,True)

		if cfg.testrun:
			teststream += self.stream[0].copy()
			testtitle.append('After antialias, downsampling')

		if cfg.instr_correction:
			if cfg.verbose:
				print('* Remove response', file = self.ofid)
			self.remove_response(
				cfg.instr_correction_prefilt,
				cfg.instr_correction_waterlevel,
				cfg.instr_correction_unit,
				cfg.verbose)

		if cfg.testrun:
			teststream += self.stream[0].copy()
			testtitle.append('After instrument correction')
			self.plot_test(teststream, testtitle)

		self.check_nan_inf(cfg.verbose)
		self.stream._cleanup()

		return()



	def plot_test(self,stream,titles):

		# I am a bit unhappy with having this plotting thing here..
		try: 
			os.mkdir('test')
		except:
			pass

		i = 0
		fig_name = stream[0].id + '.testplot.%g.png' %i
		while os.path.exists(os.path.join('test',fig_name)):
			i += 1
			fig_name = stream[0].id + '.testplot.%g.png' %i
		

		fig = plt.figure()

		for i in range(4):
			ax = fig.add_subplot(4,1,i+1)
			ax.plot(stream[i].data)
			ax.set_title(titles[i])

		plt.savefig(os.path.join('test',fig_name),format='png')




	def add_inv(self,input,unit):

		
		
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
				self.inv[id] = {'filename': file,
				'units': unit}
		
		else:
			msg = 'input must be \'resp\' or \'staxml\''
			raise ValueError(msg)





	def add_antialias(self,Fs,freq,maxorder=8):
		# From obspy
		nyquist = Fs * 0.5
		
		# rp - maximum ripple of passband, rs - attenuation of stopband
		rp, rs, order = 1, 96, 1e99
		ws = freq / nyquist  # stop band frequency
		wp = ws  # pass band frequency
		# raise for some bad scenarios
		if ws > 1:
			ws = 1.0
			msg = "** Selected corner frequency is above Nyquist. " + \
				  "** Setting Nyquist as high corner."
			warnings.warn(msg)
		while True:
			if order <= maxorder:
				break
			wp = wp * 0.99
			order, wn = cheb2ord(wp, ws, rp, rs, analog=0)
		
		z, p, k = cheby2(order, rs, wn, 
			btype='low', analog=0, output='zpk')

		self.anti_alias = zpk2sos(z,p,k)
		return()




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

	def detrend(self):
		"""
		remove linear trend
		"""


		self.stream.detrend('linear')



	def demean(self):
		"""
		remove the mean
		"""

		self.stream.detrend('demean')



	def taper(self,ttype,perc):

		    
		self.stream.taper(type=ttype,
			max_percentage=perc)

	
	def event_exclude(self,cfg):

		for trace in self.stream:
			pp.event_exclude(
			    	trace,
			    	windows=cfg.event_exclude_winsec,
			    	n_compare=cfg.event_exclude_n,
			    	freq_min=cfg.event_exclude_freqmin,
			    	freq_max =cfg.event_exclude_freqmax,
			    	factor_enrg=cfg.event_exclude_level,
			    	taper_perc=cfg.wins_taper,
			    	thresh_stdv=cfg.event_exclude_std,
			    	ofid=self.ofid,
			    	verbose=cfg.verbose)
		


	def downsampling(self,cfg,zerophase_antialias):


		# Apply antialias filter
		for trace in self.stream:
			if zerophase_antialias:
				firstpass = sosfilt(self.anti_alias,trace.data)
				trace.data = sosfilt(self.anti_alias,firstpass[::-1])[::-1]
			else:
				trace.data = sosfilt(self.anti_alias,trace.data)
		
		# Decimate if possible, otherwise interpolate
		for Fs in cfg.Fs_new:
			Fs_old = self.stream[0].stats.sampling_rate


			if Fs_old % Fs == 0:
				dec = int( Fs_old / Fs)
				self.stream.decimate(dec, no_filter=True, 
					strict_length=False)
				if cfg.verbose:
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
	    
	    

    

