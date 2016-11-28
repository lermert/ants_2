# Script to stack corr windows, by excluding those above a certain threshold
import numpy as np
from obspy import read
import matplotlib.pyplot as plt
from matplotlib import gridspec
from obspy import Trace,UTCDateTime
from glob import glob
import pandas as pd
import h5py
import os
import sys
from ants_2.tools.util import get_geoinf, rms
from ants_2.tools.windows import snratio
from ants_2.tools.measurements import log_en_ratio, energy

# #==============================================================================
# #-user input
# #==============================================================================

# # file to work on:
# #inputfile = sys.argv[1] 
# # input directory: (all .h5 in this directory will be used)
# indir = sys.argv[1]

# filt = [0.05,0.1,4]

# # set parameters for window choice
# # fixed threshold: Value or None
# threshold1 = 2e-8
# # Variable threshold: Ratio value or None
# threshold2 = 2.0#2.0
# # Variable threshold: Compare to this many other windows
# n_compare = 48
# # frequency bands: union or intersection?
# comb_freq = 'any' 
# # 'all': intersection of weights of different frequency bands is used
# # (only if all frequency bands report 0 weight, the total weight is 0)
# # 'any': if any freq. band has 0 weight, the total weight is 0
# # same thing for the two thresholds
# comb_thre = 'any'
# # the two original traces: union or intersection of weights?
# comb_trac = 'any'
# # plot?
# show_plot = False#True
# # show a grey scale plot of all windows?
# show_plot_allwin = True
# # Start time: Only windows starting at or after will be included
# t_start = '2009-12-31T00:00:00'
# # End time: Only windows ending before this time will be included
# t_end = '2010-09-15T00:00:00'
# # output directory
# stack_dir = 'stack_test'
# # Rayleigh wave group speed in m/s
# r_speed = 2900
# # Max. for color map: this factor times the maximum
# vfac = 0.25
# #==============================================================================
# #==============================================================================


def sacmeta(id1,id2,sampling_rate,nt,cnt):



	sac={}
	#==============================================================================
	#- Essential metadata  
	#==============================================================================

	maxlag = (nt-1)/2 / sampling_rate
	geoinf = get_geoinf(id1,id2)

	#==============================================================================
	sac['user0']	=	cnt

	sac['b']		=	-maxlag
	sac['e']		=	maxlag


	sac['stla']	=	geoinf[0]
	sac['stlo']	=	geoinf[1]
	sac['evla']	=	geoinf[2]
	sac['evlo']	=	geoinf[3]
	sac['dist']	=	geoinf[4]
	sac['az']		=	geoinf[5]
	sac['baz']		=	geoinf[6]

	sac['kuser0']	=	id2.split('.')[0]
	sac['kuser1']	=	id2.split('.')[2]
	sac['kuser2']	=	id2.split('.')[3]
	sac['kevnm']	=	id2.split('.')[1]


	return sac



def max_med_amp(array):

	maxarr = np.max(np.abs(array))
	med = np.median(np.abs(array))

	return maxarr / med

def diff_time(stack,r_speed):

	lag = np.linspace(-stack.stats.sac['e'],stack.stats.sac['e'],
			stack.stats.npts)
	t_pred = stack.stats.sac['dist'] / r_speed
	t_obsr = lag[np.argmax(np.abs(stack.data))]

	return abs(t_pred - abs(t_obsr))




def get_median_values(rms,n_compare):
	meds = np.zeros(rms.shape)

	if np.ndim(rms) > 1:
		cnt1 = rms.shape[1]
	else:
		cnt1 = 1

	for i_f in range(cnt1):
		for i_w in range(rms.shape[0]):

			i0 = i_w-n_compare
			i1 = i_w
			i2 = i_w + 1
			i3 = i_w + n_compare + 1


			while i0 < 0:
				i0 += 1 
				i3 += 1

			while i3 > rms.shape[0]:
				i0 -= 1
				i3 -= 1

			if cnt1 > 1:
				meds[i_w,i_f] = np.median(np.concatenate((rms[i0:i1,i_f],rms[i2:i3,i_f])))	
			else:
				meds[i_w] = np.median(np.concatenate((rms[i0:i1],rms[i2:i3])))	
	return meds

# Tasks:
def ant_stack(input_dir,synth_dir,threshold_fix,threshold_var,threshold_cor,
    n_compare,comb_freq,comb_thre,comb_trac,t_start,t_end,
    t_step,min_win,filt,plot,save,filename):

	#indir,threshold1,threshold2,filt,n_compare,comb_freq,
	#comb_thre,comb_trac,plot,plot_allwin,t_start,t_end,t_step,r_speed,vfac):

	# create output directory
	stack_dir = os.path.join(input_dir,'stack')
	if not os.path.exists(stack_dir):
		os.mkdir(stack_dir)	
	
	
	if filename is None:
		filename = os.path.join(input_dir,'measurement.csv')

	
	
	# Thresholds
	threshold_fix = float(threshold_fix)
	threshold_var = float(threshold_var)
	threshold_cor = float(threshold_cor) if threshold_cor else None

	
	# window params for measurement
	files = glob(os.path.join(input_dir,'*BHZ*BHZ*.h5'))
	g = 3300.
	window_params = {}
	window_params['hw'] = 40.
	window_params['sep_noise']          =  1.
	window_params['win_overlap']        =  False 
	window_params['wtype']              =  'hann'
	window_params['plot']               =  False
	window_params['causal_side']        =  True


	# Save the measurements
	columns = ['sta1','sta2','lat1','lon1','lat2','lon2','dist','az',
    'baz','t0','nstk','snr_c','snr_a','enr_c','enr_a','obs','syn_c','syn_a','syn']
	measurements = pd.DataFrame(columns=columns)
	msr_cnt = 0


	for inputfile in files:

		# measurement to record
		measurement = []

		# open file
		f = h5py.File(inputfile,'r')



		# initialize stack trace --> copy some info from the file stats
		# somehow also need to record the stacking strategy?
		Fs = f['stats'].attrs['sampling_rate']
		id1 = f['stats'].attrs['channel1']
		id2 = f['stats'].attrs['channel2']
		distance = f['stats'].attrs['distance']
		print id1,id2
		geoinf = get_geoinf(id1,id2)

		measurement.append(id1)
		measurement.append(id2)
		measurement.extend(geoinf)

		# try to find synthetic file
		if synth_dir is not None:
			ids_1 = "{}.{}..{}".format(id1.split('.')[0],id1.split('.')[1],
				cha_synth)
			ids_2 = "{}.{}..{}".format(id2.split('.')[0],id2.split('.')[1],
				cha_synth)
			synth_file = os.path.join(synth_dir,"{}--{}.sac".format(ids_1,ids_2))
			tr_synth = read(synth_file)[0]
		else:
			tr_synth = None


		# find nr. of correlation windows
		c_wins = f['corr_windows'].keys()
		n_corrwin = len(c_wins)
		if n_corrwin == 0:
			continue
		n_lag = len(f['corr_windows'][c_wins[0]])
		data = np.zeros(n_lag)
		
		# read out rms values (use clipped rms values? Maybe not worth the complication)
		n_freq = np.shape(f['stats'].attrs['rms_filt'])[0]

		rms1 = np.zeros((n_corrwin,n_freq))
		rms2 = np.zeros((n_corrwin,n_freq))
		rmsc = np.zeros(n_corrwin)
		correlations = np.zeros((n_corrwin,n_lag))
		# msrs1 = np.zeros(n_corrwin)
		# msrs2 = np.zeros(n_corrwin)
		t_windows = np.zeros(n_corrwin)
		if plot:
			wins_all = np.zeros((n_corrwin,n_lag))

		i = 0

###############################################################################
# READ IN DATA and trace rms; determine correlation rms and measurement; and start times
###############################################################################

		for k in c_wins:

			correlations[i,:] = f['corr_windows'][k][:]
			rms1[i,:] = f['corr_windows'][k].attrs['rms1']
			rms2[i,:] = f['corr_windows'][k].attrs['rms2']

			# starttimes
			#tstr = '{},{},{},{},{}'.format(*fh['corr_windows'].keys()[0].split('.'))
			tstr = '{},{},{},{},{}'.format(*k.split('.'))
			t_windows[i] = UTCDateTime(tstr).timestamp

			# correlation rms
			t_temp = Trace(data=correlations[i,:])
			t_temp.stats.sampling_rate = Fs
			t_temp.stats.sac = {}
			t_temp.stats.sac['dist'] = distance
			
			if filt:
				t_temp.taper(0.05)
				t_temp.filter('bandpass',freqmin=filt[0],freqmax=filt[1],
				corners=filt[2],zerophase=True)

			rmsc[i] = rms(t_temp.data)
			# msrs1[i] = energy(t_temp,g,window_params)#log_en_ratio(t_temp,3300,window_params)
			
			# # Getting the acausal energy
			# window_params['causal_side'] = False
			# msrs2[i] = energy(t_temp,g,window_params)

			# window_params['causal_side'] = True
			# msrs3[i] = log_en_ratio(t_temp,g,window_params)
			# if plot:
			# 	wins_all[i,:] = t_temp.data
			# measurements
			# Type of measurement now hardcoded, but should be option later on ToDo

			i+=1


###############################################################################
# OBTAINING THE SELECTION ARRAYS
###############################################################################
		
		# obtain the weight array
		# first strategy weights: keep windows with rms up to the threshold
		weights1_1 = (rms1 <= threshold_fix).astype(int)
		weights1_2 = (rms2 <= threshold_fix).astype(int)

		# second strategy weights
		meds1 = get_median_values(rms1,n_compare)
		meds2 = get_median_values(rms2,n_compare)
		weights2_1 = ((rms1 / meds1) <= threshold_var).astype(int)
		weights2_2 = ((rms2 / meds2) <= threshold_var).astype(int)

		# # third strategy weights
		# medsc = get_median_values(rmsc,n_compare)
		# weights3 = ((rmsc / medsc) <= threshold_cor).astype(int)

		# combine freq. bands
		if comb_freq == 'any':
			weights1_1 = np.prod(weights1_1,axis=1)
			weights1_2 = np.prod(weights1_2,axis=1)
			weights2_1 = np.prod(weights2_1,axis=1)
			weights2_2 = np.prod(weights2_2,axis=1)


		elif comb_freq == 'all':
			weights1_1 = np.sum(weights1_1,axis=1)
			weights1_1 = np.clip(weights1_1,a_min=0.,a_max=1.)
			weights2_1 = np.sum(weights2_1,axis=1)
			weights2_1 = np.clip(weights2_1,a_min=0.,a_max=1.)
			weights1_2 = np.sum(weights1_2,axis=1)
			weights1_2 = np.clip(weights1_2,a_min=0.,a_max=1.)
			weights2_2 = np.sum(weights2_2,axis=1)
			weights2_2 = np.clip(weights2_2,a_min=0.,a_max=1.)

		# combine first and second weights
		if comb_thre == 'any':
			weights1 = np.multiply(weights1_1,weights1_2)
			weights2 = np.multiply(weights2_1,weights2_2)
		elif comb_thre == 'all':
			weights1 = weights1_1 + weights1_2
			weights2 = weights2_1 + weights2_2
			weights1 = np.clip(weights1,0.,1.)
			weights2 = np.clip(weights2,0.,1.)

		# combine first and second trace
		if comb_trac == 'any':
			weights = np.multiply(weights1,weights2)
		elif comb_trac == 'all':
			weights = weights1 + weights2
			weights = np.clip(weights,0.,1.)

		
###############################################################################
# APPLY THE FIRST STAGE OF WEIGHTS
###############################################################################

	
		correlations = np.multiply(correlations.T,weights).T


		
###############################################################################
# OBTAIN CORRELATION BASED WEIGHTS
###############################################################################
		i = 0
		for k in c_names:

			if np.max(np.abs(correlations[i,:])) == 0.:
				continue
			# correlation rms
			t_temp = Trace(data=correlations[i,:])
			t_temp.stats.sampling_rate = Fs
			t_temp.stats.sac = {}
			t_temp.stats.sac['dist'] = distance
			
			if filt:
				t_temp.taper(0.05)
				t_temp.filter('bandpass',freqmin=filt[0],freqmax=filt[1],
				corners=filt[2],zerophase=True)

			rmsc[i] = rms(t_temp.data)
			# msrs1[i] = energy(t_temp,g,window_params)#log_en_ratio(t_temp,3300,window_params)
			# msrs2[i] = log_en_ratio(t_temp,g,window_params)
			

			i+=1

		# third strategy weights
		medsc = get_median_values(rmsc,int(n_compare/5))
		weights3 = ((rmsc / medsc) <= threshold_cor).astype(int)

###############################################################################
# APPLY SECOND STAGE WEIGHTS
###############################################################################

		# combine weights and correlation weights
		weights = np.multiply(weights,weights3)
		correlations = np.multiply(correlations.T,weights).T
		# msrs1 = np.array(msrs1)
		# msrs2 = np.array(msrs2)

###############################################################################
# PLOT after selection...
###############################################################################
		# if plot:

		# 	wins_all = np.multiply(wins_all.T,weights).T


		# 	maxlag = (n_lag-1)/2
		# 	v = np.max(np.max(wins_all)) * 0.3

			          
		# 	x = np.linspace(-maxlag,maxlag,n_lag)
		# 	y = np.linspace(1,n_corrwin, n_corrwin)
		# 	xv, yv = np.meshgrid(x, y)

		# 	fig = plt.figure(figsize=(15,5))
		# 	gs = gridspec.GridSpec(1, 3, width_ratios=[5, 1, 1]) 
		# 	ax0 = fig.add_subplot(gs[0])	
		# 	im = ax0.pcolormesh(xv,yv,wins_all[:-1,:-1],
		# 		cmap=plt.cm.RdGy,vmin=-v,vmax=v)
		# 	ax0.plot(np.ones(n_corrwin)*(distance/g)-window_params['hw'],
		# 		np.arange(n_corrwin))
		# 	ax0.plot(np.ones(n_corrwin)*(distance/g)+window_params['hw'],
		# 		np.arange(n_corrwin))
		# 	ax0.axis('tight')
		# 	#plt.colorbar(im,position='left')

			
		# 	ax1 = fig.add_subplot(gs[1],sharey=ax0)
			
		# 	# rmsc = [np.nan if val == 0 else val for val in rmsc*weights]
		# 	# ax1.plot(rmsc,np.arange(len(rmsc)),'x')
		# 	# ax1.set_xticks([])

		# 	ax2 = fig.add_subplot(gs[2],sharey=ax0)
			
		# 	msrs1 = [np.nan if val == 0 else val for val in msrs1*weights]
		# 	msrs2 = [np.nan if val == 0 else val for val in msrs2*weights]
		# 	ax1.plot(msrs1,np.arange(len(msrs1)),'x')
		# 	ax2.plot(msrs2,np.arange(len(msrs2)),'rx')
		# 	#ax1.set_xticks([0.5*np.max(msrs1),np.max(msrs1)])
		# 	ax1.set_yticks([])
		# 	#ax2.set_xticks([0.5*np.max(msrs2),np.max(msrs2)])
		# 	ax2.set_yticks([])
		# # show or save the plots. 
		# 	plt.tight_layout()
		# 	plt.show()

###############################################################################
# TIME LOOP FORMS INTERMEDIATE STACKS
###############################################################################


		# go through all traces
		
		stacks = []
		tstarts = []
		tends = []
		cnts_good = []
		if t_end is not None:
			t_end = UTCDateTime(t_end).timestamp
		else:
			t_end = t_windows[-1]
		if t_start is not None:
			t_start = max(t_windows[0],UTCDateTime(t_start).timestamp)
		else:
			t_start = t_windows[0]
		if t_step is not None:
			t_step = float(t_step)
		else:
			t_step = t_end - t_start

		t = t_start
		

		#for i in range(n_corrwin):
		while (t + t_step) <= t_end:

			i0 = np.argmin(np.abs(t_windows-t))
			#print t_step
			i1 = np.argmin(np.abs(t_windows-(t+t_step)))
			if t_windows[i1] > t_windows[i0] + t_step:
				i1 -= 1
			#print i0
			#print i1
			#print UTCDateTime(t_windows[i1])


			cnt_good = weights[i0:i1].sum()

			if cnt_good < min_win:
				t += t_step
				continue

			st = correlations[i0:i1].sum(axis=0)/cnt_good
			stacks.append(st)
			tstarts.append(t)
			tends.append(t+t_step)
			cnts_good.append(cnt_good)


			t += t_step
			#corr = f['corr_windows'][c_names[i]]

			#tstr = '{},{},{},{},{}'.format(*c_names[i].split('.'))

			#if UTCDateTime(tstr) < t_start:
		# 		#print 'Before start time'
		# #		continue
		# 	if UTCDateTime(tstr) >= t_end_end:
		# 		print 'after end time'
		# 		break
			
		# 	if UTCDateTime(tstr) >= t_end:
				
		# 		# weight the stack by stack length N 
		# 		if cnt_good >= min_win:			
		# 			stacks.append(data / cnt_good)
		# 			tstarts.append(t_start)
		# 			tends.append(t_end)
		# 			cnts_good.append(cnt_good)
				
		# 		# next stack
		# 		cnt_good = 0
		# 		t_end += t_step
		# 		t_start += t_step
		# 		data = np.zeros(n_lag)
				



		# 	if weights[i] == 1:

		# 		if threshold_cor is None:
					
		# 			# put onto stack
		# 			data += corr
		# 			cnt_good += 1.
					
		# 		elif threshold_cor is not None  and rms(corr) <= threshold_cor:
					
		# 			# put onto stack
		# 			data += corr
		# 			cnt_good += 1.
					
		# 	else:
		# 		pass

###############################################################################
# TAKE MEASUREMENTS ON THE INTERMEDIATE STACKS, PREPARE PLOT
###############################################################################


		wins_all = np.zeros((len(stacks),n_lag))
		
		msrs1 = np.zeros(len(stacks))
		msrs2 = np.zeros(len(stacks))
		snrs1 = np.zeros(len(stacks))
		snrs2 = np.zeros(len(stacks))
		msrs3 = np.zeros(len(stacks))

		for i in range(len(stacks)):

			
			
			t_temp = Trace(data=stacks[i])
			t_temp.stats.sampling_rate = Fs
			t_temp.stats.sac = {}
			t_temp.stats.sac['dist'] = distance
			
			if filt:
				t_temp.taper(0.05)
				t_temp.filter('bandpass',freqmin=filt[0],freqmax=filt[1],
			corners=filt[2],zerophase=True)

			# causal energy msr
			window_params['causal_side'] = True
			msrs1[i] = energy(t_temp,g,window_params)#log_en_ratio(t_temp,3300,window_params)
			snr_c = snratio(t_temp,g,window_params)
			
			# acausal energy msr
			window_params['causal_side'] = False
			msrs2[i] = energy(t_temp,g,window_params)
			snr_a = snratio(t_temp,g,window_params)
			
			window_params['causal_side'] = True # this value should not be used during log_en_ratio but just to be sure
			msrs3[i] = log_en_ratio(t_temp,g,window_params)

			wins_all[i,:] = t_temp.data

			

			stack = Trace()
			stack.stats.sampling_rate = Fs
			stack.stats.sac = {}
			stack.stats.sac.dist = distance
			stack.data = stacks[i]

		# if there is a synthetic file, get a measurement for that too
		if tr_synth is not None:
			tr_synth.stats.sac.dist = distance
			if filt:
				tr_synth.taper(0.05)
				tr_synth.filter('bandpass',freqmin=filt[0],freqmax=filt[1],
			corners=filt[2],zerophase=True)

			# causal energy msr
			window_params['causal_side'] = True
			msrs_syn1 = energy(tr_synth,g,window_params)#log_en_ratio(t_temp,3300,window_params)
		
			# acausal energy msr
			window_params['causal_side'] = False
			msrs_syn2 = energy(tr_synth,g,window_params)
			
			window_params['causal_side'] = True # this value should not be used during log_en_ratio but just to be sure
			msrs_syn3 = log_en_ratio(tr_synth,g,window_params)

		else:
			msrs_syn1 = np.nan
			msrs_syn2 = np.nan
			msrs_syn3 = np.nan
###############################################################################
# Record the measurements for this station pair
###############################################################################	
			
			if not np.nan in [msrs1[i],msrs2[i],msrs3[i],snr_a,snr_c]:
				msrinf = [str(UTCDateTime(tstarts[i])),
					  cnts_good[i],snr_c,snr_a,msrs1[i],msrs2[i],
					  msrs3[i],msrs_syn1,msrs_syn2,msrs_syn3]
			
				measurements.loc[msr_cnt] = measurement + msrinf
			
				msr_cnt += 1



			if save:
				# save the stack.
				stackfile = os.path.splitext(os.path.basename(inputfile))[0]
				stackfile = os.path.splitext(stackfile)[0]
				stackfile = os.path.splitext(stackfile)[0]
				stackfile = stackfile + '.stack.{}.{}.SAC'.format(
					UTCDateTime(tstarts[i]),
					UTCDateTime(tends[i]))
				stackfile = os.path.join(stack_dir,stackfile)
				#print cnts_good[i]
				stack.stats.sac = sacmeta(id1,id2,Fs,stack.stats.npts,cnts_good[i])

				#tdiff = diff_time(stack,r_speed)

				stack.stats.station = id1.split('.')[1]
				stack.stats.location = id1.split('.')[2]
				stack.stats.channel = id1.split('.')[3]
				stack.stats.network = id1.split('.')[0]
				#stack.stats.sac['user1'] = max_med_amp(stack.data)
				#stack.stats.sac['user2'] = len(weights)
				#stack.stats.sac['user3'] = tdiff
				stack.write(stackfile,format='SAC')

	


###############################################################################
# Plot
###############################################################################	

		if plot and len(stacks) > 0:
			maxlag = (n_lag-1)/2
			v = np.max(np.max(wins_all)) * 0.75

			x = np.linspace(-maxlag,maxlag,n_lag)
			y = np.linspace(0,len(stacks), len(stacks)+1)
			xv, yv = np.meshgrid(x, y)

			fig = plt.figure(figsize=(15,5))
			gs = gridspec.GridSpec(1, 3, width_ratios=[5, 1, 1]) 
			ax0 = fig.add_subplot(gs[0])	
			im = ax0.pcolormesh(xv,yv,wins_all,
				cmap=plt.cm.RdGy,vmin=-v,vmax=v)
			ax0.plot(np.ones(len(stacks)+1)*(distance/g)-window_params['hw'],
				np.arange(0,len(stacks)+1))
			ax0.plot(np.ones(len(stacks)+1)*(distance/g)+window_params['hw'],
				np.arange(0,len(stacks)+1))
			#ax0.axis('tight')
			#plt.colorbar(im,position='left')
			
			
			ax1 = fig.add_subplot(gs[1],sharey=ax0)
			
			# rmsc = [np.nan if val == 0 else val for val in rmsc*weights]
			# ax1.plot(rmsc,np.arange(len(rmsc)),'x')
			# ax1.set_xticks([])

			ax2 = fig.add_subplot(gs[2],sharey=ax0,sharex=ax1)
			
			m1 = [np.nan if val == 0 else val for val in msrs1]
			m2 = [np.nan if val == 0 else val for val in msrs2]
			ax1.plot(m2,np.arange(len(m1))+0.5,'d')
			ax2.plot(m1,np.arange(len(m2))+0.5,'rd')
			#ax1.set_xticks([0.5*np.max(msrs1),np.max(msrs1)])
			plt.yticks(np.arange(0,len(m1),12)+0.5,[UTCDateTime(ts).strftime("%Y.%jT%H")
			 for ts in tstarts][0::12])
			#ax2.set_xticks([0.5*np.max(msrs2),np.max(msrs2)])
			#ax2.set_yticks([])
			#ax1.set_yticks([])

			ax0.set_title(id1+'--'+id2)
			ax1.set_title('sig.enr -')
			ax2.set_title('sig.enr +')
			
		# show or save the plots. 
			plt.tight_layout()
			#plt.show()
			pltname = os.path.splitext(inputfile)[0]+'.png'
			plt.savefig(pltname)
	print 'SAVING measurements:'
	#filename = filen+'.measurements.csv'
	print filename
	measurements.to_csv(filename,index=None)

	f = os.path.splitext(filename)[0]	
	with open(os.path.join(f+'.stackinfo.txt'),'w') as fh:
		fh.write('threshold_fix: ')
		fh.write(str(threshold_fix))
		fh.write('\nthreshold_var: ')
		fh.write(str(threshold_var))
		fh.write('\nthreshold_cor: ')
		fh.write(str(threshold_cor))
		fh.write('\nn_compare: ')
		fh.write(str(n_compare))
		fh.write('\ncomb_thre: ')
		fh.write(str(comb_thre))
		fh.write('\ncomb_freq: ')
		fh.write(str(comb_freq))
		fh.write('\ncomb_trac: ')
		fh.write(str(comb_trac))
		fh.write('\nt_start: ')
		fh.write(str(t_start))
		fh.write('\nt_end: ')
		fh.write(str(t_end))
		fh.write('\nt_step (seconds): ')
		fh.write(str(t_step))



