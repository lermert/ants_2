
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import arange
import os
import h5py
from obspy import read_inventory
import numpy as np
from ants_2.tools.preprocess import bandpass as get_bandpass
from scipy.signal import sosfilt


class stainfo(object):

	def __init__(self,staid):

		self.id = staid
		self.lat = None
		self.lon = None

def plot_stations(projection='merc',data='raw',
	channels=['BHZ','LHZ'],locations = ['','00','10'],bluemarble=False):


	# figure out station IDS and their coordinates
	ids = []
	stations = []
	lats = []
	lons = []

	files = os.listdir(os.path.join('data',data))
	for f in files:
		infs = os.path.basename(f).split('.')
		if len(infs) < 4: continue
		ids.append('{}.{}'.format(*infs[0:2]))
	ids = list(set(ids))

	if ids == []:
		print 'No data found.'
		return()
	
	# look up station xml files, get coordinates
	for i in ids:
		station = stainfo(i)
		
		staxml = os.path.join('meta','stationxml',i+'.xml')
		if os.path.exists(staxml) == False:
			continue

		inv = read_inventory(staxml)

		for j in range(len(locations)):
			for k in range(len(channels)):

				if station.lat is not None: break

				try:
					c = inv.get_coordinates(i+'.'+locations[j]+
						'.'+channels[k])
					station.lat = c['latitude']
					station.lon = c['longitude']
					lats.append(c['latitude'])
					lons.append(c['longitude'])
				except:
					continue
			
		if station.lat == None:
			print 'No coordinates found for station %s\
for locations \'\',00,10.' %i
			continue
		else:
			stations.append(station)
	

	# xmin, xmax, ymin, ymax and central meridian of map
	mid_merid = (max(lons) - min(lons)) * 0.5
	mid_paral = (max(lats) - min(lats)) * 0.5
	xmin = min(lons) - 10
	xmax = max(lons) + 10
	ymin = min(lats) - 5
	ymax = max(lats) + 5

	# basemap
	m = Basemap(rsphere=6378137,
				resolution='l',
				projection=projection,
				lat_0=mid_paral,
				lon_0=mid_merid,
				llcrnrlat=ymin,
				urcrnrlat=ymax,
				llcrnrlon=xmin,
				urcrnrlon=xmax)

	
	if bluemarble:
		m.bluemarble()
		textcol = 'w'
	else:
		m.drawcoastlines()
		textcol = 'k'

	#draw the meridians and parallels
	parallels = arange(round(ymin),round(ymax),10)
	#labels = [left,right,top,bottom]
	m.drawparallels(parallels,labels=[False,True,True,False])
	meridians = arange(round(xmin),round(xmax),20)
	m.drawmeridians(meridians,labels=[True,False,False,True])

	# plot stations on map
	for sta in stations:
		print sta.id
		m.plot(sta.lon,sta.lat,'rv',markersize=12.,latlon=True)
		x, y = m(sta.lon,sta.lat)
		plt.text(x,y,'   '+sta.id,fontweight='bold',color=textcol)
	# save map in test folder
	plt.show()


def plot_converging_stack(inputfile,bandpass=None):

	f = h5py.File(inputfile,'r')
	plt.ion()

	stack = f['corr_windows'].keys()[0]
	stack = f['corr_windows'][stack][:]

	# display a counter for stacked windows
	cnt = 1

	stats = f['stats']
	Fs = stats.attrs['sampling_rate']
	cha1 = stats.attrs['channel1']
	cha2 = stats.attrs['channel2']
	max_lag = ((len(stack) - 1) / 2) / Fs
	lag = np.linspace(-max_lag,max_lag,len(stack))

	fig = plt.figure()
	ax1 = fig.add_subplot(212)

	ax1.set_title('{}--{}'.format(cha1,cha2))
	line1, = ax1.plot(lag,stack,'k')
	

	ax2 = fig.add_subplot(211)
	ax2.set_ylim([np.min(stack)*3,np.max(stack)*3])
	line2, = ax2.plot(lag,stack)
	text1 = ax2.set_title(str(cnt))

	if bandpass is not None:
		sos = get_bandpass(df=Fs,freqmin=bandpass[0],
			freqmax=bandpass[1],
			corners=bandpass[2])
		firstpass = sosfilt(sos, stack)
		stack =  sosfilt(sos, firstpass[::-1])[::-1]

	

	for key in f['corr_windows'].keys():

		cwindow = f['corr_windows'][key][:]

		if bandpass is not None:
			
			firstpass = sosfilt(sos, cwindow)
			cwindow =  sosfilt(sos, firstpass[::-1])[::-1]

		stack += cwindow
		
		ax1.set_ylim([np.min(stack)*1.5,np.max(stack)*1.5])
		text1.set_text(str(cnt))

		
		

		line1.set_ydata(stack)
		line2.set_ydata(cwindow)

		fig.canvas.draw()
		cnt += 1



def plot_window(correlation, window, measurement):
    
    
    maxlag = correlation.stats.npts * correlation.stats.delta
    lag = np.linspace(-maxlag,maxlag,correlation.stats.npts)
    
    plt.plot(lag,correlation.data/np.max(np.abs(correlation.data)))
    plt.plot(lag,window/np.max(np.abs(window)),'--')
    plt.title(correlation.id)
    plt.text(0,-0.75,'Measurement value: %g' %measurement)
    plt.xlabel('Correlation Lag in seconds.')
    plt.ylabel('Normalized correlation and window.')
    
    plt.show()