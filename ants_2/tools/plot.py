
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import os

from obspy import read_inventory

class stainfo(object):

	def __init__(self,staid):

		self.id = staid
		self.lat = None
		self.lon = None

def plot_stations(projection='cea',data='raw',
	channels=['BHZ','LHZ'],locations = ['','00','10']):


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
				resolution='i',
				projection=projection,
				lat_0=mid_paral,
				lon_0=mid_merid,
				llcrnrlat=ymin,
				urcrnrlat=ymax,
				llcrnrlon=xmin,
				urcrnrlon=xmax)

	
	m.bluemarble()

	# plot stations on map
	for sta in stations:
		print sta.id
		m.plot(sta.lon,sta.lat,'rv',markersize=12.,latlon=True)
		x, y = m(sta.lon,sta.lat)
		plt.text(x,y,'   '+sta.id,fontweight='bold',color='w')
	# save map in test folder
	plt.show()



# Some sort of loop function to plot converging stack
#if self.plot:
#			lag = np.linspace(-cfg.corr_maxlag,cfg.corr_maxlag,max_lag_samples)
#			plt.ion()
#			fig = plt.figure()
#			ax = fig.add_subplot(111)
#			line1, = ax.plot(lag, np.zeros(max_lag_samples), '-') # Returns a tuple of line objects, thus the comma
#			plot_pair = self.channel_pairs[0]
#					
#		
#					if self.plot and cpair == plot_pair:
#						line1.set_ydata(correlation)
#						fig.canvas.draw()
# 						
# 