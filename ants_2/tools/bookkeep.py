import os
from glob import glob
from obspy import UTCDateTime

def find_files(indirs, format):
        
    
    content=list()

    for indir in indirs:
        
        content.extend(glob(os.path.join(indir,'*'+format.lower())))
        content.extend(glob(os.path.join(indir,'*'+format.upper())))
        
    content.sort()
    
    return content

def name_processed_file(stats,startonly=False):
    
    inf = [
        stats.network,
        stats.station,
        stats.location,
        stats.channel
    ]
    
        
    t1=stats.starttime.strftime('%Y.%j.%H.%M.%S')
    t2=stats.endtime.strftime('%Y.%j.%H.%M.%S')
    if startonly: 
        t2 = '*'
    
    inf.append(t1)
    inf.append(t2)

    inf.append(stats._format)
    
    filenew = '{}.{}.{}.{}.{}.{}.{}'.format(*inf)
    
    return filenew
            

def name_correlation_file(sta1,sta2,corr_type,fmt='SAC'):

    name = '{}--{}.{}.{}'.format(sta1,sta2,corr_type,fmt)
    print(name)
    return(name)



class file_avail(object):

    """
    For each station id, keep a dictionary of available files within the time range requested in cfg.
    """

    def __init__(self,cfg):

        self._get_data(cfg)
        self._correlation_blocks(cfg)


    def _get_data(self,cfg):
        


        files = []
        self.stations = {}
        self.data = {}

        indirs = cfg.indirs
        t0 = UTCDateTime(cfg.time_begin) 
        t1 = UTCDateTime(cfg.time_end)
        filefmt = cfg.input_format


        files = find_files(indirs,filefmt)
        
        for f in files:
 
            fn = os.path.basename(f).split('.')
            st = UTCDateTime('{}-{}T{}:{}:{}'.format(*fn[4:9]))
            et = UTCDateTime('{}-{}T{}:{}:{}'.format(*fn[9:14]))
           
            if st > t1 or et < t0:
                continue
            else:

                station = '{}.{}'.format(*fn[0:2])
                channel = '{}.{}.{}.{}'.format(*fn[0:4])

                if station not in self.stations.keys():
                    self.stations[station] = []

                if channel not in self.stations[station]:
                    self.stations[station].append(channel)

                if channel not in self.data.keys():
                    self.data[channel] = []

                self.data[channel].append(f)
        
        return()

    def _channel_pairs(self,sta1,sta2,cfg):


        channels = []
        tensor_comp = cfg.corr_tensorcomponents


        for c1 in self.stations[sta1]:
            for c2 in self.stations[sta2]:


                if cfg.update:
                    f = name_correlation_file(c1,c2,cfg.corr_type)
                    f = os.path.join('data','correlations',f)
                    if os.path.exists(f):
                        continue

                loc1 = c1.split('.')[2]
                loc2 = c2.split('.')[2]

                if loc1 not in cfg.locations:
                    continue

                if loc2 not in cfg.locations:
                    continue

                if loc1 != loc2 and not cfg.locations_mix:
                    continue

                comp = c1[-1] + c2[-1]

                if comp in tensor_comp:
                    channels.append((c1,c2))

        return(channels)

    def _correlation_blocks(self,cfg):
        # ToDo build in updating mode again?
        staids = self.stations.keys()
        # sort alphabetically
        staids.sort()
        blcks_stations = []
        blcks_channels = []
        idprs = []

        n_ids = len(staids)
        n_blk = cfg.n_stationpairs

        for i in range(n_ids):
            for j in range(i+1,n_ids):

                if len(idprs) == n_blk:
                    blcks_stations.append(idprs)
                    idprs = []

                idprs.append((staids[i],staids[j]))

        if len(idprs) < n_blk:
            blcks_stations.append(idprs)


        idprs = []

        for blck in blcks_stations:
            idprs = []
            for pair in blck:

                idprs.extend(self._channel_pairs(pair[0],pair[1],cfg))

                

            if idprs != []:
                blcks_channels.append(idprs)

        self.blocks = blcks_channels


 
        
# def read_data(filepath,ofid,verbose):
    
#     if verbose:
#         print('===========================================================',\
#         file=ofid)
#         print('* opening file: '+filepath+'\n',file=ofid)
        
#     #- read data
#     try:
#         data=read(filepath)

#     except (TypeError, IOError):
#         if verbose: 
#             print('** file wrong type or not found, skip.',
#             file=ofid)
#         return []
#     except:
#         if cfg.verbose: 
#             print('** unexpected read error, skip.',
#             file=ofid)
#         return []
        
#     #- check if this file contains data
#     if len(data) == 0:
#         print('** file contains no data, skip.',file=ofid)
#         return []
    
 #   return data