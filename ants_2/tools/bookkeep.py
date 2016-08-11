import os
from glob import glob
from obspy import read

def find_raw_files(indirs, format):
        
    
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
    
    return data