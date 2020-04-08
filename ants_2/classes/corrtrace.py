from ants_2.tools.util import get_geoinf
from ants_2.tools.bookkeep import name_correlation_file
from obspy import Trace
try:
    import h5py
except ImportError:
    print('h5py cannot be used and if you set output="ASDF"\
    or write_intermediate > 0, processing will fail.')
import pyasdf
import os
import numpy as np
from pympler import tracker

class CorrTrace(object):

    """
    Object holds correlation data along with metainformation (station id, geographic location).
    """

    def __init__(self,cha1,cha2,sampling_rate,corr_type='ccc',
        t0=None,t1=None,stck_int=None,prepstring=None,
        window_length=None,overlap=None,corr_params=None):

        
        self.stack = Trace() # These traces get 01,01,1970 as start date, a completely random choice of start time...no better idea. 
        self.pstak = None
        self.maxlag = None # maxlag will be set the first time a correlation is added.

        # Parameters that must be set 
        self.cnt_tot = 0
        self.cnt_int = 0
        self.id1   = cha1
        self.id2   = cha2
        self.chooser = np.zeros(10)
        self.chooser[1] = 1

        if self.id1[-1] == 'E':
            cha = self.id1.split('.')[-1]
            cha = cha[0] + cha [1] + 'T'
            inf = self.id1.split('.')
            self.id1 = '{}.{}.{}.{}'.format(*(inf[0:3]+[cha]))

        if self.id1[-1] == 'N':
            cha = self.id1.split('.')[-1]
            cha = cha[0] + cha [1] + 'R'
            inf = self.id1.split('.')
            self.id1 = '{}.{}.{}.{}'.format(*(inf[0:3]+[cha]))

        if self.id2[-1] == 'E':
            cha = self.id2.split('.')[-1]
            cha = cha[0] + cha [1] + 'T'
            inf = self.id2.split('.')
            self.id2 = '{}.{}.{}.{}'.format(*(inf[0:3]+[cha]))

        if self.id2[-1] == 'N':
            cha = self.id2.split('.')[-1]
            cha = cha[0] + cha [1] + 'R'
            inf = self.id2.split('.')
            self.id2 = '{}.{}.{}.{}'.format(*(inf[0:3]+[cha]))


        self.id    = self.id1 + '--' + self.id2
        self.corr_type = corr_type

        self.stack.stats.sampling_rate = sampling_rate
        self.sampling_rate = sampling_rate
        (self.stack.stats.network, 
            self.stack.stats.station, 
            self.stack.stats.location, 
            self.stack.stats.channel) = cha1.split('.')

        try:
            geo_inf = get_geoinf(self.id1, self.id2)
            self.lat1 = geo_inf[0]
            self.lat2 = geo_inf[2]
            self.lon1 = geo_inf[1]
            self.lon2 = geo_inf[3]
            self.az   = geo_inf[5]
            self.baz  = geo_inf[6]
            self.dist = geo_inf[4]
        except FileNotFoundError:
            self.lat1 = 0.
            self.lat2 = 0.
            self.lon1 = 0.
            self.lon2 = 0.
            self.az   = 0.
            self.baz  = 0.
            self.dist = 0.


        # Parameters that are optional and will be ignored if they are set to None
        self.stck_int = stck_int
        self.params = corr_params
        self.begin = t0
        self.end   = t1
        self.window_length = window_length
        self.overlap = overlap
        self.prepstring = prepstring
        #self.mytracker = tracker.SummaryTracker()

        # open the file to dump intermediate stack results
        if self.stck_int > 0:
            int_file = '{}.{}.windows.h5'.format(self.id,self.corr_type)
            int_file = os.path.join('data','correlations',int_file)
            int_file = h5py.File(int_file,'a')

            # Save some basic information
            int_stats = int_file.create_dataset('stats',data=(0,))
            int_stats.attrs['sampling_rate']    = self.sampling_rate
            int_stats.attrs['channel1']         = self.id1
            int_stats.attrs['channel2']         = self.id2
            int_stats.attrs['distance']         = self.dist

            # Prepare a group for writing the data window
            self.int_file = int_file
            self.interm_data = int_file.create_group("corr_windows")
        else:
            self.int_file = None
            self.interm_data = None

    def _add_corr(self, corr, t):

        """
        Add one correlation window to the stack
        """

        if self.stack.stats.npts == 0:
            self.stack.data = corr 
            # set the lag
            self.nlag = self.stack.stats.npts
            self.maxlag = (self.nlag - 1) / 2 * self.sampling_rate
        else:
            self.stack.data += corr # This will cause an error if the correlations have different length.
        self.cnt_tot += 1


        if self.stck_int is not None:

            if self.pstak is not None:

                self.pstak += corr  # This will cause an error if the correlations have different length.
            else:
                self.pstak = corr.copy()

            if self.cnt_int == self.stck_int and self.stck_int > 0:
                # write intermediate result
                self.write_int(t)
                self.cnt_int = 0
                self.pstak = None
                self.int_file.flush()


            self.cnt_int += 1

        del corr
        # self.mytracker.print_diff()
        return()


    def write_stack(self, output_format):

        
        # SAC format

        if output_format.upper() == 'SAC':
            filename = os.path.join('data', 'correlations',
                                    '{}.SAC'.format(self.id))

            #- open file and write correlation function
            if self.cnt_tot > 0:
                #- Add metadata
                self.add_sacmeta()
                self.stack.write(filename,format='SAC')
            else:
                print('** Correlation stack contains no windows. Nothing written.')
                print(filename)

        #- ASDF format

        if output_format.upper() == 'ASDF':
            filename = os.path.join('data','correlations','correlations.h5')

            if self.cnt_tot > 0:
                with pyasdf.ASDFDataSet(filename) as ds:
                    info = self.add_asdfmeta()
                    ds.add_auxiliary_data(self.stack.data,
                                          data_type="CrossCorrelation",
                                          path="%s/%s" % (info["trace_id_a"].replace(".", "_"), 
                                                          info["trace_id_b"].replace(".", "_")),
                                          parameters=info)

        if self.int_file is not None:
            self.int_file.file.close()

        
        
        

    def write_int(self,t):

        tstr = t.strftime("%Y.%j.%H.%M.%S")
        #print(self.int_file)
        #print(tstr)
        try:
            self.interm_data.create_dataset(tstr, shape=self.pstak.shape, dtype=self.pstak.dtype)
            self.interm_data[tstr][:] = self.pstak
        except RuntimeError:
            pass
        self.pstak = None
        # flush was not effective.
        # self.int_file.flush()

        close_and_reopen = 1 # np.random.choice(self.chooser)
        if close_and_reopen:
            self.int_file.close()
            int_file = os.path.join('data','correlations','{}.{}.windows.h5'.format(self.id,self.corr_type))
            self.int_file = h5py.File(int_file, "a")
            self.interm_data = self.int_file["corr_windows"]
#
#
    def add_sacmeta(self):

        self.stack.stats.sac={}
        #==============================================================================
        #- Essential metadata  
        #==============================================================================



        self.stack.stats.sac['kt8']     =   self.corr_type
        self.stack.stats.sac['user0']   =   self.cnt_tot

        self.stack.stats.sac['b']       =   -self.maxlag
        self.stack.stats.sac['e']       =   self.maxlag


        self.stack.stats.sac['stla']    =   self.lat1
        self.stack.stats.sac['stlo']    =   self.lon1
        self.stack.stats.sac['evla']    =   self.lat2
        self.stack.stats.sac['evlo']    =   self.lon2
        self.stack.stats.sac['dist']    =   self.dist
        self.stack.stats.sac['az']      =   self.az
        self.stack.stats.sac['baz']     =   self.baz

        self.stack.stats.sac['kuser0']  =   self.id2.split('.')[0]
        self.stack.stats.sac['kuser1']  =   self.id2.split('.')[2]
        self.stack.stats.sac['kuser2']  =   self.id2.split('.')[3]
        self.stack.stats.sac['kevnm']   =   self.id2.split('.')[1]


        #==============================================================================
        #- Optional metadata, can be None  
        #==============================================================================


        self.stack.stats.sac['kt2']     =   self.prepstring
        self.stack.stats.sac['user1']   =   self.window_length
        self.stack.stats.sac['user2']   =   self.overlap

        if self.begin is not None:
            self.stack.stats.sac['kt0'] =   self.begin.strftime('%Y%j')

        if self.end is not None:
            self.stack.stats.sac['kt1'] =   self.end.strftime('%Y%j')

        if self.params is not None: #ToDo: this

            tr.stats.sac['user3'] = self.params[0]
            tr.stats.sac['user4'] = self.params[1]
            tr.stats.sac['user5'] = self.params[2]
            tr.stats.sac['user6'] = self.params[3]
            tr.stats.sac['user7'] = self.params[4]
            tr.stats.sac['user8'] = self.params[5]
    
    
    
    def add_asdfmeta(self):

        info = {}
        info["trace_id_a"] = self.id.split("--")[0]
        info["trace_id_b"] = self.id.split("--")[1]
        info["trace_id_a"] = info["trace_id_a"].replace(" ", "")
        info["trace_id_b"] = info["trace_id_b"].replace(" ", "")
        info["dt"] = round(self.stack.stats.sampling_rate,6)
        info["correlation_type"] = self.corr_type
        info["max_lag"] = self.maxlag
        info["lat_a"] = self.lat1
        info["lat_b"] = self.lat2
        info["lon_a"] = self.lon1
        info["lon_b"] = self.lon2
        info["dist"] = self.dist
        info["az"] = self.az
        info["baz"] = self.baz
        info["preprocess_steps"] = self.prepstring
        info["window_length"] = self.window_length
        info["window_overlap"] = self.overlap

        try:
            info["noisedata_first"] = self.begin.strftime('%Y-%m-%dT%H:%M:%S')
            info["noisedata_last"] = self.end.strftime('%Y-%m-%dT%H:%M:%S')
        except:
            pass

        print(info)

        return info



    
