
# Correlation block object:
from obspy import Stream, UTCDateTime, read
from scipy.signal import sosfilt
import numpy as np
import re

from ants_2.tools.util import get_geoinf
from ants_2.classes.corrtrace import CorrTrace
from ants_2.tools.correlations import cross_covar, interference, pcc_2, deconv_waterlevel
from ants_2.tools.treatment import ram_norm, whiten, cap, bandpass
# list of possible channels combinations
horizontals = ['RR', 'RT', 'TR', 'TT', 'TZ', 'ZT', 'RZ', 'ZR']
import os, psutil
from pympler import muppy, summary, tracker
import gc

class CorrBlock(object):

    # - initialize with station pairs
    def __init__(self, block, cfg):

        self.cfg = cfg
        self._correlations = {}
        self.inv = block.inventory
        self.readtimes = block.readtimes
        self.channels = block.channels
        self.station_pairs = block.station_pairs
        self.channel_pairs = block.channel_pairs
        t0 = UTCDateTime(self.cfg.time_begin)
        self.initialize_data_new(t0)
        self.t0 = t0
        self.sampling_rate = self.data[0].stats.sampling_rate
        self.delta = self.data[0].stats.delta
        if self.cfg.bandpass is not None:
            fmin = self.cfg.bandpass[0]
            fmax = self.cfg.bandpass[1]
            if fmax > 0.5 * self.sampling_rate:
                raise ValueError("Upper filter freq above Nyquist.")
            if fmax <= fmin:
                msg = "Bandpass upper corner frequency must be above lower corner frequency."
                raise ValueError(msg)

            order = self.cfg.bandpass[2]
            self.sos = bandpass(freqmin=fmin, freqmax=fmax,
                                df=self.sampling_rate, corners=order)

        for cp in self.channel_pairs:
            for pair in cp:
                if self.cfg.rotate:
                    pair = [re.sub('E$', 'T', str) for str in pair]
                    pair = [re.sub('N$', 'R', str) for str in pair]
                cp_name = '{}--{}'.format(*pair)
                prepstring = self.get_prepstring()

                self._correlations[cp_name] =\
                    CorrTrace(pair[0], pair[1], self.sampling_rate,
                              stck_int=cfg.interm_stack, prepstring=prepstring,
                              window_length=cfg.time_window_length,
                              corr_type=cfg.corr_type, maxlag=cfg.corr_maxlag,
                              t0=self.t0,
                              t1=UTCDateTime(cfg.time_end),
                              overlap=cfg.time_overlap, corr_params=None,
                              rotate=cfg.rotate)

        if any([i in cfg.corr_tensorcomponents for i in horizontals]):
            self.baz1 = []
            self.baz2 = []
            for pair in self.channel_pairs[0]:
                try:
                    geoinf = get_geoinf(pair[0], pair[1])
                    # find azimuth, backazimuth
                    self.baz2.append(geoinf[5])
                    self.baz1.append(geoinf[6])
                except FileNotFoundError:
                    self.baz2.append(0)
                    self.baz1.append(0)


    def run(self, output_file=None):
        process = psutil.Process(os.getpid())
        print('Working on station pairs:')
        for sta in self.station_pairs:
            print("{}--{}".format(sta[0], sta[1]))
        
        t_end = UTCDateTime(self.cfg.time_end)
        wl = self.cfg.time_window_length
        min_len_samples = int(self.cfg.time_min_window * self.sampling_rate)
        max_lag_samples = int(round(self.cfg.corr_maxlag * self.sampling_rate))

        # Time loop
        # running time:
        self.t = self.t0
        update_time = self.t
        self.updated_times = [self.t0]

        # mytracker = tracker.SummaryTracker()
        while self.t < t_end:
            # print(t, file=output_file, end="\n")
            #print("Memory usage in Gb loop begin ", process.memory_info().rss / 1.e9, 
            #      file=output_file, end="\n")
            if np.all([len(self.inv[cc])==0 for cc in self.channels]):
                break
            if len(self.data) == 0:
                break

            # slide
            # the offset is so that we stay on a fixed time stepping
            # obspy slide is untrustworthy and includes partial windows even if told not to
            windows = self.data.slide(wl - self.delta, wl - self.cfg.time_overlap,
                                      offset=self.t - self.data[0].stats.starttime, include_partial_windows=True,
                                      nearest_sample=True)
            possible_update_times = []
            ixw = 0
            for w in windows:
                # ideal starttime, not what actually happens
                w_starttime = self.t + ixw * (wl - self.cfg.time_overlap)
                if np.all([ww.stats.starttime > t_end for ww in w]):
                    self.t = t_end
                    break

                # Apply preprocessing
                w = self.preprocess(w)
                # may return a deepcopy if non-linear processing is applied.
                # - station pair loop
                # starttimes_corrwindows = []
                print(w)

                for sp_i, pair in enumerate(self.station_pairs):

                    # - select traces
                    [net1, sta1] = pair[0].split('.')
                    [net2, sta2] = pair[1].split('.')
                    str1 = w.select(network=net1, station=sta1)
                    str2 = w.select(network=net2, station=sta2)
                     
                    # - if horizontal components are involved, copy and rotate
                    if any([i in self.cfg.corr_tensorcomponents
                            for i in horizontals]):

                        str1, str2 = self.rotate(str1, str2,
                                                 self.baz1[sp_i],
                                                 self.baz2[sp_i])
                    for cpair in self.channel_pairs[sp_i]:

                        if self.cfg.rotate:
                            cpair = [re.sub('E$', 'T', str) for str in cpair]
                            cpair = [re.sub('N$', 'R', str) for str in cpair]

                        cp_name = '{}--{}'.format(*cpair)

                        loc1, cha1 = cpair[0].split('.')[2:4]
                        loc2, cha2 = cpair[1].split('.')[2:4]

                        try:
                            tr1 = str1.select(location=loc1, channel=cha1)[0]
                            tr2 = str2.select(location=loc2, channel=cha2)[0]
                        except IndexError:
                            print("Channel not found, current streams: ", file=output_file)
                            print(str1)
                            print(str2)
                            print(" channels needed: " + cha1 + "," + cha2,
                                  file=output_file)
                            print("no data in window")
                            possible_update_times.append(w_starttime)
                            continue

                        # - check minimum length requirement
                        # - Quite often not fulfilled due to data gaps
                        traces_ok = self.perform_checks(tr1, tr2, output_file,
                                                        min_len_samples)
                        if not traces_ok:
                            print("found issue with one or both traces")
                            possible_update_times.append(w_starttime)
                            continue

                        if self.cfg.corr_type == 'ccc':
                            # correlate
                            correlation = cross_covar(tr1.data, tr2.data,
                                                      max_lag_samples,
                                                      self.cfg.corr_normalize)[0]
                        elif self.cfg.corr_type == 'mic':
                            correlation = interference(tr1.data, tr2.data,
                                                       max_lag_samples)
                        elif self.cfg.corr_type == 'pcc':
                            correlation = pcc_2(tr1.data, tr2.data,
                                                max_lag_samples)[0]
                        elif self.cfg.corr_type == "dcv":
                            correlation = deconv_waterlevel(tr1.data, tr2.data,
                                max_lag_samples)


                        # add to stack
                        if len(correlation) == 2 * max_lag_samples + 1:
                            actual_window_length = tr1.stats.delta * tr1.stats.npts
                            written = self._correlations[cp_name]._add_corr(correlation,
                                                                  tr1.stats.starttime,
                                                                  actual_window_length)
                            #starttimes_corrwindows.append(tr1.stats.starttime)
                            del correlation
                        else:
                            print('Empty window.',
                                  file=output_file)
                            print("empty window")
                            possible_update_times.append(w_starttime)
                ixw += 1

            possible_update_times.append(w_starttime + wl - self.cfg.time_overlap)
            # update time
            possible_update_times.sort()
            update_time = next((ut for ut in possible_update_times if ut not in self.updated_times), None)
            print("update time: ", update_time)
            if update_time is None:
                self.t = t_end
                break

            result = self.update_data(update_time)
            if not result:
                self.t = t_end
                break

            # # check if there is a gap
            # while self.t < min([dd.stats.starttime for dd in self.data]):
            #     self.t += self.cfg.time_window_length - self.cfg.time_overlap
            #     print("jumping to t ", self.t)

        # - Write results
        for corr in self._correlations.values():
            corr.write_stack(output_format=self.cfg.format_output)

        print('Finished a correlation block.')

    def perform_checks(self, tr1, tr2, output_file, min_len_samples):
                
        if tr1.stats.starttime - tr2.stats.starttime > tr1.stats.delta:
            print("Traces are not synchronous.", file=output_file)
            return(False)

        if tr1.stats.npts < min_len_samples:
            print("Trace length < min samples\n", file=output_file)
            return(False)

        if tr2.stats.npts < min_len_samples:
            print("Trace length < min samples\n", file=output_file)
            return(False)

        if True in np.isnan(tr1.data):
            print("Trace contains nan\n", file=output_file)
            return(False)

        if True in np.isnan(tr2.data):
            print("Trace contains nan\n", file=output_file)
            return(False)

        if True in np.isinf(tr1.data):
            print("Trace contains inf\n", file=output_file)
            return(False)

        if True in np.isinf(tr2.data):
            print("Trace contains inf\n", file=output_file)
            return(False)

        if len(np.where(tr1.data == 0)) > self.cfg.max_gap_percent_of_window * tr1.stats.npts:
            print("Trace contains more 0 than desired\n", file=output_file)
            return(False)

        if len(np.where(tr2.data == 0)) > self.cfg.max_gap_percent_of_window * tr1.stats.npts:
            print("Trace contains more 0 than desired\n", file=output_file)
            return(False)

        return(True)

    # debugging @profile
    def preprocess(self, tr):

        if True in [self.cfg.cap_glitch, 
                    self.cfg.whiten,
                    self.cfg.onebit,
                    self.cfg.ram_norm] or\
            self.cfg.time_overlap != 0.0:
            tr_temp = Stream()
            for t in tr:
                tr_temp += t.copy()
            tr = tr_temp

        # linear
        if self.cfg.bandpass is not None:
            try:
                w_temp = sosfilt(self.sos, tr.data)
                tr.data = sosfilt(self.sos, w_temp[::-1])[::-1]
            except AttributeError:
                for t in tr:
                    w_temp = sosfilt(self.sos, t.data)
                    t.data = sosfilt(self.sos, w_temp[::-1])[::-1]
        # non-linear
        if self.cfg.cap_glitch:
            cap(tr, self.cfg.cap_thresh)

        if self.cfg.whiten:
            whiten(tr, self.cfg.white_freqmin,
                   self.cfg.white_freqmax,
                   self.cfg.white_taper_samples)

        if self.cfg.onebit:
            for t in tr:
                t.data = np.sign(t.data)

        if self.cfg.ram_norm:
            ram_norm(tr, self.cfg.ram_window, self.cfg.ram_prefilt)

        return(tr)


    # debugging @profile
    def rotate(self, str1, str2, baz1, baz2):
        s_temp1 = str1.copy()
        s_temp2 = str2.copy()

        try:
            s_temp1.rotate('NE->RT', back_azimuth=baz1)
        except:
            if str1[0].stats.station != str2[0].stats.station:
                print('** data not rotated for stream: ')
                print(s_temp1)
            pass

        try:
            s_temp2.rotate('NE->RT', back_azimuth=baz2)
        except:
            if str1[0].stats.station != str2[1].stats.station:
                print('** data not rotated for stream: ')
                print(s_temp2)
            pass
        del str1, str2
        return(s_temp1, s_temp2)

    # debugging @profile
    def update_data(self, update_time):
        # mytracker = tracker.SummaryTracker()
        # mytracker.print_diff()
        # add a new round of data:
        #removallist = []
        #stream_temp = self.data.trim(starttime=self.tself.t - self.cfg.time_window_length -).copy()
        #self.data = Stream()
        #self.data += stream_temp
        #gc.collect()

        for ix_c, channel in enumerate(self.channels):
            while True:
                try:
                    f = self.inv[channel].pop(0)

                    try:
                        newd = read(f)
                        self.data += newd
                        if newd[-1].stats.endtime > update_time + self.cfg.time_window_length * 3:
                            break
                        
                    except IOError:
                        print("** Could not read trace: %s" % f)

                except IndexError:
                    # No more data.
                    break

        # merge. traces with too many zeros are kicked out during quality check
        if len(self.data) > 0:
            self.data.sort(keys=["station"])
            self.data.merge(method=1, fill_value=0.0, interpolation_samples=0)        
            self.data.sort(keys=["starttime"])
            self.data.trim(starttime=update_time - self.cfg.time_window_length * 0.5 * self.cfg.max_gap_percent_of_window,
                pad=True, fill_value=0.0)
            self.t = update_time
            self.updated_times.append(update_time)
            return True
        else:
            return False

    def initialize_data_new(self, t0):
        # t0: begin time of observation
        # have at least one window in any channel
        t_min = t0 + 3 * self.cfg.time_window_length
        
        self.data = Stream()
        for channel in self.channels:
            read_successfully = False
            while not read_successfully:
                f = self.inv[channel].pop(0)
                try:
                    tr = read(f)
                    if tr[-1].stats.endtime > t_min:
                        self.data += tr
                        read_successfully = True
                    else:
                        pass
                except IOError:
                    print('** problems reading file %s'
                          % self.inv.data[channel])
        self.data._cleanup()
        return()

    def get_prepstring(self):

        prepstring = ''
        if self.cfg.bandpass:
            prepstring += 'b'
        else:
            prepstring += '-'
        if self.cfg.cap_glitch:
            prepstring += 'g'
        else:
            prepstring += '-'
        if self.cfg.whiten:
            prepstring += 'w'
        else:
            prepstring += '-'
        if self.cfg.onebit:
            prepstring += 'o'
        else:
            prepstring += '-'
        if self.cfg.ram_norm:
            prepstring += 'r'
        else:
            prepstring += '-'

        return prepstring
