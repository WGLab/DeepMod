
import os;
import sys;
import string;
import glob;
import time
import copy

import h5py
import numpy as np
import multiprocessing

from collections import defaultdict
from distutils.version import LooseVersion

import tempfile
import subprocess

import re;

from . import myCom
#import myCom

import tensorflow as tf
from tensorflow.contrib import rnn
from . import myMultiBiRNN
from . import EventTable
from . import MoveTable

rnn_pred_batch_size = 512

# default path for data in fast5 files
fast5_channel_id= 'UniqueGlobalKey/channel_id'
fast5_analysis = ''.join(['/', myCom.analyses_base]) #
fast5_events = myCom.basecall_events_base #
fast5_rawReads = ''.join(['/', myCom.raw_base, '/', myCom.reads_base]) #
fast5_basecall_fq = myCom.basecall_fastq_base #
fast5_signal = myCom.signal_base #

pre_base_str = 'rnn.pred.ind'

#
# get digitisation, offset, range, sampling_rate from fast5 files
#
def get_channel_info(moptions, sp_param):
   if not sp_param['f5status']=="": return;
   try:
      channel_info = sp_param['f5reader'][fast5_channel_id].attrs
      sp_param["channel_info"] = {'digitisation':channel_info['digitisation'], 'offset':channel_info['offset'], 'range':channel_info['range'], 'sampling_rate':channel_info['sampling_rate'], 'channel_number':channel_info['channel_number']}
   except:
      raiseError("No Channel Info", sp_param, "No Channel Info")

# print error message
def raiseError(sp_info, sp_param, errk):
   sp_param['f5status'] = errk
   print ('Error!!! %s in %s' % (sp_info, sp_param['mfile_path']))
   sys.stdout.flush()

#
# get Albacore version used
# only support v1+ and v2+
#
def getAlbacoreVersion(moptions, sp_param):
   if not sp_param['f5status']=="": return;
   try:
      ver_path = ''.join([fast5_analysis,'/', moptions['basecall_1d'] ])
      #add .decode("utf-8")  to make it compatible to py3
      try:
         used_version = LooseVersion(sp_param['f5reader'][ver_path].attrs['version'].decode("utf-8")  if 'version' in sp_param['f5reader'][ver_path].attrs else "0.0")
      except:
         used_version = LooseVersion(sp_param['f5reader'][ver_path].attrs['version']  if 'version' in sp_param['f5reader'][ver_path].attrs else "0.0")
      sp_param['get_albacore_version'] = used_version
      if used_version < LooseVersion("1.0"): #
         sp_param['used_albacore_version'] = 1;
      elif used_version < LooseVersion("2.0"): sp_param['used_albacore_version'] = 1;
      elif used_version >= LooseVersion("2.0"): sp_param['used_albacore_version'] = 2;
   except: # default verion is 1 now
      sp_param['used_albacore_version'] = 1;

# not used now.
def get_kmer_corrected_info(moptions):
   if ('kmer_model_file' not in moptions) or moptions['kmer_model_file']==None or (not os.path.isfile(moptions['kmer_model_file'])): return;

   fr = open(moptions['kmer_model_file'], 'r')
   moptions['kmer_model_dict'] = defaultdict()
   line = fr.readline();
   while line:
      line = string.strip(line);
      if len(line)>0 and (not line[0]=='#'):
         try:
            c_kmer, c_level_mean, c_level_stdv = line.split()[:3]
            c_level_mean, c_level_stdv = float(c_level_mean), float(c_level_stdv)
            moptions['kmer_model_dict'][c_kmer] = (c_level_mean, 1/(c_level_stdv*c_level_stdv))
         except:
            pass;
      line = fr.readline();
   fr.close();

# not used now
# get shift and scale values for normalization
#
def get_cur_shift_scale(moptions, sp_param):
   if not sp_param['f5status']=="": return;
   if "kmer_model_dict" not in moptions: return;

   event_key = 'm_event'

   try:
      cur_model = np.array([moptions['kmer_model_dict'][c_model_state] for c_model_state in sp_param[event_key]['model_state']], dtype=[('level_mean', np.float), ('level_stdv', np.float)]);
      c_mean_stdv = cur_model['level_mean']*cur_model['level_stdv']
      c_mean_stdv_sum = c_mean_stdv.sum()
      model_coef_matrix = np.array(( (cur_model['level_stdv'].sum(), c_mean_stdv_sum), \
                                     (c_mean_stdv_sum, (c_mean_stdv*cur_model['level_mean']).sum()) \
                                  ))
      c_event_stdv = sp_param[event_key]['mean'] * cur_model['level_stdv']
      c_event_stdv_mean = c_event_stdv * cur_model['level_mean']
      dependent_array = np.array((c_event_stdv.sum(), c_event_stdv_mean.sum()));

      sp_param['shift_scale'] = {}
      sp_param['shift_scale']['cal_shift'], sp_param['shift_scale']['cal_scale'] = np.linalg.solve(model_coef_matrix, dependent_array)
      sp_param['shift_scale']['chn_shift'], sp_param['shift_scale']['chn_scale'] = -sp_param["channel_info"]['offset'], sp_param["channel_info"]['digitisation']/sp_param["channel_info"]['range']

      sp_param['shift_scale']['shift']=sp_param['shift_scale']['chn_shift']+sp_param['shift_scale']['chn_scale']*sp_param['shift_scale']['cal_shift']
      sp_param['shift_scale']['scale']=sp_param['shift_scale']['chn_scale']*sp_param['shift_scale']['cal_scale']

      sp_param['raw_signals'] = np.round(sp_param['raw_signals']/sp_param['shift_scale']['cal_scale'] - sp_param['shift_scale']['cal_shift']/sp_param['shift_scale']['cal_scale'], 6)
   except:
      raiseError('Cannot nanopore correction', sp_param, "Cannot nanopore correction")

#
# get events from a fast5 file
#
def getEvent(moptions, sp_param):
  if not sp_param['f5status']=="": return;

  # If use move tables intead of event tables
  if moptions['move']:
      try: # get events from a fast5 file'
         mv_str = '/'.join(['', 'Analyses', moptions['basecall_1d'], moptions['basecall_2strand'], 'Move'])
         move_data = sp_param['f5reader'][mv_str][()]
         sp_param['events_data'] = move_data
      except:
         raiseError('No move data', sp_param, "No move data")
         return;
      m_event = MoveTable.getMove_Info(moptions, sp_param, move_data)
      sp_param['m_event'] = m_event
      # get sequence from events
      sp_param['m_event_basecall'] = sp_param['fq_seq']
      sp_param['left_right_skip'] = (0, 0)


      return
 # End the part of getting move tables

  try: # get events from a fast5 file
     event_path = ''.join([fast5_analysis, '/', moptions['basecall_1d'], '/', moptions['basecall_2strand'], '/', fast5_events])
     events_data = sp_param['f5reader'][event_path].value
  except:
     raiseError('No events data', sp_param, "No events data")
     return;

  convertError = False;

  if sp_param['f5status'] == "":
     sp_param['events_data'] = events_data
     if sp_param['used_albacore_version']==1:
        move0_left = 0; move0_right = len(events_data)-1;
        while move0_left<move0_right: # get the first non-stay event at the left tail
           if events_data['move'][move0_left]==0: move0_left += 1;
           else: break;
        if move0_left>move0_right-20:
           raiseError(("Too many move0 at 3'(l%d, r%d)" % (move0_left, move0_right)), sp_param, "Remove too many bases on left")
           return;
        while move0_right>move0_left: # get the last non-stay event at the right tail
           if events_data['move'][move0_right]==0: move0_right -= 1
           else: break;
        if move0_right<move0_left+20:
           raiseError(("Too many move0 at 5'(l%d, r%d)" % (move0_left, move0_right)), sp_param, 'Remove too many bases on right')
           return

        # get the starting time
        based_ind = events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"] - sp_param['raw_attributes']['start_time']
        first_base_index_in_raw_signal = np.round(events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"]).astype(np.int64) - sp_param['raw_attributes']['start_time']
        # get the potential error of the starting time
        if first_base_index_in_raw_signal<-2:
           raiseError(('The index of the first base is less than -2(%d=%.6f*%d-%d)' % (first_base_index_in_raw_signal, events_data['start'][move0_left].astype(np.float64), sp_param["channel_info"]["sampling_rate"], sp_param['raw_attributes']['start_time'])), sp_param, "The index of the first base is less than -2")
           return;
        elif first_base_index_in_raw_signal<0:
           first_base_index_in_raw_signal = 0
           if moptions['outLevel']<=myCom.OUTPUT_INFO: print ('Warning!!! first_base_index_in_raw_signal less than 0 ' + sp_param['mfile_path'])

        first_base_index_in_raw_signal = np.uint64(first_base_index_in_raw_signal)

        m_event = []; pre_i = move0_left;
        cur_length=(events_data['length'][pre_i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64');
        for i in range(move0_left+1, move0_right+1):
           if events_data['move'][i]>0: # for non-stay event
              if pre_i==move0_left:
                 m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), first_base_index_in_raw_signal, cur_length, events_data['model_state'][pre_i].upper()))
              else: # calculate starting index in raw signal
                 # calculated position
                 cal_st = (events_data['start'][pre_i]-events_data['start'][move0_left])*sp_param["channel_info"]["sampling_rate"]+based_ind
                 if cal_st<0: print("Wanging Less than 0")
                 if cal_st>0 and cal_st - (m_event[-1][2]+ m_event[-1][3])>0 and (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64')>0:
                    if (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64')>2: #
                        m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64'),  events_data['model_state'][pre_i].upper()))
                        m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), cal_st.astype('uint64'), cur_length, events_data['model_state'][pre_i].upper()))
                    else: # for a normal event
                        m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64')+ cur_length, events_data['model_state'][pre_i].upper()))
                 else:
                    m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], cur_length, events_data['model_state'][pre_i].upper()))
                 if m_event[-1][2]>np.iinfo(np.int64).max-2 or m_event[-1][2]<0:
                    if not convertError:
                       print ('ex: %.7f*%d=%.0f' % (events_data['start'][move0_left].astype(np.float64), sp_param["channel_info"]["sampling_rate"], events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"])), sp_param['raw_attributes']['start_time'], sp_param['mfile_path'], m_event[-1][2], m_event[-1][3]
                    convertError = True;
              pre_i = i;
              cur_length=(events_data['length'][i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64');
           else: # for stay event
              cur_length += (events_data['length'][i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64')
        if sp_param['f5status'] == "": # for the last event
           # calculated position
           cal_st = (events_data['start'][pre_i]-events_data['start'][move0_left])*sp_param["channel_info"]["sampling_rate"]+based_ind
           if cal_st<0: print("Wanging Less than 0")
           if cal_st>0 and cal_st - (m_event[-1][2]+ m_event[-1][3])>0 and (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64')>0:
              if (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64')>2:
                 m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64'),  events_data['model_state'][pre_i].upper()))
                 m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), cal_st.astype('uint64'), cur_length, events_data['model_state'][pre_i].upper()))
              else:
                 m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], (cal_st - (m_event[-1][2]+ m_event[-1][3])).astype('uint64')+ cur_length, events_data['model_state'][pre_i].upper()))
           else:
              m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], cur_length, events_data['model_state'][pre_i].upper()))

        # decode
        m_event = np.array(m_event, dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64), ('model_state', 'U5')]) #'S5')]
        sp_param['m_event'] = m_event
        # get sequence from events
        sp_param['m_event_basecall'] = ''.join([event_model_state[2] for event_model_state in m_event['model_state']]);
        sp_param['left_right_skip'] = (move0_left, len(events_data)-move0_right-1)
     elif sp_param['used_albacore_version']==2:
        if moptions['SignalGroup']=='simple':
           m_event = [];
           pre_i = 0; pre_length = events_data['length'][pre_i].astype('uint64');
           for cur_i in range(1, len(events_data)):
              if events_data['move'][cur_i]>0: # non-stay vents
                 m_event.append( (round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), events_data['start'][pre_i], pre_length, events_data['model_state'][pre_i]) )

                 pre_i = cur_i; pre_length = events_data['length'][pre_i].astype('uint64');
              else: # for stay events
                 pre_length += events_data['length'][cur_i].astype('uint64');
           m_event.append( (round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), events_data['start'][pre_i], pre_length, events_data['model_state'][pre_i]) )
           # format events
           m_event = np.array(m_event, dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64), ('model_state', 'U5')])
        else:
           m_event = EventTable.getEvent_Info(moptions, sp_param, events_data)

        sp_param['m_event'] = m_event
        # get sequence from events
        sp_param['m_event_basecall'] = ''.join([event_model_state[2] for event_model_state in m_event['model_state']]);
        sp_param['left_right_skip'] = (0, 0)
     else:
        raise RuntimeError ("This version of Albacore is not supported. Please use the version of Albacore 1.x or 2.x")

#
# normalize raw signals
#
def mnormalized(moptions, sp_param):

   if not sp_param['m_event']['start'][0] < (sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1]):
      print ('Fatal error signal start position is less than the end position', sp_param['mfile_path'], sp_param['m_event']['start'][0], sp_param['m_event']['start'][-1], sp_param['m_event']['length'][-1])

   # get shift and scale
   mshift = np.median(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])])
   mscale = np.median(np.abs(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])]-mshift));
   # standardize
   sp_param['raw_signals'] = (sp_param['raw_signals'] - mshift)/mscale
   # get meand
   read_med = np.median(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])])
   read_mad = np.median(np.abs(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])] - read_med))
   lower_lim = read_med - (read_mad * 5)
   upper_lim = read_med + (read_mad * 5)
   # normalize as nanoraw did.
   sp_param['raw_signals'] = np.round(np.array([upper_lim if sp_param['raw_signals'][i]>upper_lim else (lower_lim if sp_param['raw_signals'][i]<lower_lim  else sp_param['raw_signals'][i]) for i in range(np.size(sp_param['raw_signals']))]), 3)

#
# get Signal from a fast5 file
#
def getRawInfo(moptions, sp_param):
   if not sp_param['f5status']=="": return;

   try:
      # get attribute of raw signals
      for raw_data in sp_param['f5reader'][fast5_rawReads].values(): pass;
      sp_param['raw_attributes'] = dict(raw_data.attrs.items())

      sp_param['raw_signals'] = raw_data['Signal'][()]
   except:
      raiseError(("No Raw_reads/Signal data %s" % (fast5_rawReads)), sp_param, "No Raw_reads/Signal")

#
# get channel_info, AlbacoreVersion, read_id, Raw Signals, Event from a fast5 file
#
def getFast5Info(moptions, sp_param):
   # get channel info
   get_channel_info(moptions, sp_param)
   if "channel_info" not in sp_param:
      raiseError(("Channel information could not be found in %s " % fast5_channel_id), sp_param, "Channel information could not be found")
      return;
   # get albacore version
   getAlbacoreVersion(moptions, sp_param)
   if 'used_albacore_version' not in sp_param:
      return

   try:
      # get fastq attribute
      fq_path = ''.join([fast5_analysis,'/',moptions['basecall_1d'],'/',moptions['basecall_2strand'],'/',fast5_basecall_fq])
      fq_data = sp_param['f5reader'][fq_path][()]
   except:
      raiseError('No Fastq data', sp_param, "No Fastq data")
   if sp_param['f5status']=="":
      fq_data = (fq_data.decode(encoding="utf-8")).split('\n')
      sp_param['read_id'] = (fq_data[0][1:] if fq_data[0][0]=='@' else fq_data[0]).replace(" ", ":::").replace("\t", "|||")
      sp_param['fq_seq'] = fq_data[1];
   # get raw signals
   getRawInfo(moptions, sp_param)
   # get events
   if sp_param['f5status']=="":
       getEvent(moptions, sp_param)
   # normalize signals.
   if sp_param['f5status']=="":
       mnormalized(moptions, sp_param)

   if sp_param['f5status']=="":
      # get mean, std for each event
      for i in range(len(sp_param['m_event'])):
         if (len(sp_param['raw_signals'][sp_param['m_event']['start'][i]:sp_param['m_event']['start'][i]+sp_param['m_event']['length'][i]])==0):
            print ('Signal out of range {}: {}-{} {};{} for {}'.format(i, sp_param['m_event']['start'][i], sp_param['m_event']['length'][i], len(sp_param['m_event']), len(sp_param['raw_signals']), sp_param['mfile_path']))
            if i>500:
                sp_param['m_event'] = sp_param['m_event'][:i-1]
            else:
                sp_param['f5status']=="Less event"
            break;
         sp_param['m_event']['mean'][i] = round(np.mean(sp_param['raw_signals'][sp_param['m_event']['start'][i]:sp_param['m_event']['start'][i]+sp_param['m_event']['length'][i]]), 3)
         sp_param['m_event']['stdv'][i] = round(np.std(sp_param['raw_signals'][sp_param['m_event']['start'][i]:sp_param['m_event']['start'][i]+sp_param['m_event']['length'][i]]), 3)

#
# associate signals for each event in a fast5 file
#
def get_Event_Signals(moptions, sp_options, f5files):
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      start_time = time.time(); runnum = 0;

   f5data = {}
   sp_options["Error"] = defaultdict(list)
   sp_options["get_albacore_version"] = defaultdict(int)
   # for each fast5 file
   for f5f in f5files:
      try:
         with h5py.File(f5f, 'r') as mf5:
            sp_param = {}
            sp_param['mfile_path'] = f5f
            sp_param['f5reader'] = mf5
            sp_param['f5status'] = "";
            getFast5Info(moptions, sp_param)
            if 'get_albacore_version' in sp_param:
               sp_options["get_albacore_version"][str(sp_param['get_albacore_version'])] += 1
            if sp_param['f5status'] == "":
               if sp_param['read_id'] in f5data:
                  print ('Duplicate id', sp_param['read_id'], f5f)
               f5data[sp_param['read_id']] = (sp_param['m_event_basecall'], sp_param['m_event'], sp_param['raw_signals'], f5f, sp_param['left_right_skip'])
            else:
               sp_options["Error"][sp_param['f5status']].append(f5f)
            # for outputing progress
            if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
               runnum += 1;

               if runnum%500==0:
                  end_time = time.time();
                  print ("%d consuming time %d" % (runnum, end_time-start_time))
      except:
         sp_options["Error"]["Cannot open fast5 or other errors"].append(f5f)
         print("Cannot open fast5 or other errors: {}".format(f5f))
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("All consuming time %d" % (end_time-start_time))

   return f5data;

#
# get signals of events
# map bases from events to a reference genome
#
def mDetect1(moptions, sp_options, f5files):
   f5data = get_Event_Signals(moptions, sp_options, f5files)

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG: start_time = time.time();
   # for fa files of base sequences from events
   temp_fa = tempfile.NamedTemporaryFile(suffix='.fa', mode='w')
   f5keys = sorted(f5data.keys()); #f5keys.sort()
   for f5k in f5keys:
      temp_fa.write(''.join(['>', f5k, '\n', f5data[f5k][0], '\n']))
   temp_fa.flush();
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Write consuming time %d" % (end_time-start_time))

   # run alignmen tools of bwa-mem or minimap2
   temp_sam = tempfile.NamedTemporaryFile()
   if moptions['alignStr']=='bwa':
      cmd_opt = ['mem', '-x', 'ont2d', '-v', '1', '-t', '1', moptions['Ref'], temp_fa.name]
   else:
      cmd_opt = ['-ax', 'map-ont', moptions['Ref'], temp_fa.name]
   returncode = subprocess.call([moptions['alignStr'],]+cmd_opt, stdout=temp_sam)
   if not returncode==0:
      print ('Fatal Error!!! returncode is non-zero(%d) for "%s"' % (returncode, curcmd))
      errkey = "Cannot running aligment"
      for f5k in f5keys:
         sp_options["Error"][errkey].append(f5data[f5k][3])
      return;

   temp_fa.close();
   temp_sam.seek(0);
   # get content from aligned results
   align_info = temp_sam.readlines()
   align_info = [str(align_info[i], 'utf-8').strip() for i in range(len(align_info))]
   temp_sam.close();

   sp_param = defaultdict();
   sp_param['f5data'] = f5data

   f5align = defaultdict()
   f5keydict = defaultdict();
   sp_param['ref_info'] = defaultdict()

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:start_time = time.time();
   ilid = 0;
   # get alignment records
   while ilid < len(align_info):
      if len(align_info[ilid])==0 or align_info[ilid][0]=='@':
         ilid += 1
         continue;

      sp_param['f5status'] = "";
      sp_param['line'] = align_info[ilid]
      qname = handle_line(moptions, sp_param, f5align)
      if sp_param['f5status'] == "":
         f5keydict[qname] = True;
      ilid += 1

   # get unmapped reads
   for f5k in f5keys:
      if f5k not in f5keydict:
         sp_options["Error"]["Not in alignment sam"].append(f5data[f5k][3])

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Get BAM consuming time %d" % (end_time-start_time))

   # get features, prediction for each fast5 files
   sp_param['f5status']= ""
   sp_param['line'] = ""
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:start_time = time.time();
   handle_record(moptions, sp_options, sp_param, f5align, f5data)
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Analyze & annotate & save consuming time %d" % (end_time-start_time))

#
# get reference sequenceng from a reference genome
#
def getRefSeq(moptions, sp_param, rname):
   temp_seq = tempfile.NamedTemporaryFile()
   # command to be run
   cmd_opt = ['faidx', moptions['Ref'], rname]
   returncode = subprocess.call(['samtools',]+cmd_opt, stdout=temp_seq)
   if not returncode==0:
      print ('Fatal Error!!! cannot find the chrosome sequence %s' % rname)
   else:
      temp_seq.seek(0);
      seqinfo = temp_seq.readlines()
      seqinfo = [str(seqinfo[i], 'utf-8').strip() for i in range(len(seqinfo))]
      temp_seq.close();

      sp_param['ref_info'][rname] = ''.join(seqinfo[1:]).strip().upper()

#
# get mapping information and associate it with events/signals.
#
def handle_record(moptions, sp_options, sp_param, f5align, f5data):
   alignkeys = list(f5align.keys());
   numreg = re.compile('\d+')
   mdireg = re.compile('[MIDNSHPX=]{1}')

   # for each alignment record
   for readk_ind in range(len(alignkeys)):
     sp_param['f5status']= ""
     readk = alignkeys[readk_ind]
     sp_param['mfile_path'] = f5data[readk][3]

     mapq, flag, rname, pos, cigar, readseq = f5align[readk]

     # output chromosome of interest
     if (not moptions['ConUnk']) and ((not rname.find('_')==-1) or (not rname.find('-')==-1) or (not rname.find('/')==-1) or (not rname.find(':')==-1)):
        continue;
     isinreg = False;
     # check the region of interest
     for cur_mr in moptions['region']:
        if (cur_mr[0] in ['', None, rname]):
           isinreg = True;
           break;
     if not isinreg:
        continue;

     # get reference information
     if rname not in sp_param['ref_info']:
        getRefSeq(moptions, sp_param, rname)
     refseq = sp_param['ref_info'][rname]

     # mapped position and strand
     pos = pos - 1
     forward_reverse = '-' if flag&0x10 else '+'

     numinfo = numreg.findall(cigar);
     mdiinfo = mdireg.findall(cigar)
     numinfo = [int(numinfo[i]) for i in range(len(numinfo))] #map(int, numinfo)

     leftclip = 0; rightclip = 0;
     # remove left clip
     while mdiinfo[0] in ['I', 'D', 'N', 'S', 'H', 'P', 'X']:
         if mdiinfo[0] in ['I', 'S', 'X']:
            leftclip += numinfo[0];  readseq = readseq[numinfo[0]:]
         if mdiinfo[0] in ['H']: leftclip += numinfo[0]
         if mdiinfo[0] in ['D', 'N', 'X']:
            pos += numinfo[0]
         numinfo = numinfo[1:];  mdiinfo = mdiinfo[1:]
     # remove right clip
     while mdiinfo[-1] in ['I', 'D', 'N', 'S', 'H', 'P', 'X']:
         if mdiinfo[-1] in ['I', 'S', 'X']:
            rightclip += numinfo[-1]; readseq = readseq[:-numinfo[-1]]
         if mdiinfo[-1] in ['H']: rightclip += numinfo[-1]
         numinfo = numinfo[:-1]; mdiinfo = mdiinfo[:-1]
     if forward_reverse=='+': # remove clipped events
        if rightclip>0: m_event = f5data[readk][1][leftclip:-rightclip]
        else: m_event = f5data[readk][1][leftclip:]
     else:
        if leftclip>0: m_event = f5data[readk][1][rightclip:-leftclip]
        else: m_event = f5data[readk][1][rightclip:]

     isinreg = False;
     # for specify regions
     for cur_mr in moptions['region']:
        if (cur_mr[0] in ['', None, rname]) and \
           (cur_mr[1] in ['', None] or pos>cur_mr[1]) and \
           (cur_mr[2] in ['', None] or pos+len(m_event)<cur_mr[2]):
           isinreg = True;
           break;
     ### for check consistency
     if not isinreg:
        continue;

     lastmatch = None; firstmatch = None;
     first_match_pos = None; last_match_pos = None
     last_al_match = None;  first_al_match = None
     lasmtind = 0;
     base_map_info = []; #indel_groups = defaultdict()
     nummismatch = 0; numinsert = 0; numdel = 0;
     read_ind = 0;
     # get map detail: reference base, read base, ref position, read-position
     for n1ind in range(len(numinfo)):
        mdi = mdiinfo[n1ind];
        for n1i in range(numinfo[n1ind]):
           # for each cigar type
           if mdi=='M':
              base_map_info.append((refseq[pos], readseq[read_ind], pos, read_ind, 0))
              if refseq[pos]==readseq[read_ind]:
                 if firstmatch==None: firstmatch = read_ind
                 if lastmatch==None or lastmatch<read_ind: lastmatch = read_ind; lasmtind=n1ind
                 if first_al_match==None: first_al_match=len(base_map_info)-1
                 if last_al_match==None or last_al_match<len(base_map_info): last_al_match=len(base_map_info)-1
                 if first_match_pos==None: first_match_pos = pos
                 if last_match_pos==None or last_match_pos<pos: last_match_pos = pos
              else: nummismatch += 1
              pos += 1; read_ind += 1;
           elif mdi =='I':
              base_map_info.append(('-', readseq[read_ind], pos, read_ind, 0))
              read_ind += 1;
              numinsert += 1
           elif mdi == 'D':
              base_map_info.append((refseq[pos], '-', pos, read_ind, 0))
              pos += 1;
              numdel += 1
           elif mdi == 'N':
              base_map_info.append((refseq[pos], '-', pos, read_ind, 0))
              pos += 1;
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print ('CIGAR-Error N exist', f5data[readk][3])
           elif mdi == 'S':
              read_ind += 1;
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print ('CIGAR-Error!!! S in the middle of the sequence', f5data[readk][3])
           elif mdi == 'H':
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print ('CIGAR-Error!!! H in the middle of the sequence', f5data[readk][3])
           elif mdi == 'P':
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print ('CIGAR-Error!!! P exist', f5data[readk][3])
           elif mdi == '=':
             base_map_info.append((refseq[pos], readseq[read_ind], pos, read_ind, 0))
             if first_match_pos==None: first_match_pos  = pos
             if last_match_pos==None or last_match_pos<pos: last_match_pos = pos
             pos += 1; read_ind += 1;
             if firstmatch==None: firstmatch = read_ind - 1
             if lastmatch==None or lastmatch<read_ind-1: lastmatch = read_ind - 1; lasmtind=n1ind
             if last_al_match==None or last_al_match<len(base_map_info): last_al_match=len(base_map_info)-1
             if first_al_match==None: first_al_match=len(base_map_info)-1
           elif mdi == 'X':
             base_map_info.append((refseq[pos], readseq[read_ind], pos, read_ind, 0))
             pos += 1; read_ind += 1;
             nummismatch += 1
           else:
             if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                print ('CIGAR-Error!!!', 'Warning unknow CIGAR element ' + str(numinfo[n1ind]) + ' ' + mdi, f5data[readk][3])
     if firstmatch==None or lastmatch==None or firstmatch<0 or lastmatch<0:
        if moptions['outLevel']<=myCom.OUTPUT_WARNING:
           print ("Errorfast5 "+f5data[readk][3])
           print('match-Error!!! no first and/or last match',f5data[readk][3],('firstmatch=%d' % firstmatch) if not (firstmatch==None) else "N", ('lastmatch%d' % lastmatch) if not (lastmatch==None) else "N", str(flag), rname, str(pos));
           print('\tf=%d, chr=%s, p=%d, c=%s, s=%s' % (flag, rname, pos, cigar, readseq))
           continue;

     # remove more un-matched events
     if forward_reverse=='+':
        if not firstmatch==None: leftclip += firstmatch
        if (not lastmatch==None) and len(m_event)-lastmatch>1: rightclip += len(m_event)-lastmatch-1
     else:
        if not firstmatch==None: rightclip += firstmatch
        if (not lastmatch==None) and len(m_event)-lastmatch>1: leftclip += len(m_event)-lastmatch-1

     if forward_reverse=='+':
        if len(m_event)-lastmatch>1:
           m_event = m_event[firstmatch:(lastmatch+1-len(m_event))]
        elif firstmatch>0: m_event = m_event[firstmatch:]
     else:
        if firstmatch>0: m_event = m_event[(len(m_event)-1-lastmatch):-firstmatch]
        elif len(m_event)-lastmatch>1: m_event = m_event[(len(m_event)-1-lastmatch):]
     # check potential error
     if firstmatch>0 or len(base_map_info)-last_al_match>1:
        if moptions['outLevel']<=myCom.OUTPUT_WARNING and ((firstmatch>0) or (len(base_map_info)-last_al_match>1 and refseq[last_match_pos+1] not in ['N'])):
           print ("Errorfast5"+f5data[readk][3])
           print ('Warning!!! first not match', firstmatch, lastmatch, first_al_match, last_al_match, len(base_map_info), numinfo[lasmtind-2:(lasmtind+5)], mdiinfo[lasmtind-2:(lasmtind+5)], lasmtind, len(numinfo))
           print('\tref='+refseq[last_match_pos:last_match_pos+20]+"\n\tred="+readseq[lastmatch:lastmatch+20])
           if firstmatch>0:
              print('\tref='+refseq[(first_match_pos-20 if first_match_pos-20>0 else 0):first_match_pos]+"\n\tred="+readseq[(firstmatch-20 if firstmatch-20>0 else 0):firstmatch])
           print('\tf=%d, chr=%s, p=%d, c=%s, s=%s' % (flag, rname, pos, cigar, readseq)) # flag, rname, pos, cigar, readseq

        if len(base_map_info)-last_al_match>1:
           base_map_info = base_map_info[first_al_match:(last_al_match+1-len(base_map_info))]
        elif first_al_match>0:
           base_map_info = base_map_info[first_al_match:]

     # format base
     base_map_info = np.array(base_map_info, dtype=[('refbase', 'U1'), ('readbase', 'U1'), ('refbasei', np.uint64), ('readbasei', np.uint64), ('mod_pred', np.int)])
     if forward_reverse=='-':
        base_map_info = np.flipud(base_map_info)
        for bmii in range(len(base_map_info)):
            base_map_info['refbase'][bmii]  = get_complement(base_map_info['refbase'][bmii])
            base_map_info['readbase'][bmii] = get_complement(base_map_info['readbase'][bmii])
        leftclip, rightclip = rightclip, leftclip
     if False: #True: # for test base_map_info  ### for check consistency
        ref_align_key = '/Analyses/NanomoCorrected_000/BaseCalled_template/Alignment/genome_alignment'
        read_align_key = '/Analyses/NanomoCorrected_000/BaseCalled_template/Alignment/read_alignment'
        with h5py.File(f5data[readk][3], 'r') as mf5:
           read_align_list = [bt.decode(encoding="utf-8") for bt in mf5[read_align_key]]
           ref_align_list = [bt.decode(encoding="utf-8") for bt in mf5[ref_align_key]]
           for rali in range(len(read_align_list)):
              if not read_align_list[rali]==base_map_info['readbase'][rali]:
                 print ("Error not equal1! %s %s %d %s" % (read_align_list[rali], base_map_info['readbase'][rali], rali, f5data[readk][3]))
              if not ref_align_list[rali]==base_map_info['refbase'][rali]:
                 print ("Error not equal2! %s %s %d %s" % (ref_align_list[rali], base_map_info['refbase'][rali], rali, f5data[readk][3]))

        # ## for check consistency
     if True:
        # for special alignment
        # CCG     or CGG
        # C-G        C-G
        for ali in range(len(base_map_info)):
           if base_map_info['refbase'][ali]=='C' and base_map_info['readbase'][ali]=='C':
              if ali+1<len(base_map_info) and base_map_info['readbase'][ali+1]=='-' and base_map_info['refbase'][ali+1]=='G':
                 addali = 2;
                 while ali + addali < len(base_map_info):
                     if base_map_info['readbase'][ali+addali]=='-' and base_map_info['refbase'][ali+addali]=='G': addali += 1;
                     else: break;
                 if ali + addali < len(base_map_info) and base_map_info['readbase'][ali+addali]=='G' and base_map_info['refbase'][ali+addali]=='G':
                    base_map_info['readbase'][ali+1], base_map_info['readbase'][ali+addali] = base_map_info['readbase'][ali+addali], base_map_info['readbase'][ali+1]
           if base_map_info['refbase'][ali]=='G' and base_map_info['readbase'][ali]=='G':
              if ali-1>-1 and base_map_info['readbase'][ali-1]=='-' and base_map_info['refbase'][ali-1]=='C':
                 addali = 2;
                 while ali - addali >-1:
                     if base_map_info['readbase'][ali-addali]=='-' and base_map_info['refbase'][ali-addali]=='C': addali += 1;
                     else: break;
                 if ali - addali>-1 and base_map_info['readbase'][ali-addali]=='C' and base_map_info['refbase'][ali-addali]=='C':
                     base_map_info['readbase'][ali-1], base_map_info['readbase'][ali-addali] = base_map_info['readbase'][ali-addali], base_map_info['readbase'][ali-1]

     if len(m_event)<50:
         raiseError("Less Event", sp_param, "Less Event");
         sp_options["Error"]["Less Event"].append(f5data[readk][3])
         continue;

     # get feature
     mfeatures,isdif = get_Feature(moptions, sp_options, sp_param, f5align, f5data, readk, leftclip, rightclip, base_map_info, forward_reverse, rname, first_match_pos, numinsert, numdel)
     if isdif and moptions['outLevel']<=myCom.OUTPUT_WARNING:
        print("Dif is true")
        print([lastmatch, firstmatch, first_match_pos, last_match_pos, first_al_match, last_al_match, lasmtind, len(base_map_info), nummismatch, numinsert, numdel, len(base_map_info)-nummismatch-numinsert-numdel])
     if not sp_param['f5status']=="": continue

     # generate/save prediction information
     pred_mod_num = mPredict1(moptions, sp_options, sp_param, mfeatures, base_map_info, readk, leftclip, rightclip)
     predfile = (sp_options['ctfolder'] if sp_options['ctfolder'][-1] not in ['/', '\\'] else sp_options['ctfolder'][:-1])+'/rnn.pred.detail.fast5'+'.'+str(sp_options['batchid'])
     pred_f5_key = 'pred_'+str(readk_ind)
     sp_options['Mod'].append([rname, forward_reverse, f5align[readk][3]-1, pred_f5_key, f5data[readk][3][len(moptions['wrkBase'])+1:], predfile[len(moptions['outFolder']+moptions['FileID'])+1:]])
     # save the prediction
     with h5py.File(predfile, 'a') as save_data:
         if "pred" in save_data: base_group = save_data["pred"]
         else: base_group = save_data.create_group("pred")

         if pred_f5_key in base_group:
            del base_group[pred_f5_key]
            save_data.flush()
         pred_group = base_group.create_group(pred_f5_key)

         # save mapped chr, strand, positions
         pred_group.attrs['mapped_chr'] = rname
         pred_group.attrs['mapped_strand'] = forward_reverse
         pred_group.attrs['mapped_start'] = base_map_info['refbasei'][0] if forward_reverse=='+' else base_map_info['refbasei'][-1]
         pred_group.attrs['mapped_end'] = base_map_info['refbasei'][-1] if forward_reverse=='+' else base_map_info['refbasei'][0]

         if forward_reverse=='+':
            pred_group.attrs['clipped_bases_start'] = leftclip
            pred_group.attrs['clipped_bases_end'] = rightclip
         else:
            pred_group.attrs['clipped_bases_start'] = rightclip
            pred_group.attrs['clipped_bases_end'] = leftclip

         # save indel, matches/mismatches
         pred_group.attrs['num_insertions'] = numinsert
         pred_group.attrs['num_deletions'] = numdel
         pred_group.attrs['num_matches'] = len(base_map_info)-nummismatch-numinsert-numdel
         pred_group.attrs['num_mismatches'] = nummismatch

         # save the fast5 file and prediction
         pred_group.attrs['pred_mod_num'] = pred_mod_num
         pred_group.attrs['f5file'] = f5data[readk][3]
         pred_group.attrs['readk'] = readk
         base_map_info = np.array(base_map_info, dtype=[('refbase', 'S1'), ('readbase', 'S1'), ('refbasei', np.uint64), ('readbasei', np.uint64), ('mod_pred', np.int)])
         pred_group.create_dataset('predetail', data=base_map_info, compression="gzip")

         try:
            save_data.flush();
            save_data.close();
         except:
            sp_options["Error"]['Cannot save data'].append(f5data[readk][3])
            print ('Error!!! %s in %s' % ("Cannot save data", f5data[readk][3]))

   # save index information
   sp_options['Mod'] = sorted(sp_options['Mod'])
   # index file
   pred_ind_file =  (sp_options['ctfolder'] if sp_options['ctfolder'][-1] not in ['/', '\\'] else sp_options['ctfolder'][:-1])+'/%s.' + pre_base_str + '.' + str(sp_options['batchid'])
   if len(sp_options['Mod'])>0:
      cur_chr = None; cur_writer = None;
      for mfi in sp_options['Mod']:
         if cur_chr==None or (not cur_chr == mfi[0]):
            if not cur_chr==None:
               cur_writer.flush();
               cur_writer.close()
            cur_chr = mfi[0]
            cur_writer = open((pred_ind_file % cur_chr), 'w')
         cur_m_f = []
         for mfidetail in mfi:
            cur_m_f.append(str(mfidetail))
         cur_m_f.append('\n')
         cur_writer.write(' '.join(cur_m_f))
      if not cur_writer==None:
         cur_writer.flush();
         cur_writer.close()

#
# make modificatoin prediction for a long read
#
def mPredict1(moptions, sp_options, sp_param, mfeatures, base_map_info, readk, start_clip, end_clip):
   #
   modevents = sp_param['f5data'][readk][1]
   # get features. labels might be all zero
   t0, ty, tx = np.split(mfeatures, [1,3], axis=1);
   t0 = t0.astype(int)
   m_data = []; m_y = [];
   for ie in range(start_clip-100, len(modevents)-end_clip+100):
      mind = ie - (start_clip-100)
      if ie>=start_clip and ie<len(modevents)-end_clip:
          m_y.append(ty[mind])
          # format to input with windoe size
          m_data.append(tx[(mind-int(moptions['windowsize']/2)):(mind+int(moptions['windowsize']/2)+1)])

   # for input
   test_feature = np.reshape(m_data, (len(m_data), len(m_data[0]), len(m_data[0][0])))
   test_label = np.reshape(m_y, (len(m_y), len(m_y[0]))).astype(int)

   sp_options['rnn'][0].run(sp_options['rnn'][3])

   # split into small group
   if len(test_feature) > rnn_pred_batch_size*1.2:
      x_sub_group = np.array_split(test_feature, int(len(test_feature)/rnn_pred_batch_size))
      y_sub_group = np.array_split(test_label,   int(len(test_feature)/rnn_pred_batch_size))
   else:
      x_sub_group = [test_feature]; y_sub_group = [test_label]
   # make prediction on each small groups
   for subi in range(len(x_sub_group)):
      if subi==0:
         mfpred_output = (sp_options['rnn'][0].run([sp_options['rnn'][4]], \
               feed_dict={sp_options['rnn'][1]:x_sub_group[subi], sp_options['rnn'][2]:y_sub_group[subi]}))[0];
      else:
         mfpred_output = np.concatenate((mfpred_output, (sp_options['rnn'][0].run([sp_options['rnn'][4]], \
               feed_dict={sp_options['rnn'][1]:x_sub_group[subi], sp_options['rnn'][2]:y_sub_group[subi]}))[0]), axis=0);

   # associate the prediction with reference positions and read positions
   modevents = sp_param['f5data'][readk][1]
   aligni = 0; pred_mod_num = 0;
   for ie in range(start_clip, len(modevents)-end_clip):
      while base_map_info['readbase'][aligni]=='-': aligni += 1
      if not base_map_info['readbase'][aligni] == modevents['model_state'][ie][2]:
         print ('Error Does not match', base_map_info['readbase'][aligni], modevents['model_state'][ie][2], aligni, ie)
      if mfpred_output[ie-start_clip]==1:
         base_map_info['mod_pred'][aligni] = 1;
         pred_mod_num += 1;

      aligni += 1
   return pred_mod_num

#
# get feature for a long read
#
def get_Feature(moptions, sp_options, sp_param, f5align, f5data, readk, start_clip, end_clip, base_map_info, forward_reverse, rname, mapped_start_pos, num_insertions, num_deletions):
   modevents = sp_param['f5data'][readk][1]
   # class num, bin num, and bin size
   clnum = 2; binnum = 50; binlen = 0.2;
   if forward_reverse=='+':
      align_ref_pos = mapped_start_pos
   else:
      align_ref_pos = mapped_start_pos + len(base_map_info) - num_insertions - 1

   # initialize feature matrix
   if moptions['fnum']==57:
      mfeatures = np.zeros((len(modevents)-end_clip+100-(start_clip-100), (binnum+3+3+4)));
   else: mfeatures = np.zeros((len(modevents)-end_clip+100-(start_clip-100), (3+3+4)));

   aligni = 0; isdif = False;
   # get feature for each event; each event is a row
   for ie in range(start_clip-100, len(modevents)-end_clip+100):
      cur_row_num = ie - (start_clip-100); cur_base = ''
      if ie>=start_clip and ie<len(modevents)-end_clip:
         if align_ref_pos<mapped_start_pos:
            print ('ERRRR align_ref_pos(%d)<mapped_start_pos(%d)' % (align_ref_pos, mapped_start_pos))
         # get non-indel events
         while base_map_info['readbase'][aligni]=='-':
            if not align_ref_pos==base_map_info['refbasei'][aligni]:
               print ('ERRRR align_ref_pos(%d) not equal to %d' % (align_ref_pos, base_map_info['refbasei'][aligni] ))
            if not base_map_info['refbase'][aligni]=='-':
               if forward_reverse=='+': align_ref_pos += 1
               else: align_ref_pos -= 1
            aligni += 1
         if not base_map_info['readbase'][aligni] == modevents['model_state'][ie][2]:
            print ('Error Does not match', base_map_info['readbase'][aligni], modevents['model_state'][ie][2], aligni, ie)
            sp_param['f5status']= "Error Does not match"
            if f5data[readk][3] not in sp_options["Error"]['Error Does not match']:
               sp_options["Error"]['Error Does not match'].append(f5data[readk][3])
            if aligni>50: break;
            isdif = True;
         mfeatures[cur_row_num][0] = align_ref_pos
         cur_base = base_map_info['refbase'][aligni]
         # both positive/negative labels is zero
         mfeatures[cur_row_num][1] = 0; mfeatures[cur_row_num][2] = 0
         if not base_map_info['refbase'][aligni]=='-':
            if forward_reverse=='+': align_ref_pos += 1
            else: align_ref_pos -= 1
         aligni += 1

      # for bin features
      if ie>=0 and ie<len(modevents) and moptions['fnum']==57:
         for currs in sp_param['f5data'][readk][2][modevents['start'][ie]:int(modevents['start'][ie]+int(modevents['length'][ie]+0.5))]:
             if currs>10 or currs<-10: print ('Error raw signal', currs, ie, modevents['start'][ie], modevents['length'][ie])
             curbin = int((currs+5)/binlen)
             if curbin<0: curbin = 0
             elif not curbin<binnum: curbin = binnum-1
             mfeatures[cur_row_num][curbin+3] += 1
      if ie>=0 and ie<len(modevents):
         # for reference base type
         if cur_base in myCom.g_ACGT:
            mfeatures[cur_row_num][moptions['fnum']-3+3-4+myCom.g_ACGT.index(cur_base)] = 1
         cur_index_add = moptions['fnum'] - 3 + 3
         # for mean, std and length of signals
         mfeatures[cur_row_num][cur_index_add + 0] = modevents["mean"][ie]
         mfeatures[cur_row_num][cur_index_add + 1] = modevents["stdv"][ie]
         mfeatures[cur_row_num][cur_index_add + 2] = modevents["length"][ie]


   return (mfeatures, isdif)


#
# get complementary base of a given base
#
def get_complement(na):
   if na in myCom.acgt: return myCom.na_bp[na]
   else: return na;

#
# get mean/std of signals
#
def calculate_mean_std(m_event, event_ind, forward_reverse, raw_pv, moptions, sp_param):
   if forward_reverse=='-':
      pvsignal = raw_pv[m_event[-event_ind-1][2]:(m_event[-event_ind-1][2]+m_event[-event_ind-1][3])]
   else:
      pvsignal = raw_pv[m_event[event_ind][2]:(m_event[event_ind][2]+m_event[event_ind][3])]
   # get mean/std for an event
   c_mean = round(np.mean(pvsignal), 3)
   c_std = round(np.std(pvsignal), 3)
   return (c_mean, c_std)

#
# get required information from a mapped record
#
def handle_line(moptions, sp_param, f5align):
   lsp = sp_param['line'].split('\t')
   qname, flag, rname, pos, mapq, cigar, _, _, _, seq, _ = lsp[:11]
   # check query name, map quality, reference position, cigar and reference name
   if qname=='*': sp_param['f5status'] = "qname is *"
   elif int(mapq)==255: sp_param['f5status'] = "mapq is 255"
   elif int(pos)==0: sp_param['f5status'] = "pos is 0"
   elif cigar=='*': sp_param['f5status'] = "cigar is *"
   elif rname=='*': sp_param['f5status'] = "rname is *"
   if not sp_param['f5status']=="": return qname

   if (qname not in f5align) or f5align[qname][0]<int(mapq):
      f5align[qname] = (int(mapq), int(flag), rname, int(pos), cigar, seq)

   return qname

#
# the worker of the detection in a multiprocessing way
#
def detect_handler(moptions, h5files_Q, failed_Q, file_map_info_q):

   _, init_l, _, _, _, X, Y, _, _, _, _, mfpred = myMultiBiRNN.mCreateSession(moptions['fnum'], moptions['hidden'], moptions['windowsize'], moptions)
   config = tf.ConfigProto()
   config.gpu_options.allow_growth = True
   sess = tf.Session(config=config)
   # load module
   new_saver = tf.train.import_meta_graph(moptions['modfile'][0]+'.meta')
   new_saver.restore(sess,tf.train.latest_checkpoint(moptions['modfile'][1]))

   while not h5files_Q.empty():
      cur_start_time = time.time()
      try:
         # get fast5 file
         f5files, ctfolderid, batchid = h5files_Q.get(block=False)
      except:
         break;

      sp_options = defaultdict();
      sp_options['ctfolderid'] = ctfolderid
      sp_options['ctfolder'] = moptions['outFolder']+moptions['FileID']+'/'+str(ctfolderid)
      if not os.path.isdir(sp_options['ctfolder']):
         os.system('mkdir '+sp_options['ctfolder'])
      #if moptions['testrnn']:
      sp_options['rnn'] = (sess, X, Y, init_l, mfpred)
      sp_options['batchid'] = batchid

      sp_options['Mod'] = [];
      # make modification prediction for each fast5
      mDetect1(moptions, sp_options, f5files)
      # outputing errors
      for errtype, errfiles in sp_options["Error"].items():
         failed_Q.put((errtype, errfiles));

      print ("Cur Prediction consuming time %d for %d %d" % (time.time() - cur_start_time, ctfolderid, batchid))

   sess.close()

#
# read index files for each output prediction file
#
def read_file_list(cur_cif, cur_chr, cur_strand, sp_options):
   cur_list = []
   with open(cur_cif, 'r') as mr:
       line = mr.readline();
       while line: # get where is prediction file from each line
           line = line.strip();
           if len(line)>0:
              lsp = line.split();
              if line[0]=='#':
                 if lsp[1][0] not in ['/', '\\']:
                    lsp[1] = lsp[1] + '/'
                 if lsp[0]=='#base_folder_fast5': sp_options['base_folder_fast5'] = lsp[1];
                 elif lsp[0]=='#base_folder_output': sp_options['base_folder_output'] = lsp[1];
              else:
                 if lsp[1] == cur_strand:
                    cur_list.append(lsp)
                 if not lsp[0]==cur_chr:
                    print ('Warning!!! The chr should be %s but % is found.' % (cur_chr, lsp[0]))
           line = mr.readline();
   sp_options['handlingList'] = cur_list

#
# get prediction detail
#
def read_pred_detail(moptions, sp_options, f5info):
   f5pred_file = sp_options['base_folder_output'] + '/' + f5info[5]
   f5_pred_key = ('/pred/%s/predetail' % f5info[3])
   # get prediction detail from saved prediction file
   # each file contains predictions for multiple fast5
   with h5py.File(f5pred_file, 'r') as mr:
      m_pred = mr[f5_pred_key].value;
      mapped_chrom = mr['/pred/%s' % f5info[3]].attrs['mapped_chr'] #.decode(encoding="utf-8")
      mapped_strand= mr['/pred/%s' % f5info[3]].attrs['mapped_strand'] #.decode(encoding="utf-8")
   m_pred = np.array(m_pred, dtype=[('refbase', 'U1'), ('readbase', 'U1'), ('refbasei', np.uint64), ('readbasei', np.uint64), ('mod_pred', np.int)])
   return (m_pred, mapped_chrom, mapped_strand)

#
# summarize modification for each genome position of interest
#
def sum_handler(moptions, chr_strand_Q):
   while not chr_strand_Q.empty():
      try:
          # get setting for summarization of predicted modifications
          cur_cif, cur_chr, cur_strand = chr_strand_Q.get(block=False)
      except:
         break;

      sp_options = {}
      # get prediction files
      read_file_list(cur_cif, cur_chr, cur_strand, sp_options)
      sp_options['4NA'] = {moptions['Base']:defaultdict()}
      sp_options['4NAfile'] = {}
      for nak in sp_options['4NA']:
         if not moptions['mod_cluster']:
            sp_options['4NAfile'][nak] = ('%s/mod_pos.%s%s.%s.bed' % (moptions['outFolder'], cur_chr, cur_strand, nak))
         else:
            sp_options['4NAfile'][nak] = ('%s/cluster_mod_pos.%s%s.%s.bed' % (moptions['outFolder'], cur_chr, cur_strand, nak))

      cur_start_time = time.time(); hlnum = 0;
      for hl in sp_options['handlingList']:
         # read prediction detail for each fast5
         m_pred, mapped_chrom, mapped_strand = read_pred_detail(moptions, sp_options, hl)
         if not (mapped_chrom==cur_chr and mapped_strand==cur_strand):
            print("ERRoR not the same chr (real=%s vs expect=%s) and strand (real=%s VS expect=%s)" % (mapped_chrom, cur_chr, mapped_strand, cur_strand))
         #####################################################
         if moptions['mod_cluster']:  # revised; should not used now
            from numpy.lib.recfunctions import append_fields
            m_pred = append_fields(m_pred, 'mod_pred2', m_pred['mod_pred']+0, usemask=False)
            for mi in range(len(m_pred)):
               if not m_pred['mod_pred2'][mi] == m_pred['mod_pred'][mi]:
                  print('mod_pred_Error %d %d %d' % (mi, m_pred['mod_pred2'][mi], m_pred['mod_pred'][mi]))
               if m_pred['mod_pred2'][mi]==1 or m_pred['refbase'][mi] not in ['C']: continue;

               m_3 = []
               m_5 = []
               mj = mi-1
               while mj>-1 and len(m_3)<12:
                  if m_pred['refbase'][mj] in ['N', 'n']: break;
                  if m_pred['refbase'][mj] not in ['-']:
                     m_3.append((m_pred['refbase'][mj], m_pred['mod_pred2'][mj]))
                  mj -= 1;
               if len(m_3)>0: m_3 = m_3[::-1]
               mj = mi + 1;
               while mj < len(m_pred) and len(m_5)<12:
                  if m_pred['refbase'][mj] in ['N', 'n']: break;
                  if m_pred['refbase'][mj] not in ['-']:
                     m_5.append((m_pred['refbase'][mj], m_pred['mod_pred2'][mj]))
                  mj += 1
               cpgnum = 0; meth_cpgnum = 0;
               m_3and5 = [m_3, m_5]
               for m_53 in m_3and5:
                  for mj in range(len(m_53)-1):
                     if m_53[mj][0]=='C' and m_53[mj+1][0]=='G':
                        cpgnum +=1
                        if -0.1 < m_53[mj][1]-1 < 0.1:
                           meth_cpgnum += 1

               if cpgnum>0 and (meth_cpgnum>0 and meth_cpgnum/float(cpgnum)>0.5):
                  m_pred['mod_pred'][mi] = 1
         #####################################################################################
         for mi in range(len(m_pred)):
            # get prediction for each base type
            if m_pred['refbase'][mi] not in sp_options['4NA']: continue;
            if m_pred['refbase'][mi] in ['-', 'N', 'n']: continue;
            if (cur_chr, cur_strand, m_pred['refbasei'][mi]) not in sp_options['4NA'][m_pred['refbase'][mi]]:
               sp_options['4NA'][m_pred['refbase'][mi]][(cur_chr, cur_strand, int(m_pred['refbasei'][mi]) )] = [0, 0, m_pred['refbase'][mi]]
            if not (m_pred['refbase'][mi] == sp_options['4NA'][m_pred['refbase'][mi]][(cur_chr, cur_strand, int(m_pred['refbasei'][mi]) )][2]):
               print ('Error !!!! NA not equal %s == %s' % (m_pred['refbase'][mi], sp_options['4NA'][m_pred['refbase'][mi]][(cur_chr, cur_strand, int(m_pred['refbasei'][mi]) )][2]))
            if not m_pred['readbase'][mi]=='-':
               sp_options['4NA'][m_pred['refbase'][mi]][(cur_chr, cur_strand, int(m_pred['refbasei'][mi]) )][0] += 1
               if -0.1 < m_pred['mod_pred'][mi]-1 < 0.1:
                  sp_options['4NA'][m_pred['refbase'][mi]][(cur_chr, cur_strand, int(m_pred['refbasei'][mi]) )][1] += 1
         hlnum += 1
         if hlnum % 1000==0:
            print ("\tCurrent time consuming %d for %d" % (time.time() - cur_start_time, hlnum))
            cur_start_time = time.time()

      print ('====sum done! To save')
      for nak in sp_options['4NA']:
         print ('\tSave %s' % sp_options['4NAfile'][nak])
         if len(sp_options['4NA'][nak])>0:
            with open(sp_options['4NAfile'][nak], 'w') as mw:
                # save the summary for each base type
                pos_keys = sp_options['4NA'][nak].keys();
                pos_keys = sorted(pos_keys);
                for pk in pos_keys:
                    neighna = [sp_options['4NA'][nak][pk][2]]
                    mw.write(' '.join([ pk[0], str(pk[2]), str(pk[2]+1), ''.join(neighna), \
                                     str(1000 if sp_options['4NA'][nak][pk][0]>1000 else sp_options['4NA'][nak][pk][0]), \
                                     pk[1], str(pk[2]), str(pk[2]+1), '0,0,0', str(sp_options['4NA'][nak][pk][0]), \
                                     ('%d' % (100*sp_options['4NA'][nak][pk][1]/(sp_options['4NA'][nak][pk][0] if sp_options['4NA'][nak][pk][0]>0 else 1))), \
                                     str(sp_options['4NA'][nak][pk][1]), '\n' ]))
#
# prediction manager of a multiprocess process
#
def mDetect_manager(moptions):
   pmanager = multiprocessing.Manager();
   # get input folder
   while (not moptions['wrkBase']==None) and len(moptions['wrkBase'])>0 and moptions['wrkBase'][-1] in ['/', '\\']:
      moptions['wrkBase'] = moptions['wrkBase'][:-1]

   # need to make prediction of modification
   if moptions['predDet']==1:

      # get well-trained model
      if moptions['modfile'].rfind('/')==-1:
         moptions['modfile'] = [moptions['modfile'], './']
      else:
         moptions['modfile'] = [moptions['modfile'], moptions['modfile'][:moptions['modfile'].rfind('/')+1]]

      start_time = time.time();

      # get fast5 files in a recurisive way
      f5files = glob.glob(os.path.join(moptions['wrkBase'],"*.fast5" ))
      if moptions['recursive']==1:
         f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*.fast5" )))
         f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*/*.fast5" )))
         f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*/*/*.fast5" )))

      print('Total files=%d' % len(f5files))

      # output folder
      if not os.path.isdir(moptions['outFolder']+moptions['FileID']):
         os.system('mkdir -p '+moptions['outFolder']+moptions['FileID'])

      # prepare multiprocessing
      h5files_Q = pmanager.Queue();
      file_map_info_q = pmanager.Queue();
      failed_Q = pmanager.Queue()

      # spliting fast5 files into different lists
      h5_batch = []; h5batchind = 0;
      sub_folder_size = 100; sub_folder_id = 0;
      for f5f in f5files:
         h5_batch.append(f5f);
         if len(h5_batch)==moptions['files_per_thread']:
            # each batch
            h5files_Q.put((h5_batch, sub_folder_id, h5batchind))
            h5_batch = []; h5batchind += 1
            if h5batchind % sub_folder_size ==0:
               sub_folder_id += 1
      if len(h5_batch)>0:
         h5files_Q.put((h5_batch, sub_folder_id, h5batchind))
         h5_batch = []; h5batchind += 1

      # start multiprocessing
      share_var = (moptions, h5files_Q, failed_Q, file_map_info_q)
      handlers = []
      for hid in range(moptions['threads']):
         p = multiprocessing.Process(target=detect_handler, args=share_var);
         p.start();
         handlers.append(p);

      # check potential errors
      failed_files = defaultdict(list);
      while any(p.is_alive() for p in handlers):
         try:
            errk, fns = failed_Q.get(block=False);
            failed_files[errk].extend(fns)

         except:
            time.sleep(1);
            continue;

      # prepare modificatoin summary for reference positions of interest
      moptions['predpath'] = moptions['outFolder'] + '/'+moptions['FileID']
      pred_ind_pref = moptions['outFolder'] + '/'+moptions['FileID']+'/'+pre_base_str
      pred_chr_files = glob.glob(os.path.join(moptions['outFolder']+moptions['FileID'], '*/*.'+pre_base_str+'.*'))
      chr_dict = defaultdict(list);
      for pcf in pred_chr_files:
         chr_dict[ pcf.split('/')[-1].split('.'+pre_base_str)[0] ].append(pcf)
      chr_keys = chr_dict.keys();
      for ck in chr_keys:
         cur_ind_pred_f = pred_ind_pref + '.' + ck;
         cur_list = [ ['#base_folder_fast5', moptions['wrkBase']], ['#base_folder_output', os.path.abspath(moptions['outFolder']+moptions['FileID'])] ]
         for sub_c_f in chr_dict[ck]:
            with open(sub_c_f, 'r') as mr:
                line = mr.readline()
                while line:
                   line = line.strip()
                   if len(line)>0:
                      lsp = line.split();
                      lsp[2] = int(lsp[2])
                      cur_list.append(lsp)
                   line = mr.readline()
         cur_list = sorted(cur_list)
         with open(cur_ind_pred_f, 'w') as indf_writer:
             for mfi in cur_list:
                cur_m_f = []
                for mfidetail in mfi:
                   cur_m_f.append(str(mfidetail))
                cur_m_f.append('\n')
                indf_writer.write(' '.join(cur_m_f))
      # error info
      if len(failed_files)>0:
         print ('Error information for different fast5 files:')
         for errtype, errfiles in failed_files.items():
            print ('\t'+errtype, len(errfiles))

      moptions['outFolder'] = moptions['outFolder']+moptions['FileID']
      end_time = time.time();
      print ("Per-read Prediction consuming time %d" % (end_time-start_time))

   ### for summarizing modificatoin prediction
   start_time = time.time();
   # get all index files of prediction
   all_chr_ind_files = glob.glob(os.path.join(moptions['predpath'], pre_base_str+'.*'))
   print('Find: %s %d %s' % (moptions['predpath'], len(all_chr_ind_files), pre_base_str))
   print (all_chr_ind_files)

   # for each chromosome, a thread will be initialized for multiprocessing summarization of modifications
   chr_strand_Q = pmanager.Queue(); jobnum = 0;
   for cur_cif in all_chr_ind_files:
      chr_strand_Q.put((cur_cif, cur_cif.split(pre_base_str)[-1][1:], '+'))
      chr_strand_Q.put((cur_cif, cur_cif.split(pre_base_str)[-1][1:], '-'))
      jobnum +=2

   # star to summarize modificaiton prediction of reference genomes of interest
   share_var = (moptions, chr_strand_Q)
   handlers = []
   for hid in range(moptions['threads'] if moptions['threads']<jobnum else jobnum):
      p = multiprocessing.Process(target=sum_handler, args=share_var);
      p.start();
      handlers.append(p);
   while any(p.is_alive() for p in handlers):
      try:
         time.sleep(1);
      except:
         time.sleep(1);
         continue;

   end_time = time.time();
   print ("Genomic-position Detection consuming time %d" % (end_time-start_time))

   os.system('touch '+moptions['outFolder']+'.done')

# for independent testing of code
if __name__=='__main__':
#   if len(sys.argv)>4:
      moptions = {}
      moptions['basecall_1d'] = 'Basecall_1D_000'
      moptions['basecall_1d'] = ['Basecall_1D_000']
      moptions['basecall_2strand'] = 'BaseCalled_template'

      moptions['outLevel'] = myCom.OUTPUT_WARNING
      moptions['outLevel'] = myCom.OUTPUT_INFO

      moptions['modfile'] = '../../mod_output/train1/2/mod_train'

      moptions['fnum'] = 3;
      moptions['hidden'] = 100;
      moptions['windowsize'] = 21;

      moptions['threads'] = 8
      moptions['threads'] = 1
      moptions['files_per_thread'] = 500

      mDetect_manager(moptions)
