
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
from . import myDetect

#
# map long reads
# then call another function to get feature for each base of interest
#
def mGetFeature1(moptions, sp_options, f5files):
   # associate signals to events
   f5data = myDetect.get_Event_Signals(moptions, sp_options, f5files)

   # save all sequences information
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG: start_time = time.time();
   temp_fa = tempfile.NamedTemporaryFile(suffix='.fa', mode='w')
   f5keys = sorted(f5data.keys()); #f5keys.sort()
   for f5k in f5keys:
      temp_fa.write(''.join(['>', f5k, '\n', f5data[f5k][0], '\n']))
   temp_fa.flush();
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Write consuming time %d" % (end_time-start_time))

   # alignment using bwa-mem or minimap2
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
   # get sam information
   align_info = temp_sam.readlines()
   align_info = [str(align_info[i], 'utf-8').strip() for i in range(len(align_info))]
   temp_sam.close();

   sp_param = defaultdict();
   sp_param['f5data'] = f5data

   # for alignment
   f5align = defaultdict()
   f5keydict = defaultdict();
   sp_param['ref_info'] = defaultdict()

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:start_time = time.time();
   ilid = 0;
   # for each record in sam, get alignment information
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

   # for unmapped reads
   for f5k in f5keys:
      if f5k not in f5keydict:
         sp_options["Error"]["Not in alignment sam"].append(f5data[f5k][3])

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Get BAM consuming time %d" % (end_time-start_time))

   sp_param['f5status']= ""
   sp_param['line'] = ""
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:start_time = time.time();
   # handle each alignment
   handle_record(moptions, sp_options, sp_param, f5align, f5data)
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Analyze & annotate & save consuming time %d" % (end_time-start_time))

#
# get mapping information
# then call another function to get feature of each base in a long read
#
def handle_record(moptions, sp_options, sp_param, f5align, f5data):
   alignkeys = list(f5align.keys());
   # alignment detail
   numreg = re.compile('\d+')
   mdireg = re.compile('[MIDNSHPX=]{1}')

   feat_file_ind_dict = []
   feat_list = None; feat_file_ind = 0
   start_c_time = time.time();

   for readk in alignkeys:
     if len(feat_file_ind_dict)>0 and feat_list.nbytes > moptions['size_per_batch']:
        # save features when the size is larger than the defined size
        cur_feat_file_base = sp_options['ctfolder'] + '/'+str(feat_file_ind)
        np.savetxt(cur_feat_file_base+'.xy.gz', feat_list, fmt='%.3f')
        with open(cur_feat_file_base+'.xy.ind', 'w') as ind_mw:
            for f_ind in feat_file_ind_dict:
               ind_mw.write('%d %s\n' % (f_ind[1], f_ind[0]))
        print ("\t%s-%d Total consuming time %d" % (sp_options['ctfolder'][sp_options['ctfolder'].rfind('/'):], feat_file_ind, time.time()-start_c_time)); sys.stdout.flush()
        feat_file_ind_dict = []
        feat_list = None;
        feat_file_ind += 1

     # get alignment detail
     mapq, flag, rname, pos, cigar, readseq = f5align[readk]

     if not ( (rname in moptions['fulmodlist'] and len(moptions['fulmodlist'][rname])>0) or \
        ((not moptions['anymodlist']==None) and rname in moptions['anymodlist'] and len(moptions['anymodlist'][rname])>0) or \
        ((not moptions['nomodlist']==None) and rname in moptions['nomodlist'] and len(moptions['nomodlist'][rname])>0) ):
        continue;

     # get reference sequece
     if rname not in sp_param['ref_info']:
        myDetect.getRefSeq(moptions, sp_param, rname)
     refseq = sp_param['ref_info'][rname]

     # mapped starting position and strand
     pos = pos - 1
     forward_reverse = '-' if flag&0x10 else '+'

     numinfo = numreg.findall(cigar);
     mdiinfo = mdireg.findall(cigar)
     numinfo = [int(numinfo[i]) for i in range(len(numinfo))] #map(int, numinfo)

     # remove clip from both tails
     leftclip = 0; rightclip = 0;
     while mdiinfo[0] in ['I', 'D', 'N', 'S', 'H', 'P', 'X']:
         if mdiinfo[0] in ['I', 'S', 'X']:
            leftclip += numinfo[0];  readseq = readseq[numinfo[0]:]
         if mdiinfo[0] in ['H']: leftclip += numinfo[0]
         if mdiinfo[0] in ['D', 'N', 'X']:
            pos += numinfo[0]
         numinfo = numinfo[1:];  mdiinfo = mdiinfo[1:]
     while mdiinfo[-1] in ['I', 'D', 'N', 'S', 'H', 'P', 'X']:
         if mdiinfo[-1] in ['I', 'S', 'X']:
            rightclip += numinfo[-1]; readseq = readseq[:-numinfo[-1]]
         if mdiinfo[-1] in ['H']: rightclip += numinfo[-1]
         numinfo = numinfo[:-1]; mdiinfo = mdiinfo[:-1]
     if forward_reverse=='+':
        if rightclip>0: m_event = f5data[readk][1][leftclip:-rightclip]
        else: m_event = f5data[readk][1][leftclip:]
     else:
        if leftclip>0: m_event = f5data[readk][1][rightclip:-leftclip]
        else: m_event = f5data[readk][1][rightclip:]

     # is in region of interest if provided
     isinreg = False;
     if (moptions['region'][0] in ['', None, rname]) and \
        (moptions['region'][1] in ['', None] or pos>moptions['region'][1]) and \
        (moptions['region'][2] in ['', None] or pos+len(m_event)<moptions['region'][2]):
        isinreg = True;
     if not isinreg:
        continue;

     # associate mapped reference positions with read positions
     lastmatch = None; firstmatch = None;
     first_match_pos = None; last_match_pos = None
     last_al_match = None;  first_al_match = None
     lasmtind = 0;
     base_map_info = []; #indel_groups = defaultdict()
     nummismatch = 0; numinsert = 0; numdel = 0;
     read_ind = 0;
     for n1ind in range(len(numinfo)):
        mdi = mdiinfo[n1ind];
        # for each mapped types
        for n1i in range(numinfo[n1ind]):
           if mdi=='M':
              base_map_info.append((refseq[pos], readseq[read_ind], pos, read_ind))
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
              base_map_info.append(('-', readseq[read_ind], pos, read_ind))
              read_ind += 1;
              numinsert += 1
           elif mdi == 'D':
              base_map_info.append((refseq[pos], '-', pos, read_ind))
              pos += 1;
              numdel += 1
           elif mdi == 'N':
              base_map_info.append((refseq[pos], '-', pos, read_ind))
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
             base_map_info.append((refseq[pos], readseq[read_ind], pos, read_ind))
             if first_match_pos==None: first_match_pos  = pos
             if last_match_pos==None or last_match_pos<pos: last_match_pos = pos
             pos += 1; read_ind += 1;
             if firstmatch==None: firstmatch = read_ind - 1
             if lastmatch==None or lastmatch<read_ind-1: lastmatch = read_ind - 1; lasmtind=n1ind
             if last_al_match==None or last_al_match<len(base_map_info): last_al_match=len(base_map_info)-1
             if first_al_match==None: first_al_match=len(base_map_info)-1
           elif mdi == 'X':
             base_map_info.append((refseq[pos], readseq[read_ind], pos, read_ind))
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

     # remove unmatch in both tails
     if not firstmatch==None: leftclip += firstmatch
     if (not lastmatch==None) and len(m_event)-lastmatch>1: rightclip += len(m_event)-lastmatch-1
     # remove events whose bases are not mapped.
     if forward_reverse=='+':
        if len(m_event)-lastmatch>1:
           m_event = m_event[firstmatch:(lastmatch+1-len(m_event))]
        elif firstmatch>0: m_event = m_event[firstmatch:]
     else:
        if firstmatch>0: m_event = m_event[(len(m_event)-1-lastmatch):-firstmatch]
        elif len(m_event)-lastmatch>1: m_event = m_event[(len(m_event)-1-lastmatch):]
     # print detail if unexpected error occurs
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

     # post-process mapping information
     base_map_info = np.array(base_map_info, dtype=[('refbase', 'U1'), ('readbase', 'U1'), ('refbasei', np.uint64), ('readbasei', np.uint64)])
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
     #
     # handle map like
     # CCG    or CGG
     # C-G       C-G
     #
     if 'motif' in moptions and moptions['motif'][0]=='CG':
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
     # too short reads
     if len(m_event)<500:
         sp_options["Error"]["Less(<500) events"].append(f5data[readk][3])
         continue;

     # get feautre
     mfeatures,isdif = get_Feature(moptions, sp_options, sp_param, f5align, f5data, readk, leftclip, rightclip, base_map_info, forward_reverse, rname, first_match_pos, numinsert, numdel)
     if isdif and moptions['outLevel']<=myCom.OUTPUT_WARNING:
        print("Dif is true")
        print([lastmatch, firstmatch, first_match_pos, last_match_pos, first_al_match, last_al_match, lasmtind, len(base_map_info), nummismatch, numinsert, numdel, len(base_map_info)-nummismatch-numinsert-numdel])

     # merge to previously handled features of other fast5 files
     if len(mfeatures)>0:
        if len(feat_file_ind_dict)==0:
           feat_file_ind_dict.append((f5data[readk][3], 0));
           feat_list = mfeatures
        else:
           feat_file_ind_dict.append((f5data[readk][3], len(feat_list)))
           feat_list = np.concatenate((feat_list, mfeatures), axis=0)

   # store the last feature data.
   if len(feat_file_ind_dict)>0:
      cur_feat_file_base = sp_options['ctfolder'] + '/'+str(feat_file_ind)
      np.savetxt(cur_feat_file_base+'.xy.gz', feat_list, fmt='%.3f')
      with open(cur_feat_file_base+'.xy.ind', 'w') as ind_mw:
          for f_ind in feat_file_ind_dict:
             ind_mw.write('%d %s\n' % (f_ind[1], f_ind[0]))
      print ("\t%s-%d Total consuming time %d" % (sp_options['ctfolder'][sp_options['ctfolder'].rfind('/'):], feat_file_ind, time.time()-start_c_time)); sys.stdout.flush()
      feat_file_ind_dict = []
      feat_list = None;
      feat_file_ind += 1

#
# get feature for each base of interest in long reads according to raw signals and mapping information
#
def get_Feature(moptions, sp_options, sp_param, f5align, f5data, readk, start_clip, end_clip, base_map_info, forward_reverse, rname, mapped_start_pos, num_insertions, num_deletions):
   # event information
   modevents = sp_param['f5data'][readk][1]
   # class number, bin num and bin length
   clnum = 2; binnum = 50; binlen = 0.2;
   if forward_reverse=='+':
      align_ref_pos = mapped_start_pos
   else:
      align_ref_pos = mapped_start_pos + len(base_map_info) - num_insertions - 1

   # initialize feature matrix for all events.
   if moptions['fnum']==57:
      #mfeatures = np.zeros((len(modevents)-end_clip+100-(start_clip-100), (binnum+3+3+4)));
      mfeatures = np.zeros((len(modevents)-end_clip+100-(start_clip-100), (binnum+3+3+4)));
   else: mfeatures = np.zeros((len(modevents)-end_clip+100-(start_clip-100), (3+3+4)));

   # filter poor alignment
   checkneighbornums = [3,6]
   checkratios = {3:[6,5,4,2], 6:[11,10,9,3]}
   checkratios = {3:[6,5,4,2], 6:[12,10,9,3]}
   cgpos = [[], []]
   affectneighbor = 1; # 2;
   for aligni in range(len(base_map_info)):
      # for methylated positions and not-used adjacent positions
      if 'motif' in moptions and base_map_info['readbase'][aligni]==moptions['motif'][0][moptions['motif'][1]]:
         m_a_st = aligni-moptions['motif'][1]; m_a_end = aligni+len(moptions['motif'][0])-moptions['motif'][1]
         if m_a_st>-1 and m_a_end<=len(base_map_info) and ''.join(base_map_info['readbase'][m_a_st:m_a_end])==moptions['motif'][0] and (not ''.join(base_map_info['refbase'][m_a_st:m_a_end])==moptions['motif'][0]):
            cgpos[1].extend([(forward_reverse, base_map_info['refbasei'][addi]) for addi in range(aligni-affectneighbor if aligni-affectneighbor>-1 else 0, aligni+affectneighbor+1 if aligni+affectneighbor+1<len(base_map_info) else len(base_map_info))])
      if (not base_map_info['refbase'][aligni]=='-') and \
         (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['fulmodlist'][rname]:
         if not base_map_info['readbase'][aligni]=='-':
            nextnogap = aligni + 1;
            while nextnogap<len(base_map_info):
               if not base_map_info['refbase'][nextnogap]=='-': break;
               nextnogap += 1
            iscg = False;
            # find gaps
            for checkneighbornum in checkneighbornums:
               if not nextnogap<len(base_map_info): continue;
               matchnum = 0; gapnum = 0;
               # get gaps for two window sizes
               for checki in range(aligni-checkneighbornum, aligni+checkneighbornum+1):
                  if checki>-1 and checki<len(base_map_info):
                     if base_map_info['refbase'][checki]==base_map_info['readbase'][checki]: matchnum += 1
                     if base_map_info['refbase'][checki]=='-' or base_map_info['readbase'][checki]=='-': gapnum += 1
               if gapnum<=checkratios[checkneighbornum][3]:
                  for addi in range(aligni-affectneighbor if aligni-affectneighbor>-1 else 0, nextnogap+affectneighbor if nextnogap+affectneighbor<len(base_map_info) else len(base_map_info)):
                     if addi==aligni: # for methylated positions
                        cgpos[0].append((forward_reverse, base_map_info['refbasei'][addi]))
                     else: # for non-used positions
                        cgpos[1].append((forward_reverse, base_map_info['refbasei'][addi]))
                  iscg = True; break;
            if iscg: continue;
         # add more not-used positions if more gaps exist
         if not base_map_info['readbase'][aligni]=='-':
            nextnogap = aligni
            for _ in range(affectneighbor):
               nextnogap += 1;
               while nextnogap<len(base_map_info['refbase']):
                 if not base_map_info['refbase'][nextnogap]=='-': break;
                 nextnogap += 1
            prenogap = aligni
            for _ in range(affectneighbor):
               prenogap -= 1;
               while prenogap>-1:
                  if not base_map_info['refbase'][prenogap]=='-': break;
                  prenogap -= 1

            read0 = aligni; read1 = aligni
            for _ in range(affectneighbor):
               read0 -= 1
               while read0>-1:
                  if base_map_info['readbase'][read0]=='-': read0 -= 1
                  else: break;
               read1 += 1
               while read1<len(base_map_info['readbase']):
                  if base_map_info['readbase'][read1]=='-': read1 += 1
                  else: break;

            if read0<prenogap:
               if read0>-1: prenogap = read0
               else: prenogap = 0
            if read1>nextnogap:
               if read1<len(base_map_info['readbase']): nextnogap = read1
               else: nextnogap = len(base_map_info['readbase'])-1
            if prenogap<0: prenogap = 0
            if not nextnogap<len(base_map_info['readbase']): nextnogap=len(base_map_info['readbase'])-1
            if not prenogap<len(base_map_info['readbase']): prenogap=len(base_map_info['readbase'])-1
            for excldi in range(prenogap, nextnogap+1):
               cgpos[1].append((forward_reverse, base_map_info['refbasei'][excldi]))

   print ('%s%s %d, %d >> %d %d, %d-%d=%d' % (forward_reverse, f5data[readk][3], len(cgpos[0]), len(cgpos[1]), len(modevents)-end_clip-start_clip, start_clip, len(modevents), end_clip, len(modevents)-end_clip))

   aligni = 0; isdif = False;
   for ie in range(start_clip-100, len(modevents)-end_clip+100):
      cur_row_num = ie - (start_clip-100); cur_base = ''
      # for aligned bases
      if ie>=start_clip and ie<len(modevents)-end_clip:
         if align_ref_pos<mapped_start_pos:
            print ('ERRRR align_ref_pos(%d)<mapped_start_pos(%d)' % (align_ref_pos, mapped_start_pos))
         while base_map_info['readbase'][aligni]=='-':
            if not align_ref_pos==base_map_info['refbasei'][aligni]:
               print ('ERRRR align_ref_pos(%d) not equal to %d' % (align_ref_pos, base_map_info['refbasei'][aligni] ))
            if not base_map_info['refbase'][aligni]=='-':
               if forward_reverse=='+': align_ref_pos += 1
               else: align_ref_pos -= 1
            aligni += 1
         if not base_map_info['readbase'][aligni] == modevents['model_state'][ie][2]:
            print ('Error Does not match', base_map_info['readbase'][aligni], modevents['model_state'][ie][2], aligni, ie)
            isdif = True;
         # the first column is the aligned reference position
         mfeatures[cur_row_num][0] = align_ref_pos
         cur_base = base_map_info['refbase'][aligni]
         # the second/third column is for negative/positive labels of methylation
         if moptions['posneg'] == 0: # for a data without any modification
            if ( (not moptions['anymodlist']==None) and rname in moptions['nomodlist'] and (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['nomodlist'][rname] ):
                mfeatures[cur_row_num][1] = 1; mfeatures[cur_row_num][2] = 0
            elif (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['fulmodlist'][rname]:
                mfeatures[cur_row_num][1] = 1; mfeatures[cur_row_num][2] = 0
            elif ((not moptions['anymodlist']==None) and rname in moptions['anymodlist'] and (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['anymodlist'][rname] ):
                mfeatures[cur_row_num][1] = 1; mfeatures[cur_row_num][2] = 0
         else: # for a data with both modified and un-modified positions
            if (forward_reverse, base_map_info['refbasei'][aligni]) in cgpos[0] and (not base_map_info['refbase'][aligni]=='-'):
               mfeatures[cur_row_num][1] = 0; mfeatures[cur_row_num][2] = 1
            else:
               if (forward_reverse, base_map_info['refbasei'][aligni]) not in cgpos[1]:
                  if moptions['anymodlist']==None:
                      if moptions['nomodlist']==None or ( rname in moptions['nomodlist'] and (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['nomodlist'][rname] ):
                         mfeatures[cur_row_num][1] = 1; mfeatures[cur_row_num][2] = 0
                  elif rname in moptions['anymodlist'] and (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['anymodlist'][rname]:
                      pass
                  else:
                      if moptions['nomodlist']==None or ( rname in moptions['nomodlist'] and (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['nomodlist'][rname] ):
                         mfeatures[cur_row_num][1] = 1; mfeatures[cur_row_num][2] = 0
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
         # for reference base type feature
         if cur_base in myCom.g_ACGT:
            mfeatures[cur_row_num][moptions['fnum']-3+3-4+myCom.g_ACGT.index(cur_base)] = 1
         cur_index_add = moptions['fnum'] - 3 + 3
         # for signal mean std and length.
         mfeatures[cur_row_num][cur_index_add + 0] = modevents["mean"][ie]
         mfeatures[cur_row_num][cur_index_add + 1] = modevents["stdv"][ie]
         mfeatures[cur_row_num][cur_index_add + 2] = modevents["length"][ie]

   # truncated too much not-used positions
   nbkeys = defaultdict();
   for mfi in range(len(mfeatures)):
      if mfeatures[mfi][1] + mfeatures[mfi][2] > 0.9:
         for ini in range(mfi-25, mfi+26):
            if ini<0 or ini>len(mfeatures)-1:
               print("Warning wrong del mfeatures id %d for %s" % (ini, f5data[readk][3]))
            else:
               nbkeys[ini] = True;
   keepInd = sorted(list(nbkeys.keys()));
   if len(keepInd)>0:
      if not len(keepInd)>len(mfeatures)*0.9:
         mfeatures = mfeatures[np.array(keepInd)]
   else:
      mfeatures = []

   return (mfeatures, isdif)


#
# get the complementary bases
#
def get_complement(na):
   if na in myCom.acgt: return myCom.na_bp[na]
   else: return na;

#
# get required information for reach mapping records.
#
def handle_line(moptions, sp_param, f5align):
   lsp = sp_param['line'].split('\t')
   qname, flag, rname, pos, mapq, cigar, _, _, _, seq, _ = lsp[:11]
   # checked query name
   if qname=='*': sp_param['f5status'] = "qname is *"
   # check mapping quality
   elif int(mapq)==255: sp_param['f5status'] = "mapq is 255"
   # check mapped positions
   elif int(pos)==0: sp_param['f5status'] = "pos is 0"
   # check mapped string
   elif cigar=='*': sp_param['f5status'] = "cigar is *"
   # check reference name
   elif rname=='*': sp_param['f5status'] = "rname is *"
   if not sp_param['f5status']=="": return qname

   if (qname not in f5align) or f5align[qname][0]<int(mapq):
      f5align[qname] = (int(mapq), int(flag), rname, int(pos), cigar, seq)

   return qname

#
# feature handler/workder for multiprocessing
#
def getFeature_handler(moptions, h5files_Q, failed_Q, version_Q):
   while not h5files_Q.empty():
      try:
         # get a list of files
         f5files, ctfolderid = h5files_Q.get(block=False)
      except:
         break;

      sp_options = defaultdict();
      sp_options['ctfolder'] = moptions['outFolder']+str(ctfolderid)
      if not os.path.isdir(sp_options['ctfolder']):
         os.system('mkdir '+sp_options['ctfolder'])
      # get features
      mGetFeature1(moptions, sp_options, f5files)
      # output errors
      for errtype, errfiles in sp_options["Error"].items():
         failed_Q.put((errtype, errfiles));
      # double check albacore version
      for vk in sp_options["get_albacore_version"]:
         version_Q.put((vk, sp_options["get_albacore_version"][vk]))

#
# read sequence information from a reference genome
#
def readFA(mfa, t_chr=None):
   fadict = defaultdict();
   with open(mfa, 'r') as mr:
      cur_chr = None;
      line = mr.readline();
      while line:
         # remove empty spaces
         line = line.strip();
         if len(line)>0:
            if line[0]=='>': # for each chromosome line
               if (not cur_chr==None) and (t_chr in [None, cur_chr]):
                  fadict[cur_chr] = ''.join(fadict[cur_chr])
               cur_chr = line[1:].split()[0]
               if t_chr in [None, cur_chr]:
                  fadict[cur_chr] = []
            else: # for sub-sequence line in a reference file
               if t_chr in [None, cur_chr]:
                  fadict[cur_chr].append(line.upper())
         line = mr.readline();
      # for the last chromosome in the file
      if (not cur_chr==None) and (t_chr in [None, cur_chr]):
         fadict[cur_chr] = ''.join(fadict[cur_chr])
   return fadict

#
# get reference positions for motif-based modifications
#
def readMotifMod(fadict, mpat='Cg', mposinpat=0, t_chr=None, t_start=None, t_end=None):
   pos_dict = defaultdict(int)

   # get motif and complementary motif
   pat3 = copy.deepcopy(mpat.upper())
   comp_pat3 = ''.join([get_complement(curna) for curna in pat3][::-1])
   comp_mposinpat = len(comp_pat3)-1-mposinpat

   fakeys = fadict.keys();
   cpgdict = defaultdict(int);
   all_a = defaultdict()
   for fak in fakeys:
       cpgnum = [0, 0]
       # motif-based reference positions
       cpgdict[fak] = defaultdict()
       # position of bases of interest
       all_a[fak] = defaultdict()
       for i in range(len(fadict[fak])):
          if (t_start==None or i>=t_start) and (t_end==None or i<=t_end):
             if fadict[fak][i]==mpat[mposinpat]: # for forward strand
                all_a[fak][('+', i)] = True;
             elif get_complement(fadict[fak][i])==mpat[mposinpat]: # for reverse strand
                all_a[fak][('-', i)] = True;

             # check motif in forward strand
             if i-mposinpat>=0 and i+len(comp_pat3)-1-mposinpat<len(fadict[fak]) and ''.join(fadict[fak][i-mposinpat:(i+len(comp_pat3)-1-mposinpat+1)])==pat3:
                cpgdict[fak][('+', i)] = [1, fadict[fak][i]]; cpgnum[0] += 1
             elif i-comp_mposinpat>=0 and i+len(comp_pat3)-1-comp_mposinpat<len(fadict[fak]) and ''.join(fadict[fak][i-comp_mposinpat:(i+len(comp_pat3)-1-comp_mposinpat+1)])==comp_pat3: # check motif in reverse strand
                cpgdict[fak][('-', i)] = [1, fadict[fak][i]]; cpgnum[1] += 1
             else:
                pass
       print('%s%d site: %d(+) %d(-) for %s' % (pat3, mposinpat, cpgnum[0], cpgnum[1], fak))
   return (cpgdict, all_a)


#
# a multiprocessing manager to get features from all long reads.
#
def getFeature_manager(moptions):
   start_time = time.time();
   # multipprocessing manager
   pmanager = multiprocessing.Manager();

   # prepare output folder
   if os.path.isdir(moptions['outFolder']):
      os.system('rm -dr '+moptions['outFolder'])
   if not os.path.isdir(moptions['outFolder']):
      os.system('mkdir '+moptions['outFolder'])

   moptions['size_per_batch'] = moptions['size_per_batch'] * (10**7)

   # read reference information
   fadict = readFA(moptions['Ref'],moptions['region'][0])
   if moptions['motifORPos']==1: # get motif-based positions for modifications
      moptions['fulmodlist'], moptions['nomodlist'] = readMotifMod(fadict, moptions['motif'][0], moptions['motif'][1], moptions['region'][0], moptions['region'][1], moptions['region'][2])
      moptions['anymodlist'] = None
      moptions['nomodlist'] = None; # add for simple process
   elif moptions['motifORPos']==2: # modification position is specified by the files
      fuldfiles = glob.glob(moptions["fulmod"]);
      moptions['fulmodlist'] = defaultdict(lambda: defaultdict());
      if not moptions["anymod"]==None: # partially modified positions
         anydfiles = glob.glob(moptions["anymod"])
         moptions['anymodlist'] = defaultdict(lambda: defaultdict());
      else:
         moptions['anymodlist'] = None
      if not moptions["nomod"]==None: # completely un-modified positions
         nodfiles = glob.glob(moptions["nomod"])
         moptions['nomodlist'] = defaultdict(lambda: defaultdict());
      else:
         moptions['nomodlist'] = None
      mthreadin = [moptions['fulmodlist'], moptions['anymodlist'], moptions['nomodlist']]
      mthfiles = [fuldfiles, anydfiles, nodfiles]
      # read completely modified positions, partially modified positions, completely un-modified positions from files
      for mthi in range(len(mthreadin)):
         curmeth = mthreadin[mthi]; curfilelist = mthfiles[mthi]
         if curmeth==None or curfilelist==None: continue;
         for curmthf in curfilelist:
             with open(curmthf, 'r') as mreader:
                line = mreader.readline();
                while line:
                   if len(line)>0:
                      tchr, tstrand, tpos = line.split()[:3]
                      curmeth[tchr][(tstrand, int(tpos))] = [1-mthi, fadict[tchr][int(tpos)]];
                   line = mreader.readline();
   for tchr in moptions['fulmodlist'] if moptions['anymodlist']==None else moptions['anymodlist']:
      if len(moptions['fulmodlist'][tchr])>0 or ((not moptions['anymodlist']==None) and len(moptions['anymodlist'][tchr])>0):
          print ('%s fulmod=%d anymod=%d nomod=%d' % (tchr, len(moptions['fulmodlist'][tchr]), len(moptions['anymodlist'][tchr]) if (not moptions['anymodlist']==None) else -1, len(moptions['nomodlist'][tchr]) if (not moptions['nomodlist']==None) else -1))

   if True: #False:
      # get all input fast5 files
      f5files = glob.glob(os.path.join(moptions['wrkBase'],"*.fast5" ))
      if moptions['recursive']==1:
         f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*.fast5" )))
         f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*/*.fast5" )))
         f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*/*/*.fast5" )))


   print('Total files=%d' % len(f5files))
   h5files_Q = pmanager.Queue();
   failed_Q = pmanager.Queue()
   version_Q = pmanager.Queue()

   # split input fast5 files into different batch
   h5_batch = []; h5batchind = 0;
   for f5f in f5files:
      h5_batch.append(f5f);
      if len(h5_batch)==moptions['files_per_thread']:
         h5files_Q.put((h5_batch, h5batchind))
         h5batchind += 1
         h5_batch = []; #break; ### feature500
   if len(h5_batch)>0:
      h5files_Q.put((h5_batch, h5batchind))

   # each thread handle a batch a time and repeat for all batches.
   share_var = (moptions, h5files_Q, failed_Q, version_Q)
   handlers = []
   for hid in range(moptions['threads']):
      p = multiprocessing.Process(target=getFeature_handler, args=share_var);
      p.start();
      handlers.append(p);

   # get failed files.
   failed_files = defaultdict(list);
   version_default = defaultdict(lambda: defaultdict(int));
   while any(p.is_alive() for p in handlers):
      try:
         errk, fns = failed_Q.get(block=False);
         failed_files[errk].extend(fns)
         curv, curv_num = version_Q.get(block=False);
         version_default[curv] += curv_num
      except:
         time.sleep(1);
         continue;

   # output failure information
   if len(failed_files)>0:
      print ('Error information for different fast5 files:')
      for errtype, errfiles in failed_files.items():
         print ('\t%s %d' % (errtype, len(errfiles)))
   print("abversion info {}".format(str(version_default)))
   sys.stdout.flush()
   end_time = time.time();
   print ("Total consuming time %d" % (end_time-start_time))



# for indepdent testing of code
if __name__=='__main__':
#   if len(sys.argv)>4:
      moptions = {}
      moptions['basecall_1d'] = 'Basecall_1D_000'
      moptions['basecall_1d'] = ['Basecall_1D_000']
      moptions['basecall_2strand'] = 'BaseCalled_template'

      moptions['outLevel'] = myCom.OUTPUT_WARNING
      moptions['outLevel'] = myCom.OUTPUT_INFO

      moptions['modfile'] = '../../mod_output/train1/2/mod_train'

      moptions['fnum'] = 53;
      moptions['hidden'] = 100;
      moptions['windowsize'] = 21;

      moptions['threads'] = 8
      moptions['threads'] = 1
      moptions['files_per_thread'] = 500

      mDetect_manager(moptions)
