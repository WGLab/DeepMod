
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
#import myCom


def mGetFeature1(moptions, sp_options, f5files):
   f5data = myDetect.get_Event_Signals(moptions, sp_options, f5files)

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG: start_time = time.time();
   temp_fa = tempfile.NamedTemporaryFile(suffix='.fa', mode='w')
   f5keys = sorted(f5data.keys()); #f5keys.sort()
   for f5k in f5keys:
      temp_fa.write(''.join(['>', f5k, '\n', f5data[f5k][0], '\n']))
   temp_fa.flush();   
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Write consuming time %d" % (end_time-start_time))
   
   temp_sam = tempfile.NamedTemporaryFile()
   #cmd_opt = ['mem', '-x', 'ont2d', '-v', '1', '-t', '1', moptions['Ref'], temp_fa.name]
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
   while ilid < len(align_info):
      if len(align_info[ilid])==0 or align_info[ilid][0]=='@':
         ilid += 1
         continue;

      sp_param['f5status'] = "";
      sp_param['line'] = align_info[ilid]
      qname = myDetect.handle_line(moptions, sp_param, f5align)
      if sp_param['f5status'] == "":
         f5keydict[qname] = True;
      ilid += 1

   for f5k in f5keys:
      if f5k not in f5keydict:
         sp_options["Error"]["Not in alignment sam"].append(f5data[f5k][3])

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Get BAM consuming time %d" % (end_time-start_time))

   sp_param['f5status']= ""
   sp_param['line'] = ""
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:start_time = time.time();
   handle_record(moptions, sp_options, sp_param, f5align, f5data)
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Analyze & annotate & save consuming time %d" % (end_time-start_time))

def handle_record(moptions, sp_options, sp_param, f5align, f5data):
   alignkeys = f5align.keys();
   numreg = re.compile('\d+')
   mdireg = re.compile('[MIDNSHPX=]{1}')

   feat_file_ind_dict = []
   feat_list = None; feat_file_ind = 0
   start_c_time = time.time();

   for readk in alignkeys:
     if len(feat_file_ind_dict)>0 and feat_list.nbytes > moptions['size_per_batch']:
        cur_feat_file_base = sp_options['ctfolder'] + '/'+str(feat_file_ind)
        np.savetxt(cur_feat_file_base+'.xy.gz', feat_list, fmt='%.3f')
        with open(cur_feat_file_base+'.xy.ind', 'w') as ind_mw:
            for f_ind in feat_file_ind_dict:
               ind_mw.write('%d %s\n' % (f_ind[1], f_ind[0]))
        print ("\t%s-%d Total consuming time %d" % (sp_options['ctfolder'][sp_options['ctfolder'].rfind('/'):], feat_file_ind, time.time()-start_c_time)); sys.stdout.flush()
        feat_file_ind_dict = []
        feat_list = None;
        feat_file_ind += 1

     #print (f5data[readk][3]);
     _, flag, rname, pos, cigar, readseq = f5align[readk]
     if (rname not in moptions['anymodlist']) or \
        len(moptions['fulmodlist'][rname])==0 or \
        len(moptions['anymodlist'][rname])==0:
        continue;

     if rname not in sp_param['ref_info']:
        myDetect.getRefSeq(moptions, sp_param, rname)
     refseq = sp_param['ref_info'][rname]

     pos = pos - 1
     numinfo = numreg.findall(cigar);
     mdiinfo = mdireg.findall(cigar)
     
     forward_reverse = '-' if flag&0x10 else '+'
     numinfo = [int(numinfo[i]) for i in range(len(numinfo))] #map(int, numinfo)

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

     lastmatch = None; firstmatch = None; 
     first_match_pos = None; last_match_pos = None
     last_al_match = None;  first_al_match = None
     lasmtind = 0;
     base_map_info = []; #indel_groups = defaultdict()
     nummismatch = 0; numinsert = 0; numdel = 0;
     read_ind = 0;
     for n1ind in range(len(numinfo)):
        mdi = mdiinfo[n1ind];
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

     if not firstmatch==None: leftclip += firstmatch
     if (not lastmatch==None) and len(m_event)-lastmatch>1: rightclip += len(m_event)-lastmatch-1
      
     if forward_reverse=='+':
        if len(m_event)-lastmatch>1:
           m_event = m_event[firstmatch:(lastmatch+1-len(m_event))]
        elif firstmatch>0: m_event = m_event[firstmatch:]
     else:
        if firstmatch>0: m_event = m_event[(len(m_event)-1-lastmatch):-firstmatch]
        elif len(m_event)-lastmatch>1: m_event = m_event[(len(m_event)-1-lastmatch):]
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

     base_map_info = np.array(base_map_info, dtype=[('refbase', 'U1'), ('readbase', 'U1'), ('refbasei', np.uint64), ('readbasei', np.uint64)])
     if forward_reverse=='-':
        base_map_info = np.flipud(base_map_info)
        for bmii in range(len(base_map_info)):
            base_map_info['refbase'][bmii]  = get_complement(base_map_info['refbase'][bmii])
            base_map_info['readbase'][bmii] = get_complement(base_map_info['readbase'][bmii])
        leftclip, rightclip = rightclip, leftclip
     mfeatures,isdif = get_Feature(moptions, sp_options, sp_param, f5align, f5data, readk, leftclip, rightclip, base_map_info, forward_reverse, rname, first_match_pos, numinsert, numdel)
     #mfeatures = get_Feature(moptions, sp_options, sp_param, f5align, m_event, readk, base_map_info)
     if isdif and moptions['outLevel']<=myCom.OUTPUT_WARNING:
        print("Dif is true")
        print([lastmatch, firstmatch, first_match_pos, last_match_pos, first_al_match, last_al_match, lasmtind, len(base_map_info), nummismatch, numinsert, numdel, len(base_map_info)-nummismatch-numinsert-numdel])

     
     if len(mfeatures)>0:
        if len(feat_file_ind_dict)==0:
           feat_file_ind_dict.append((sp_param['f5f'], 0));
           feat_list = cur_feat_list
        else:
           feat_file_ind_dict.append((sp_param['f5f'], len(feat_list)))
           feat_list = np.concatenate((feat_list, cur_feat_list), axis=0)
       
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

 
def get_Feature(moptions, sp_options, sp_param, f5align, f5data, readk, start_clip, end_clip, base_map_info, forward_reverse, rname, mapped_start_pos, num_insertions, num_deletions):
   modevents = sp_param['f5data'][readk][1]
   clnum = 3; binnum = 50; binlen = 0.2;
   #align_ref_pos = mapped_start_pos
   if forward_reverse=='+':
      align_ref_pos = mapped_start_pos
   else:
      align_ref_pos = mapped_start_pos + len(base_map_info) - num_insertions - 1

  
   if moptions['fnum']==53: 
      mfeatures = np.zeros((len(modevents)-end_clip+100-(start_clip-100), (binnum+3+3)));
   else: mfeatures = np.zeros((len(modevents)-end_clip+100-(start_clip-100), (3+3)));
   
   read_align_list = base_map_info['readbase']
   ref_align_list = base_map_info['refbase']
   for ali in range(len(read_align_list)):
      if ref_align_list[ali]=='C' and read_align_list[ali]=='C':
         if ali+1<len(read_align_list) and read_align_list[ali+1]=='-' and ref_align_list[ali+1]=='G':
            addali = 2;
            while ali + addali < len(read_align_list):
               if read_align_list[ali+addali]=='-' and ref_align_list[ali+addali]=='G': addali += 1;
               else: break;
            if ali + addali < len(read_align_list) and read_align_list[ali+addali]=='G' and ref_align_list[ali+addali]=='G':
               read_align_list[ali+1], read_align_list[ali+addali] = read_align_list[ali+addali], read_align_list[ali+1]
      if ref_align_list[ali]=='G' and read_align_list[ali]=='G':
         if ali-1>-1 and read_align_list[ali-1]=='-' and ref_align_list[ali-1]=='C':
            addali = 2;
            while ali - addali >-1:
               if read_align_list[ali-addali]=='-' and ref_align_list[ali-addali]=='C': addali += 1;
               else: break;
            if ali - addali>-1 and read_align_list[ali-addali]=='C' and ref_align_list[ali-addali]=='C':
               read_align_list[ali-1], read_align_list[ali-addali] = read_align_list[ali-addali], read_align_list[ali-1]
   read_align = ''.join(read_align_list)
   ref_align = ''.join(ref_align_list)
   cgpos = [[], []]
   affectneighbor = 2;
   for aligni in range(len(ref_align)):
      if (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['anymodlist'][rname]:
         if (not ref_align_list[aligni]=='-') and \
            (forward_reverse, base_map_info['refbasei'][aligni]) in moptions['fulmodlist'][rname]:
            if not read_align_list[ali]=='-':
               cgpos[0].append((forward_reverse, base_map_info['refbasei'][aligni])) 
         nextnogap = aligni + 1;
         for _ in range(affectneighbor):
            while nextnogap<len(ref_align):
               if not ref_align[nextnogap]=='-': break;
               nextnogap += 1
         prenogap = aligni - 1;
         for _ in range(affectneighbor):
            while prenogap>-1:
               if not ref_align[prenogap]=='-': break;
               prenogap -= 1

         read0 = aligni-1; read1 = aligni+1
         for _ in range(affectneighbor):
            while read0>-1:
               if read_align[read0]=='-': read0 -= 1
               else: break;
            while read1<len(read_align):
               if read_align[read1]=='-': read1 += 1
               else: break;
         
         if read0<prenogap:
            if read0>-1: prenogap = read0
            else: prenogap = 0
         if read1>nextnogap:
            if read1<len(read_align): nextnogap = read1
            else: nextnogap = len(read_align)-1
         if prenogap<0: prenogap = 0
         for excldi in range(prenogap, nextnogap+1):
            cgpos[1].append((forward_reverse, base_map_info['refbasei'][excldi]))

   if fortest:
            outprint = [[], [], []]; mxlsize = [0, 0, 0]
            for aligni in range(len(ref_align)):
                outprint[0].append(ref_align[aligni])
                outprint[1].append(read_align[aligni])
                if (forward_reverse, base_map_info['refbasei'][aligni]) in cgpos[0]:
                   outprint[2].append('+')
                   mxlsize[1] += 1
                elif (forward_reverse, base_map_info['refbasei'][aligni]) in cgpos[1]:
                   outprint[2].append('X')
                   mxlsize[0] += 1
                else:
                   outprint[2].append('0')
                   mxlsize[2] += 1
            perlnum = 100;
            printline = len(outprint[0])/perlnum;
            for pli in range(1, printline+2):
                if pli == printline+1:
                   print (''.join(outprint[0][(pli-1)*perlnum:]))
                   print (''.join(outprint[1][(pli-1)*perlnum:]))
                   print (''.join(outprint[2][(pli-1)*perlnum:]))
                else:
                   print (''.join(outprint[0][(pli-1)*perlnum:pli*perlnum]))
                   print (''.join(outprint[1][(pli-1)*perlnum:pli*perlnum]))
                   print (''.join(outprint[2][(pli-1)*perlnum:pli*perlnum]))
            print (len(ref_align), mxlsize, len(read_align)) #, featfile)


   print ('%s %d, %d >> %d %d, %d-%d=%d' % (f5data[readk][3], len(cgpos[0]), len(cgpos[1]), len(modevents)-end_clip-start_clip, start_clip, len(modevents), end_clip, len(modevents)-end_clip))

   if not len(modevents)==(start_clip+end_clip+len(read_align)-num_deletions):
      print ('ERRRR %d==%d(%d+%d+%d-%d)' % (len(modevents), start_clip+end_clip+len(read_align)-num_deletions, start_clip, end_clip, len(read_align), num_deletions))

   aligni = 0; isdif = False;
   for ie in range(start_clip-100, len(modevents)-end_clip+100):
      cur_row_num = ie - (start_clip-100)
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
         mfeatures[cur_row_num][0] = align_ref_pos
         if (forward_reverse, base_map_info['refbasei'][aligni]) in cgpos[0]:
            mfeatures[cur_row_num][1] = 0; mfeatures[cur_row_num][2] = 1
         elif (forward_reverse, base_map_info['refbasei'][aligni]) not in cgpos[1]:
            mfeatures[cur_row_num][1] = 1; mfeatures[cur_row_num][2] = 0
         else: mfeatures[cur_row_num][1] = 0; mfeatures[cur_row_num][2] = 0
         if not base_map_info['refbase'][aligni]=='-':
            if forward_reverse=='+': align_ref_pos += 1
            else: align_ref_pos -= 1
         aligni += 1

      if ie>=0 and ie<len(modevents) and moptions['fnum']==53:
         for currs in sp_param['f5data'][readk][2][modevents['start'][ie]:int(modevents['start'][ie]+int(modevents['length'][ie]+0.5))]:
             if currs>10 or currs<-10: print ('Error raw signal', currs, ie, modevents['start'][ie], modevents['length'][ie])
             curbin = int((currs+5)/binlen)
             if curbin<0: curbin = 0
             elif not curbin<binnum: curbin = binnum-1
             mfeatures[cur_row_num][curbin+3] += 1
      if ie>=0 and ie<len(modevents):
         cur_index_add = moptions['fnum'] - 3 + 3
         mfeatures[cur_row_num][cur_index_add + 0] = modevents["mean"][ie]
         mfeatures[cur_row_num][cur_index_add + 1] = modevents["stdv"][ie]
         mfeatures[cur_row_num][cur_index_add + 2] = modevents["length"][ie]
   return (mfeatures, isdif)


#

def getFeature_handler(moptions, h5files_Q, failed_Q):
   while not h5files_Q.empty():
      try:
         f5files, ctfolderid = h5files_Q.get(block=False)
      except:
         break; 

      sp_options = defaultdict();
      sp_options['ctfolder'] = moptions['outFolder']+str(ctfolderid)
      if not os.path.isdir(sp_options['ctfolder']):
         os.system('mkdir '+sp_options['ctfolder'])

      mGetFeature1(moptions, sp_options, f5files)

      for errtype, errfiles in sp_options["Error"].items():
         failed_Q.put((errtype, errfiles));
  
def mGetFeature_manager(moptions):
   #moptions['testrnn'] = True; #moptions['testrnn'] = False;
   #moptions['alignStr'] = 'bwa'
   
   start_time = time.time();

   f5files = glob.glob(os.path.join(moptions['wrkBase'],"*.fast5" ))
   if moptions['recursive']==1:
      f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*.fast5" )))
      f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*/*.fast5" )))
      f5files.extend(glob.glob(os.path.join(moptions['wrkBase'],"*/*/*/*.fast5" )))

   #if moptions['recursive']==1:
   #   f5files = glob.glob(os.path.join(moptions['wrkBase'],"*/*.fast5" ))
   #else: f5files = glob.glob(os.path.join(moptions['wrkBase'],"*.fast5" ))
   # f5files = ['/scr1/users/liuq1/project/nanopore/nanodeepmod/na12878albacore126/chr7/09_albacore/workspace/24/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_88738_ch24_read125_strand.fast5']

   '''
   allmchrs = []
   for ami in range(1, 23):
      allmchrs.append('chr'+str(ami))
   allmchrs.append('chrX'); allmchrs.append('chrY');
   moptions['fulmeth'] = {mchr:defaultdict() for mchr in allmchrs};
   moptions['anymeth'] = {mchr:defaultdict() for mchr in allmchrs};
   bs_folder = '/home/liuq1/project/nanopore/NanoDeepMod_v0.1.1/bisultfiteseq/'
   target_chr = 'chr7'; fulper = '0.95'
   fumfiles = [bs_folder+target_chr+'_CpG_'+fulper+'.txt', \
               bs_folder+target_chr+'_CHG_'+fulper+'.txt', \
               bs_folder+target_chr+'_CHH_'+fulper+'.txt',] 
   anyfulfiles = [bs_folder+target_chr+'_any_'+fulper+'.txt']
   mthreadin = [moptions['fulmeth'], moptions['anymeth']]
   mthfiles = [fumfiles, anyfulfiles]
   for mthi in range(len(mthreadin)):
      curmeth = mthreadin[mthi]; curfilelist = mthfiles[mthi]
      for curmthf in curfilelist:
         with open(curmthf, 'r') as mreader:
             line = mreader.readline();
             while line:
                 if len(line)>0:
                    tchr, tstrand, tpos = line.split()[:3]
                    curmeth[tchr][(tstrand, int(tpos))] = True;
                 line = mreader.readline();
   for tchr in allmchrs:
      if len(moptions['fulmeth'][tchr])>0 or len(moptions['anymeth'][tchr])>0:
         print ('%s fulmeth=%d anymeth=%d' % (tch, len(moptions['fulmeth'][tchr]), len(moptions['anymeth'][tchr])))
   '''
   fuldfiles = [glob.glob(moptions["fulmod"])]; 
   anydfiles = [glob.glob(moptions["anymod"])]
   moptions['fulmodlist'] = defaultdict(lambda: defaultdict()); 
   moptions['anymodlist'] = defaultdict(lambda: defaultdict());
   mthreadin = [moptions['fulmodlist'], moptions['anymodlist']]
   mthfiles = [fuldfiles, anydfiles]
   for mthi in range(len(mthreadin)):
      curmeth = mthreadin[mthi]; curfilelist = mthfiles[mthi]
      for curmthf in curfilelist:
          with open(curmthf, 'r') as mreader:
              line = mreader.readline();
              while line:
                 if len(line)>0:
                    tchr, tstrand, tpos = line.split()[:3]
                    curmeth[tchr][(tstrand, int(tpos))] = True;
                 line = mreader.readline();
   for tchr in moptions['anymodlist']:
      if len(moptions['fulmodlist'][tchr])>0 or len(moptions['anymodlist'][tchr])>0:
         print ('%s fulmod=%d anymod=%d' % (tchr, len(moptions['fulmodlist'][tchr]), len(moptions['anymodlist'][tchr]))) 

   pmanager = multiprocessing.Manager();
   h5files_Q = pmanager.Queue();
   failed_Q = pmanager.Queue()

   h5_batch = []; h5batchind = 0;
   for f5f in f5files:
      h5_batch.append(f5f);
      if len(h5_batch)==moptions['files_per_thread']:
         h5files_Q.put((h5_batch, h5batchind))
         h5_batch = []; h5batchind += 1
   if len(h5_batch)>0:
      h5files_Q.put((h5_batch, h5batchind))

   share_var = (moptions, h5files_Q, failed_Q)
   handlers = []
   for hid in range(moptions['threads']):
      p = multiprocessing.Process(target=getFeature_handler, args=share_var);
      p.start();
      handlers.append(p);

   failed_files = defaultdict(list);
   while any(p.is_alive() for p in handlers):
      try:
         errk, fns = failed_Q.get(block=False);
         failed_files[errk].extend(fns)
      except:
         time.sleep(1);
         continue;

   if len(failed_files)>0:
      print ('Error information for different fast5 files:')
      #for errtype, errfiles in failed_files.iteritems():
      for errtype, errfiles in failed_files.items():
         print ('\t'+errtype, len(errfiles))

   end_time = time.time();
   print ("Total consuming time %d" % (end_time-start_time))

if __name__=='__main__':
      moptions = {}
      moptions['basecall_1d'] = ['Basecall_1D_000']
      moptions['basecall_2strand'] = 'BaseCalled_template'

      moptions['wrkBase'] = "/scr1/users/liuq1/project/nanopore/nanodeepmod/na12878albacore126/chr7/09_albacore/workspace/"
      moptions['recursive'] = 1

      moptions['outFolder'] = '/scr1/users/liuq1/project/nanopore/nanodeepmod/na12878albacore126/chr7/test'
      moptions['FileID'] = 'test'
      moptions['Ref'] = "/mnt/isilon/wang_lab/shared/ref_human/hg38/hg38.fa"
      #moptions['kmer_model_file'] = 'scripts/kmer_model/r9.4_450bps.nucleotide.5mer.template.model'

      moptions['outLevel'] = myCom.OUTPUT_WARNING
      moptions['outLevel'] = myCom.OUTPUT_INFO

      moptions['fnum'] = 53; 
      moptions['hidden'] = 100;
      moptions['windowsize'] = 21;
      moptions['windowsize'] = 3
  
      moptions['threads'] = 8
      moptions['threads'] = 1
      moptions['files_per_thread'] = 500

      moptions['size_per_batch'] = 7

      mGetFeature_manager(moptions)


