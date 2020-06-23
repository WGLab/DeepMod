#!/usr/bin/env python

import os, sys, string
import multiprocessing
import time;
import glob
from collections import defaultdict;


def mprint(mstr):
   print(mstr); sys.stdout.flush()

if len(sys.argv)<4:
   print ("Usage: python {} pred_folder-of-DeepMod Base-of-interest unique-fileid-in-sum-file [chr-list]".format(sys.argv[0]))
   print ("       pred_folder-of-DeepMod: the prediction must in its sub-folder.")
   sys.exit(1)

pred_folder = sys.argv[1]
baseofint = sys.argv[2]
sum_fileid = sys.argv[3]


if len(sys.argv)>4:
   chrkeys = ["%s" % cid for cid in sys.argv[4].split(',')]
else:
   chrkeys = []
   for i in range(1, 23):
      chrkeys.append("chr%d" % i)
   chrkeys.append("chrX")
   chrkeys.append("chrY")
   chrkeys.append("chrM")

chrkeys = set(chrkeys)

#####################################
def readbed2(bedf):
   cur_ad = defaultdict();
   with open(bedf, 'r') as mr:
      while True:
         line = mr.readline();
         if not line: break;
         line = line.strip();
         lsp = line.split();
         cur_ad[(lsp[0], int(lsp[1]), lsp[5])] = [int(lsp[9]), int(lsp[11])]
   return cur_ad

def mergeMod(g_ad, cur_ad):
    for k in cur_ad:
       if k in g_ad:
            g_ad[k][0] += cur_ad[k][0]
            g_ad[k][1] += cur_ad[k][1]
       else: g_ad[k] = cur_ad[k]
      
def save_mod(res_file, amod_dict, baseOfin):
   poskeys = sorted(list(amod_dict.keys()))
   for pk in poskeys:
      if amod_dict[pk][1]==0: del amod_dict[pk]

   poskeys = list(amod_dict.keys())
   poskeys = sorted(poskeys);
   with open(res_file, 'w') as mw:
      for pk in poskeys:
          mw.write('%s %d %d %s %d %s  %d %d 0,0,0 %d %d %d\n' % (pk[0], pk[1],pk[1]+1, baseOfin,amod_dict[pk][0] if amod_dict[pk][0]<1000 else 1000,pk[2], pk[1],pk[1]+1, amod_dict[pk][0], int(amod_dict[pk][1]*100/amod_dict[pk][0]) if amod_dict[pk][0]>0 else 0, amod_dict[pk][1]  ))
       

def sum_amod_handler(run_Q):
   handli = 0;
   while not run_Q.empty():
      try:
          ck, pred_folder, baseOfin = run_Q.get(block=False)
          #print (ck, pred_folder, baseOfin); continue
      except:
         break;

      allbedfiles = glob.glob(os.path.join(pred_folder, ("*/*/*/*.%s-.%s.bed" % (ck, baseOfin)) ))
      allbedfiles.extend(glob.glob(os.path.join(pred_folder, ("*/*/*.%s-.%s.bed" % (ck, baseOfin)) )))
      allbedfiles.extend(glob.glob(os.path.join(pred_folder, ("*/*.%s-.%s.bed" % (ck, baseOfin)) )))
      mprint ("%s - %s: %d" % (ck, baseOfin, len(allbedfiles) )) 
      allbedfiles.extend(glob.glob(os.path.join(pred_folder, ("*/*/*/*.%s+.%s.bed" % (ck, baseOfin)) )))
      allbedfiles.extend(glob.glob(os.path.join(pred_folder, ("*/*/*.%s+.%s.bed" % (ck, baseOfin)) )))
      allbedfiles.extend(glob.glob(os.path.join(pred_folder, ("*/*.%s+.%s.bed" % (ck, baseOfin)) )))
      mprint ("%s -+ %s: %d" % (ck, baseOfin, len(allbedfiles) ))

      #  0    1     2     3 4 5   6      7      8   9 0 1
      #chr1 949802 949803 T 1 - 949802 949803 0,0,0 1 0 0 
      amod_dict = defaultdict();      
      res_file = "%s/%s.%s.%s.bed" % (pred_folder, sum_fileid, ck, baseOfin)
      for bedf_ind in range(len(allbedfiles)):
          mprint("\t %s %s %d/%d" % (ck, baseOfin, bedf_ind+1, len(allbedfiles)))
          cur_ad = readbed2(allbedfiles[bedf_ind])
          mergeMod(amod_dict, cur_ad)

      save_mod(res_file, amod_dict, baseOfin)

##############################
pmanager = multiprocessing.Manager();
run_Q = pmanager.Queue();
for ck in chrkeys:
   run_Q.put((ck, pred_folder, baseofint))

mhandlers = [];
share_var = (run_Q, )
m_thread_num = len(chrkeys);
for i in range(m_thread_num+1):
   p = multiprocessing.Process(target=sum_amod_handler, args=share_var)
   p.start();
   mhandlers.append(p);
while any(p.is_alive() for p in mhandlers):
   try:
      time.sleep(1);
   except:
      time.sleep(1);
      continue;





