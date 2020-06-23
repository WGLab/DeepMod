#!/usr/bin/env python

import os, sys, time
from collections import defaultdict
import glob

import numpy as np

from scipy import stats

import locale
locale.setlocale(locale.LC_ALL, 'en_US')

import tensorflow as tf

batch_size = 4096

cov_thrd = 5

def readBed(bedfile, t_chr=None, t_start=None, t_end=None):
   print('read {}'.format(bedfile)); sys.stdout.flush()
   beddict = defaultdict()
   with open(bedfile, 'r') as bedreader:
      start_time = time.time();
      line = bedreader.readline();
      while True:
         line = bedreader.readline();
         if not line: break;

         line = line.strip();
         if len(line)>20:
            mchr, start_pos, end_pos, _, _, m_strand, _, _, _, true_cov, meth_perc = line.split()
            start_pos = int(start_pos)
            true_cov = int(true_cov)
            if true_cov < cov_thrd: continue;
            meth_perc = round(int(meth_perc)/100.0, 3)
            if (t_chr not in [None, mchr]) or (not ((t_start==None or start_pos>=t_start) and (t_end==None or start_pos<=t_end))):
                 continue;
            if true_cov==0: continue
            beddict[(mchr, m_strand, start_pos)] = meth_perc
   return beddict

def readpredmod(predmodf, preddict, t_chr=None, t_start=None, t_end=None, cgCposdict=None):
   print('read {}'.format(predmodf)); sys.stdout.flush()
   with open(predmodf, 'r') as mr:
      while True:
          line = mr.readline();
          if not line: break;
          line = line.strip();
          if len(line)>0:
             lsp = line.split();
             cur_chr = lsp[0];
             cur_pos = int(lsp[1]);
             cur_strand = lsp[5];

             if not (cgCposdict==None or (cur_chr, cur_strand, cur_pos) in cgCposdict): continue

             cur_cov = int(lsp[9]);
             cur_m_p = int(lsp[10]);
             cur_m_c = int(lsp[11]);
             if (t_chr not in [None, cur_chr]) or (not ((t_start==None or cur_pos>=t_start) and (t_end==None or cur_pos<=t_end))):
                continue;
             if cur_cov==0: continue;
 
             if (cur_chr, cur_strand, cur_pos) not in preddict:
                 preddict[(cur_chr, cur_strand, cur_pos)] = [cur_cov, round(cur_m_p/100.0, 3), cur_m_c, line]
             else:
                 print("Warning_duplicate {}".format(predmodf))
                 preddict[(cur_chr, cur_strand, cur_pos)][0] += cur_cov
                 preddict[(cur_chr, cur_strand, cur_pos)][2] += cur_m_c
                 if preddict[(cur_chr, cur_strand, cur_pos)][0]>0:
                    preddict[(cur_chr, cur_strand, cur_pos)][1] = round(preddict[(cur_chr, cur_strand, cur_pos)][2]/float(preddict[(cur_chr, cur_strand, cur_pos)][0]), 3)



pred_file = sys.argv[1]+'.%s.C.bed'
save_file = sys.argv[1]+'_clusterCpG.%s.C.bed'
gmotfolder = sys.argv[2]

mpat = 'Cg'; mposinpat=0
stposofinterest = None; edposofinterest = None;

nbsize = 25;
train_mod = 'DeepMod/train_mod/na12878_cluster_train_mod-keep_prob0.7-nb25-chr1/{}.cov{}.nb{}'.format(mpat, cov_thrd, nbsize)

chrkeys = []
for i in range(1, 23):
   chrkeys.append("chr%d" % i)
chrkeys.append("chrX")
chrkeys.append("chrY")
chrkeys.append("chrM")


new_saver = tf.train.import_meta_graph(train_mod+'.meta')
print(new_saver); sys.stdout.flush()
with tf.Session() as sess:
   print("restore model: {} {}".format(train_mod+'.meta', train_mod[:train_mod.rindex('/')+1]))
   print(new_saver.restore(sess,tf.train.latest_checkpoint(train_mod[:train_mod.rindex('/')+1]))); sys.stdout.flush()

   mgraph = tf.get_default_graph()
   output = mgraph.get_tensor_by_name('output:0')
   X = mgraph.get_tensor_by_name('X:0')
   keep_prob = mgraph.get_tensor_by_name('keep_prob:0')
 
   for chrofinterest in chrkeys:
      #read pred
      preddict = defaultdict()

      cur_cg_pos = '%s/motif_%s_C.bed' % (gmotfolder, chrofinterest)
      if not os.path.isfile(cur_cg_pos): 
         print("Warning_motif!!! no file {}".format(cur_cg_pos))
         continue;
      if not os.path.isfile(pred_file % chrofinterest):
         print("Warning_pred!!! no file {}".format(pred_file % chrofinterest))
         continue;
 
      cgposdict = defaultdict();
      with open(cur_cg_pos, 'r') as mr:
         while True:
            line = mr.readline();
            if not line: break;
            lsp = line.split();
            cgposdict[ (lsp[0], lsp[2], int(lsp[1]) ) ] = True
      print("{}: read {} done! ".format(len(cgposdict), cur_cg_pos)); sys.stdout.flush()
      readpredmod(pred_file % chrofinterest, preddict, chrofinterest, cgCposdict=cgposdict)
      print("size={} vs ".format(len(preddict), len(cgposdict) )); sys.stdout.flush()

      train_data = []
      pdkeys = sorted(list( preddict.keys() ))
      for cspk in pdkeys: # preddict:
         if cspk not in cgposdict: 
            print("not in cpg warning!!! {} {}".format(chrofinterest, cspk))

         partner_pos = (cspk[0], '-' if cspk[1]=='+' else '+', cspk[2]+1 if cspk[1]=='+' else cspk[2]-1)
         cur_x = [preddict[cspk][1], preddict[partner_pos][1] if partner_pos in preddict else 0]
         for pdis in range(11):
            cur_x.append(0)
         cur_x.append(0)
         if len(train_data)<10: print("test")
         for rpos in range(cspk[2]-nbsize, cspk[2]+nbsize+1):
            if rpos in [cspk[2], partner_pos[2]]: continue;
            
            if (cspk[0], '+', rpos) in cgposdict and (cspk[0], '+', rpos) in preddict:
                cur_x[int(preddict[(cspk[0], '+', rpos)][1]/0.1+0.5) + 3] += 1
                cur_x[2] += 1
                if len(train_data)<10: print("\t\t{}: {}".format((cspk[0], '+', rpos), preddict[(cspk[0], '+', rpos)]))
            elif (cspk[0], '-', rpos) in cgposdict and (cspk[0], '-', rpos) in preddict:
                cur_x[int(preddict[(cspk[0], '-', rpos)][1]/0.1+0.5) + 3] += 1
                cur_x[2] += 1
                if len(train_data)<10: print("\t\t{}: {}".format((cspk[0], '-', rpos), preddict[(cspk[0], '-', rpos)]))
         for i in range(3, len(cur_x)):
            if cur_x[2]>0: cur_x[i] = round(cur_x[i]/float(cur_x[2]), 3)
         if len(train_data)<10: print('\t{}'.format(cur_x)); sys.stdout.flush()
         train_data.append(cur_x)

      print("format data: data={}; {}".format(len(train_data), len(train_data[0]))); sys.stdout.flush()
      
      batch_data = np.array_split(train_data, int(len(train_data)/batch_size) if len(train_data)>batch_size else 2)
      m_pred_new_per = []
      for i in range(len(batch_data)):
          moutp = sess.run([output], feed_dict={X:batch_data[i], keep_prob:1})
          for mpind in moutp:
              for curpd in mpind:
                 m_pred_new_per.append(curpd)
      print("new per: {}, {}  {} {}".format(len(pdkeys), len(train_data), len(m_pred_new_per), curpd ))
      for wind in range(10):
         print("'{}' <{}> {}".format(preddict[pdkeys[wind]][-1], m_pred_new_per[wind], train_data[wind]))
      with open(save_file % chrofinterest, 'w') as mwriter:
         for wind in range(len(pdkeys)):
            mwriter.write("{} {}\n".format(preddict[pdkeys[wind]][-1], int(m_pred_new_per[wind]*100)))
 
 
