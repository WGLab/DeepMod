#!/usr/bin/env python

import os, sys, string, time
from collections import defaultdict;
import multiprocessing



def read_genome(mfafile):
   ref_genome = defaultdict();
   with open(mfafile, 'r') as mr:
      cur_chr = None;
      while True:
         line = mr.readline();
         if not line: break;
         line = line.strip();
         if len(line)==0: continue;
         if line[0]=='>':
            if not cur_chr==None:
               ref_genome[cur_chr] = ''.join(seqlist);
            cur_chr = line[1:].split()[0];
            seqlist = []
         else:
            seqlist.append(line.upper())
      ref_genome[cur_chr] = ''.join(seqlist);
   print("Total chr: {}".format(len(ref_genome))); sys.stdout.flush()
   return ref_genome
#ref_genome = read_genome(ref_fa);

def handle_motif_pos(run_Q):
   handli = 0;
   while not run_Q.empty():
      try:
         rgkey, ref_genome, res_folder, na_bp, curna, curmotif = run_Q.get(block=False)
         #print(rgkey, ref_genome, res_folder, na_bp, curna, curmotif); continue
      except:
         break;

      #curna_dict = defaultdict();
      #curmotif_dict = defaultdict();
      nafile = '%sna_%s_%s.bed' % (res_folder, rgkey, curna)
      motiffile = '%smotif_%s_%s.bed' % (res_folder, rgkey, curna)
      mw_na = open(nafile, 'w')
      mw_motif = open(motiffile, 'w')

      batchsize = 500000
      print("get motif for {}={}".format(rgkey, len(ref_genome))); sys.stdout.flush()
      cur_hi = 0; start_time = time.time();
      for na_ind in range(len(ref_genome)):
          cur_hi += 1
          if cur_hi % batchsize == 0:
             print('\t time consuming ({})= {} {}'.format( rgkey, cur_hi, time.time() - start_time ) )
             sys.stdout.flush()
             start_time = time.time()

          if (ref_genome[na_ind]==curna or na_bp[ref_genome[na_ind]]==curna):
             #curna_dict[(rgkey, na_ind, '+' if ref_genome[rgkey][na_ind]==curna else '-')] = True;
             mw_na.write('%s\t%s\t%s\n' % (rgkey, na_ind, '+' if ref_genome[na_ind]==curna else '-'))
          if ref_genome[na_ind]==curna and (not curmotif==None):
             for cur_mot in curmotif:
                 is_mot = True; mot_ind = 0;
                 for cur_sub_r_ind in range(na_ind - curmotif[cur_mot], na_ind + len(cur_mot) - curmotif[cur_mot] ):
                     if cur_sub_r_ind<0 or cur_sub_r_ind>len(ref_genome)-1: 
                        is_mot = False; break;
                     if not ref_genome[cur_sub_r_ind] == cur_mot[mot_ind]: 
                        is_mot = False; break;
                     mot_ind += 1
                 if is_mot:
                    mw_motif.write('%s\t%s\t%s\n' % (rgkey, na_ind, '+'))
                    mw_motif.write('%s\t%s\t%s\n' % (rgkey, na_ind+1, '-'))
                    break;
      mw_na.close(); mw_motif.close()

                   

ref_fa = 'ref/hg38.fa'
ref_fa = sys.argv[1]
res_folder = 'genome.motif/C/'
res_folder = sys.argv[2]+'/'
if not os.path.isdir(res_folder):
   os.system('mkdir -p '+res_folder)

curna='C'
curmotif={'CG':0}
curna = sys.argv[3];
curmotif = {sys.argv[4]:int(sys.argv[5])}

if len(sys.argv)>6:
   chrkeys = ["chr%s" % cid for cid in sys.argv[6].split(',')]
else:
   chrkeys = []
   for i in range(1, 23):
      chrkeys.append("chr%d" % i)
   chrkeys.append("chrX")
   chrkeys.append("chrY")
   chrkeys.append("chrM")

chrkeys = set(chrkeys)


na_bp = {"A":"T", \
         "C":"G", \
         "G":"C", \
         "T":"A", \
         "a":"t", \
         "c":"g", \
         "g":"c", \
         "t":"a", \
         "N":"N", \
         "n":"n" \
         }


ref_genome = read_genome(ref_fa);

##############################
pmanager = multiprocessing.Manager();
run_Q = pmanager.Queue();
for curk in chrkeys:
   run_Q.put((curk, ref_genome[curk], res_folder, na_bp, curna, curmotif))

mhandlers = [];
share_var = (run_Q, )
m_thread_num = len(chrkeys);
for i in range(m_thread_num):
   p = multiprocessing.Process(target=handle_motif_pos, args=share_var)
   p.start();
   mhandlers.append(p);
while any(p.is_alive() for p in mhandlers):
   try:
      time.sleep(1);
   except:
      time.sleep(1);
      continue;






