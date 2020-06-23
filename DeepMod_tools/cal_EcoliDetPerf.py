#!/usr/bin/env python

import os, sys, time
from collections import defaultdict
import glob
import copy
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.metrics import matthews_corrcoef

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from pkg_resources import resource_string

from scipy.stats import binom
import copy

ggplot = importr('ggplot2')
importr('gridExtra')
importr('plyr')

na4com = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

def readFA(mfa, mpat='Cg', mposinpat=0, t_chr=None, t_start=None, t_end=None):
   pos_dict = defaultdict(int)
  
   pat3 = copy.deepcopy(mpat.upper())
   comp_pat3 = ''.join([na4com[curna] for curna in pat3][::-1])
   comp_mposinpat = len(comp_pat3)-1-mposinpat
 
   fadict = defaultdict();
   with open(mfa, 'r') as mr:
      cur_chr = None;

      line = mr.readline();
      while line:
         line = line.strip();
         if len(line)>0:
            if line[0]=='>': 
               if not cur_chr==None:
                  fadict[cur_chr] = ''.join(fadict[cur_chr])
               cur_chr = line[1:].split()[0]
               if t_chr in [None, cur_chr]:
                  fadict[cur_chr] = []
            else:
               if t_chr in [None, cur_chr]: 
                  fadict[cur_chr].append(line)
         line = mr.readline();
      if not cur_chr==None:
         fadict[cur_chr] = ''.join(fadict[cur_chr]) 
   fakeys = fadict.keys();
   cpgdict = defaultdict(int); cpgnum = [0, 0]
   for fak in fakeys:
       cpgdict[fak] = defaultdict()
       for i in range(len(fadict[fak])):
          if (t_start==None or i>=t_start) and (t_end==None or i<=t_end):
             if i-mposinpat>=0 and i+len(comp_pat3)-1-mposinpat<len(fadict[fak]) and ''.join(fadict[fak][i-mposinpat:(i+len(comp_pat3)-1-mposinpat+1)])==pat3:
                cpgdict[fak][('+', i)] = [1, fadict[fak][i]]; cpgnum[0] += 1
                cpgdict[fak][('-', i)] = [0, fadict[fak][i]]
             elif i-comp_mposinpat>=0 and i+len(comp_pat3)-1-comp_mposinpat<len(fadict[fak]) and ''.join(fadict[fak][i-comp_mposinpat:(i+len(comp_pat3)-1-comp_mposinpat+1)])==comp_pat3:
                cpgdict[fak][('+', i)] = [0, fadict[fak][i]]
                cpgdict[fak][('-', i)] = [1, fadict[fak][i]]; cpgnum[1] += 1
             else:
                cpgdict[fak][('+', i)] = [0, fadict[fak][i]]
                cpgdict[fak][('-', i)] = [0, fadict[fak][i]]
   print('%s%d site: %d(+) %d(-)' % (pat3, mposinpat, cpgnum[0], cpgnum[1]))
   return cpgdict

def readmodf_dict(cpgdict, modf, pred_dict, mna, t_start=None, t_end=None):
   with open(modf, 'r') as mr:
      while True:
          line = mr.readline();
          if not line: break;
          line = line.strip();
          if len(line)>0:
             lsp = line.split();
             cur_chr = lsp[0];
             cur_pos = int(lsp[1]);
             cur_strand = lsp[5];

             cur_cov = int(lsp[9]);
             cur_m_p = int(lsp[10]);
             cur_m_c = int(lsp[11]);

             if not ((t_start==None or cur_pos>=t_start) and (t_end==None or cur_pos<=t_end)):
                line = mr.readline();
                continue;

             if not (mna==lsp[3] and lsp[3]==(cpgdict[cur_chr][(cur_strand, cur_pos)][1] if cur_strand=='+' else na4com[cpgdict[cur_chr][(cur_strand, cur_pos)][1]])):
                print ('Error !! NA not equal %s == %s == %s %s' % (mna, lsp[3], cpgdict[cur_chr][(cur_strand, cur_pos)][1], modf))

             if (cur_chr, cur_pos, cur_strand) not in pred_dict: 
                pred_dict[(cur_chr, cur_pos, cur_strand)] = [cur_cov, cur_m_p, cur_m_c, lsp[3]]
             else:
                pred_dict[(cur_chr, cur_pos, cur_strand)][0] += cur_cov
                pred_dict[(cur_chr, cur_pos, cur_strand)][2] += cur_m_c
                pred_dict[(cur_chr, cur_pos, cur_strand)][1] = int(pred_dict[(cur_chr, cur_pos, cur_strand)][2]*100/pred_dict[(cur_chr, cur_pos, cur_strand)][0]) if pred_dict[(cur_chr, cur_pos, cur_strand)][0]>0 else 0

def add_from_dict(cpgdict, pred_dict, label, pred_list, mna, tp_fp_tn_fn, mpat='Cg', mposinpat=0):
   for posk in pred_dict:
             cur_chr, cur_pos, cur_strand = posk
             cur_cov, cur_m_p, cur_m_c, lsp3 = pred_dict[posk]

             iscpg = False;
             if cpgdict[cur_chr][(cur_strand, cur_pos)][0]==1:
                 iscpg = True;
                 pred_list.append((label, cur_cov, cur_m_p, cur_m_c, mpat, np.log(binom.pmf(cur_m_c, cur_cov, 0.05)) ))
                 if (lsp3==mpat[mposinpat]): pass
                 else: print ('Error not methylated pos %s %s %s' % (mna, cur_strand))
             if not iscpg:
                isclosec = False;
                for i in range(-3, 4):
                   if (cur_strand, cur_pos+i) in cpgdict[cur_chr] and cpgdict[cur_chr][(cur_strand, cur_pos+i)][0]==1:
                      isclosec = True; break;
                if lsp3==mpat[mposinpat]:
                   pred_list.append((0, cur_cov, cur_m_p, cur_m_c, mpat+'_n'+str(abs(i))+mpat[mposinpat] if isclosec else 'Other'+mpat[mposinpat], np.log(binom.pmf(cur_m_c, cur_cov, 0.05)) ))
                else:
                   pred_list.append((0, cur_cov, cur_m_p, cur_m_c, mpat+'_nb' if isclosec else 'Other', np.log(binom.pmf(cur_m_c, cur_cov, 0.05)) ))
             if pred_list[-1][0]==0:
                tp_fp_tn_fn[2] += cur_cov - cur_m_c
                tp_fp_tn_fn[1] += cur_m_c
             else:
                tp_fp_tn_fn[0] += cur_m_c
                tp_fp_tn_fn[3] += cur_cov - cur_m_c


def readmodf(cpgdict, modf, label, pred_list, mna, tp_fp_tn_fn, mpat='Cg', mposinpat=0, t_start=None, t_end=None):
   with open(modf, 'r') as mr:
      line = mr.readline();
      while line:
          line = line.strip();
          if len(line)>0:
             lsp = line.split();
             cur_chr = lsp[0];
             cur_pos = int(lsp[1]);
             cur_strand = lsp[5];
             
             cur_cov = int(lsp[9]);
             cur_m_p = int(lsp[10]);
             cur_m_c = int(lsp[11]); 

             if not ((t_start==None or cur_pos>=t_start) and (t_end==None or cur_pos<=t_end)): 
                line = mr.readline();
                continue;

             if not (mna==lsp[3] and lsp[3]==(cpgdict[cur_chr][(cur_strand, cur_pos)][1] if cur_strand=='+' else na4com[cpgdict[cur_chr][(cur_strand, cur_pos)][1]])):
                print ('Error !! NA not equal %s == %s == %s %s' % (mna, lsp[3], cpgdict[cur_chr][(cur_strand, cur_pos)][1], modf))
             iscpg = False;
             if cpgdict[cur_chr][(cur_strand, cur_pos)][0]==1:
                 iscpg = True;
                 pred_list.append((label, cur_cov, cur_m_p, cur_m_c, mpat, np.log(binom.pmf(cur_m_c, cur_cov, 0.05)) ))
                 if (lsp[3]==mpat[mposinpat]): pass
                 else: print ('Error not methylated pos %s %s %s' % (mna, cur_strand, modf))
             if not iscpg:
                isclosec = False;
                for i in range(-3, 4):
                   if (cur_strand, cur_pos+i) in cpgdict[cur_chr] and cpgdict[cur_chr][(cur_strand, cur_pos+i)][0]==1:
                      isclosec = True; break;
                if lsp[3]==mpat[mposinpat]:
                   pred_list.append((0, cur_cov, cur_m_p, cur_m_c, mpat+'_n'+str(abs(i))+mpat[mposinpat] if isclosec else 'Other'+mpat[mposinpat], np.log(binom.pmf(cur_m_c, cur_cov, 0.05)) ))
                else:
                   pred_list.append((0, cur_cov, cur_m_p, cur_m_c, mpat+'_nb' if isclosec else 'Other', np.log(binom.pmf(cur_m_c, cur_cov, 0.05)) ))
             if pred_list[-1][0]==0:
                tp_fp_tn_fn[2] += cur_cov - cur_m_c
                tp_fp_tn_fn[1] += cur_m_c
             else:
                tp_fp_tn_fn[0] += cur_m_c
                tp_fp_tn_fn[3] += cur_cov - cur_m_c 
          line = mr.readline();   


sssfolder = sys.argv[1];   # 
mreffile = sys.argv[2];    # 
mpat=sys.argv[3];          # Cg
mposinpat=int(sys.argv[4]);# 0 

chrofinterest = sys.argv[5];
if chrofinterest=='': chrofinterest = None;
stposofinterest = int(sys.argv[6]);
if stposofinterest<0: stposofinterest = None;
edposofinterest = int(sys.argv[7]);
if edposofinterest<0: edposofinterest = None;

basefig = sys.argv[8]
hastwoclass = 1;

sssfiles = {mpat[mposinpat]:glob.glob(os.path.join(sssfolder, 'mod_pos.*.'+mpat[mposinpat]+'.bed'))}
sssfiles[mpat[mposinpat]].extend(glob.glob(os.path.join(sssfolder, '*/mod_pos.*.'+mpat[mposinpat]+'.bed')))
sssfiles[mpat[mposinpat]].extend(glob.glob(os.path.join(sssfolder, '*/*/mod_pos.*.'+mpat[mposinpat]+'.bed')))
print(str(len(sssfiles[mpat[mposinpat]])) + " " + str(sssfolder))

## for negative;
umrfiles = []
for cur_umr_f in sys.argv[9].split(','):
   if not os.path.isdir(cur_umr_f):
       print("No prediction folder {}".format(cur_umr_f))
       sys.exit(1);
   umrfiles.extend(glob.glob(os.path.join(cur_umr_f, '*/*/mod_pos.*.'+mpat[mposinpat]+'.bed')))
   umrfiles.extend(glob.glob(os.path.join(cur_umr_f, '*/mod_pos.*.'+mpat[mposinpat]+'.bed')))
   umrfiles.extend(glob.glob(os.path.join(cur_umr_f, 'mod_pos.*.'+mpat[mposinpat]+'.bed')))
print(str(len(umrfiles)) + "  " + str(sys.argv[9].split(',')))
sys.stdout.flush()

for sa in sssfiles:
   print (sa)
   for nf in sssfiles[sa]:
      print ('\t'+nf)

cpgdict = readFA(mreffile, mpat, mposinpat, chrofinterest, stposofinterest, edposofinterest)

pred_dict = defaultdict();
for modf in umrfiles:
   readmodf_dict(cpgdict, modf, pred_dict, mpat[mposinpat], stposofinterest, edposofinterest)

baseinfo = [mpat, mpat+'_n1'+mpat[mposinpat], mpat+'_n2'+mpat[mposinpat], mpat+'_n3'+mpat[mposinpat], 'Other'+mpat[mposinpat], mpat+'_nb', 'Other']

classify_m = ['Methylation_Percentage']
classify_types = [baseinfo, [mpat]]
filename = [['all_mp','motif_mp'] ]
cov_thr = [1, 5]
mlinestyle = {1:'bo-', 3:'gx--', 5:'r*-.', 7:'cs-', 10:'md--', 15:'k+-.'}

pred_list = []; tp_fp_tn_fn = [0, 0, 0, 0]

add_from_dict(cpgdict, pred_dict, 0, pred_list, mpat[mposinpat], tp_fp_tn_fn, mpat, mposinpat)

if True:
   for na4 in sssfiles:
      for cur_f in sssfiles[na4]:
         print('%s %s' % (na4, cur_f)); sys.stdout.flush();
         readmodf(cpgdict, cur_f, hastwoclass, pred_list, na4, tp_fp_tn_fn, mpat, mposinpat, stposofinterest, edposofinterest);
   pred_list = np.array(pred_list, dtype=[('Methylation', np.uint), ('Coverage', np.uint64), ('Methylation_Percentage', np.uint64), ('Methylation_Coverage', np.uint64), ('BaseInfo', 'U20'), ('logp', np.float64)])
   
   if hastwoclass==1:
      cov_plot_thr = [1, 5]
      for ct_ind in range(len(classify_types)):
         ct = classify_types[ct_ind]
         cur_ct_data = pred_list[np.isin(pred_list['BaseInfo'], ct)]
         for cm_ind in range(len(classify_m)):
             print('basetype={} classify_measure={}'.format(ct, classify_m[cm_ind]))
             cm = classify_m[cm_ind]
           
             # 1 for roc, 2: pr;
             roc_or_pr = 2; roc_or_pr=0
             for roc_or_pr in range(1,3):
              if roc_or_pr>0: 
                mfig= plt.figure()
              if roc_or_pr==2:
                cur_fn = basefig+'/ap_plot_met_pr_'+filename[cm_ind][ct_ind]+'.png'
                xylab = ['Recall', 'Precision'];  leg_mpos = "lower left"
                for covt in cov_plot_thr:
                   precision, recall, thresholds = precision_recall_curve(cur_ct_data['Methylation'][cur_ct_data['Coverage']>=covt], cur_ct_data[cm][cur_ct_data['Coverage']>=covt])
                   ap_pr = average_precision_score(cur_ct_data['Methylation'][cur_ct_data['Coverage']>=covt], cur_ct_data[cm][cur_ct_data['Coverage']>=covt])
                   plt.plot(recall, precision, mlinestyle[covt], lw=2, label='Coverage>=%d (AP=%0.3f)' % (covt, ap_pr))
                   print('\t\t %s %d ap=%.5f' % (cur_fn, covt, ap_pr))
              elif roc_or_pr==1:
                xylab = ['False Positive Rate', 'True Positive Rate']; leg_mpos = "lower right"
                cur_fn = basefig+'/roc_plot_met_roc_'+filename[cm_ind][ct_ind]+'.png'
                prev = 0; prev_ind = -1
                for covt in cov_plot_thr:
                   fpr, tpr, mthr = roc_curve(cur_ct_data['Methylation'][cur_ct_data['Coverage']>=covt], cur_ct_data[cm][cur_ct_data['Coverage']>=covt])
                   #print(','.join([str(np.round(t1, 5)) for t1 in mthr]))
                   roc_auc = auc(fpr, tpr)
                   if (not np.isnan(roc_auc)) and (abs(roc_auc - prev)>=0.02 or (covt>10 and abs(roc_auc - prev)>=0.005) or (cov_plot_thr.index(covt)-prev_ind>1 and abs(roc_auc - prev)>=0.005)):
                      plt.plot(fpr, tpr, mlinestyle[covt], lw=2, label='Coverage>=%d (AUC=%0.3f)' % (covt, roc_auc)) 
                      prev = roc_auc; prev_ind = cov_plot_thr.index(covt)
                   if not np.isnan(roc_auc):
                      print ('\t\t %s %d %.7f' % (cur_fn, covt, roc_auc))
                plt.plot([0, 1], [0, 1])
              if roc_or_pr>0:
                plt.xlim([0.0, 1.0]);              plt.ylim([0.0, 1.0])
                plt.xlabel(xylab[0]);              plt.ylabel(xylab[1])
                plt.legend(loc=leg_mpos)
                mfig.savefig(cur_fn, dpi=300);              plt.close(mfig)


