
import os, sys


#                                             0  1  2  3  4     5   6  7  8  9
cmd_temp = "nohup python DeepMod/tools/cal_EcoliDetPerf.py %s %s %s %d  '%s' %d %d %s %s &"

base_ecoli_res_folder = 'ecoli_pred/'

refbase = 'ref/'

mrun_param = [\
              [base_ecoli_res_folder+'gAtc', \
               refbase + 'Ecoli_k12_mg1655.fasta', \
               'gAtc', 1, '', 1000000, 2000000, base_ecoli_res_folder+'gAtc/', base_ecoli_res_folder+'con1a;'+base_ecoli_res_folder+'con2a'    ], \
              [base_ecoli_res_folder+'tcgA', \
               refbase + 'Ecoli_k12_mg1655.fasta', \
               'tcgA', 3, '', 1000000, 2000000, base_ecoli_res_folder+'tcgA/', base_ecoli_res_folder+'con1a;'+base_ecoli_res_folder+'con2a'     ], \
              [base_ecoli_res_folder+'gaAttc', \
               refbase + 'Ecoli_k12_mg1655.fasta', \
               'gaAttc', 2, '', 1000000, 2000000, base_ecoli_res_folder+'gaAttc/', base_ecoli_res_folder+'con1a;'+base_ecoli_res_folder+'con2a' ], \

              [base_ecoli_res_folder+'Cgmpe',\
               refbase + 'Ecoli_k12_mg1655.fasta', \
               'Cg', 0, '', -1, -1, base_ecoli_res_folder+'Cgmpe/', base_ecoli_res_folder+'con1;'+base_ecoli_res_folder+'con2'  ],\
              [base_ecoli_res_folder+'Cgsss', \
               refbase + 'Ecoli_k12_mg1655.fasta', \
               'Cg', 0, '', -1, -1, base_ecoli_res_folder+'Cgsss/', base_ecoli_res_folder+'con1;'+base_ecoli_res_folder+'con2'  ], \
              [base_ecoli_res_folder+'gCgc',\
               refbase + 'Ecoli_k12_mg1655.fasta', \
               'gCgc', 1, '', -1, -1, base_ecoli_res_folder+'gCgc/', base_ecoli_res_folder+'con1;'+base_ecoli_res_folder+'con2'   ],\
             ]



for curRunP in mrun_param:
   print(curRunP[0])
   #if not os.path.isdir(curRunP[7]):
   #   os.system('mkdir -p '+curRunP[7])
   cur_cmd = (cmd_temp % (tuple(curRunP)))
   print ('\t'+cur_cmd)
   if len(sys.argv)>2: 
      os.system(cur_cmd)
      #print ()
