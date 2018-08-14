
import os;
import sys;

import string;

from collections import defaultdict

import argparse;
from argparse import RawTextHelpFormatter

from scripts.myCom import *

parser = argparse.ArgumentParser(description="Detect nucleotide modification from nanopore signals data.", epilog="For example, \n \
\tpython %(prog)s train: Training a modification classifier.\n \
\tpython %(prog)s predict: Predict modification positions for each long reads in FAST5 files.\n \
\tpython %(prog)s detect: Detect modification by integrating all long reads. \n \
", formatter_class=RawTextHelpFormatter);


def non_negative(i, mstr):
   if i<1: return (("\n\tError %d could not be negative(%d)" % (mstr, i)))
   else: return ''


def printParameters(moptions):
   mpkeys = moptions.keys(); #mpkeys.sort()
   sorted(mpkeys)
   print('%30s: %s' % ('Current directory', os.getcwd()))
   for mpk in mpkeys:
      print ('%30s: %s' % (mpk, str(moptions[mpk])))
   sys.stdout.flush()

def mCommonParam(margs):

   ErrorMessage = ""
   moptions = defaultdict()
   moptions['outLevel'] = margs.outLevel
   moptions["wrkBase"] = margs.wrkBase
   #if moptions["wrkBase"]==None or (not os.path.isdir(moptions["wrkBase"])):
   #   ErrorMessage = ErrorMessage + ("\n\tThe directory (%s) does not exist" % moptions["wrkBase"])

   moptions["FileID"] = margs.FileID
   moptions['outFolder'] = margs.outFolder
   moptions['outFolder'] = format_last_letter_of_folder(moptions['outFolder'])
   if moptions['outFolder']==None or (not os.path.isdir(moptions['outFolder'])):
      try:
         os.system('mkdir -p '+moptions['outFolder']);
      except: 
         ErrorMessage = ErrorMessage + ("\n\tThe output folder (%s) does not exist and cannot be created." % moptions['outFolder'])

   moptions['recursive'] = margs.recursive
   moptions['files_per_thread'] = margs.files_per_thread
   if moptions['files_per_thread']<2: moptions['files_per_thread'] = 2
   moptions['threads'] = margs.threads
   if moptions['threads']<1: moptions['threads'] = 1

   moptions['windowsize'] = margs.windowsize
   ErrorMessage = ErrorMessage + non_negative(moptions['windowsize'], 'windowsize')
   if moptions['windowsize']<1: moptions['windowsize'] = 1

   moptions['alignStr'] = margs.alignStr;

   return [moptions, ErrorMessage]

def mDetect(margs):
   moptions, ErrorMessage = mCommonParam(margs)

   moptions['basecall_1d'] = margs.basecall_1d
   moptions['basecall_2strand'] = margs.basecall_2strand
   moptions['ConUnk'] = margs.ConUnk

   moptions['predDet'] = margs.predDet
   if moptions['predDet']:
      moptions['Ref'] = margs.Ref
      if moptions['Ref']==None or (not os.path.isfile(moptions['Ref'])):
         ErrorMessage = ErrorMessage + ("\n\t reference file does not exist (%s)" % moptions['Ref'])

      moptions['fnum'] = margs.fnum
      ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
      moptions['hidden'] = margs.hidden
      ErrorMessage = ErrorMessage + non_negative(moptions['hidden'], 'hidden')
      moptions['modfile'] = margs.modfile
      if moptions['modfile']==None:
         print("No mod file is provided. The default one is used")
         moptions['modfile'] = ('train_mod/rnn_P90wd%d_f53/mod_train_P90wd%d_f53' % (moptions['windowsize'], moptions['windowsize']))
      if (not os.path.isfile(moptions['modfile']+'.meta')):
         ErrorMessage = ErrorMessage + ("\n\tThe meta file (%s) does not exist" % (moptions['modfile']+'.meta' if not moptions['modfile']==None else ""))
      #moptions['predpath'] = moptions['outFolder'] #+moptions['FileID']+'_mpred.txt'
   else:
      moptions['predpath'] = margs.predpath
      if moptions['predpath']==None or (not os.path.isdir(moptions['predpath'])):
         ErrorMessage = ErrorMessage + ("\n\tThe predpath does not exist")

   moptions['region'] = [ ]
   if margs.region == None or len(margs.region)==0:
      moptions['region'].append([None, None, None])
   else:
      mregionlist = margs.region.split(';')
      for mr in mregionlist:
         mr_sp = mr.split(':')
         moptions['region'].append([mr_sp[0], int(mr_sp[1]), int(mr_sp[2]) ])

   print("\nNanopore sequencing data analysis is resourece-intensive and time consuming. ")
   print("Some potential strong recommendations are below:")
   print("\tIf your reference genome is large as human genome and your Nanopore data is huge,")
   print("\tIt would be faster to run this program parallelly to speed up.")
   print("\tYou might run different input folders of your fast5 files and ")
   print("\tgive different output names (--FileID) or folders (--outFolder)")
   print("\tA good way for this is to run different chromosome individually.\n")

   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['detect', '--help']);
      sys.exit(1)
   
   #if moptions['predDet']:
   #   from scripts import myMultiBiRNN
   #   myMultiBiRNN.pred_entry(moptions)

   from scripts import myDetect
   myDetect.mDetect_manager(moptions)


def mPredict(margs):
   from scripts import myMultiBiRNN

   moptions, ErrorMessage = mCommonParam(margs)

   moptions['fnum'] = margs.fnum
   ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
   moptions['hidden'] = margs.hidden
   ErrorMessage = ErrorMessage + non_negative(moptions['hidden'], 'hidden')
   moptions['modfile'] = margs.modfile
   if moptions['modfile']==None or (not os.path.isfile(moptions['modfile']+'.meta')):
      ErrorMessage = ErrorMessage + ("\n\tThe meta file (%s) does not exist" % (moptions['modfile']+'.meta'))

   moptions['basecall_1d'] = margs.basecall_1d
   moptions['basecall_2strand'] = margs.basecall_2strand

   if not margs.test==None:
      moptions['test'] = margs.test.split(',')
      if moptions['test'][0] == 'I': moptions['test'][0] = '+'
      elif moptions['test'][0] == 'P': moptions['test'][0] = '0'
      else:
         ErrorMessage = ErrorMessage + "Unknown option for test: the first character must be I or P "+margs.test
      if moptions['test'][0] in ['+']:
         moptions['test'][1] = int(moptions['test'][1]) * (10**6)
         moptions['test'][2] = int(moptions['test'][2]) * (10**6)
      else: moptions['test'][1] = int(moptions['test'][1])/100.0
   else: moptions['test'] = ['N', '100']

   print("Predict")
   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['predict', '--help']);
      sys.exit(2)
  
 
   myMultiBiRNN.pred_entry(moptions)


def mTrain(margs):
   from scripts import myMultiBiRNN

   moptions, ErrorMessage = mCommonParam(margs)

   #moptions["wrkBase2"] = margs.wrkBase2
   #if moptions["wrkBase2"]==None or (not os.path.isdir(moptions["wrkBase2"])):
   #   ErrorMessage = ErrorMessage + ("\n\tThe directory (%s) does not exist" % moptions["wrkBase2"])

   moptions['fnum'] = margs.fnum
   ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
   moptions['hidden'] = margs.hidden
   ErrorMessage = ErrorMessage + non_negative(moptions['hidden'], 'hidden')

   if not margs.test==None:
      moptions['test'] = margs.test.split(',')
      if moptions['test'][0] == 'E': moptions['test'][0] = '-'
      elif moptions['test'][0] == 'P': moptions['test'][0] = '0'
      else:
         ErrorMessage = ErrorMessage + "Unknown option for test: the first character must be E or P "+margs.test
      if moptions['test'][0] in ['-']:
         moptions['test'][1] = int(moptions['test'][1]) * (10**6)
         moptions['test'][2] = int(moptions['test'][2]) * (10**6)
      else: moptions['test'][1] = int(moptions['test'][1])/100.0 
   else: moptions['test'] = ['N', '100']

   print("Train")
   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['train', '--help']);
      sys.exit(2)

   myMultiBiRNN.mMult_RNN_LSTM_train(moptions)

def mGetFeatures(margs):
   from scripts import myGetFeatureBasedPos

   moptions, ErrorMessage = mCommonParam(margs)
   moptions['posneg'] = margs.posneg
   moptions['fnum'] = margs.fnum
   ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
   moptions['size_per_batch'] = margs.size_per_batch
   if moptions['size_per_batch'] < 0.001: moptions['size_per_batch'] = 0.001

   moptions['basecall_1d'] = margs.basecall_1d
   moptions['basecall_2strand'] = margs.basecall_2strand

   moptions['region'] = [None, None, None]
   if not (margs.region==None or margs.region.strip()==''):
      rsp = margs.region.split(':')
      for rv_ind in range(len(rsp)):
         rsp[rv_ind] = rsp[rv_ind].strip();
         if not rsp[rv_ind]=='': 
            moptions['region'][rv_ind] = rsp[rv_ind]
   #moptions['region'] = [ ]
   #if margs.region == None or len(margs.region.strip())==0:
   #   moptions['region'].append([None, None, None])
   #else:
   #   mregionlist = margs.region.split(';')
   #   for mr in mregionlist:
   #      mr = mr.strip();
   #      if len(mr)==0: continue;
   #      mr_sp = mr.split(':')
   #      cur_r = [None, None, None]
   #      for rv_ind in range(len(mr_sp)):
   #         mr_sp[rv_ind] = mr_sp[rv_ind].strip();
   #         if not mr_sp[rv_ind]=='':
   #            cur_r[rv_ind] = mr_sp[rv_ind]
   #            if rv_ind>0: cur_r[rv_ind] = int(cur_r[rv_ind])
   #      moptions['region'].append(cur_r)

   moptions['Ref'] = margs.Ref
   if moptions['Ref']==None or (not os.path.isfile(moptions['Ref'])):
      ErrorMessage = ErrorMessage + ("\n\t reference file does not exist (%s)" % moptions['Ref'])
   moptions['motifORPos'] = margs.motifORPos
   if margs.motifORPos==1:
      moptions['motif'] = [margs.motif.upper(), margs.ModinMotif]
   elif margs.motifORPos==2:
      moptions['fulmod'] = margs.fulmod
      if moptions['fulmod']==None:
         ErrorMessage = ErrorMessage + ("\t There is no parameter for --fulmod.")
      moptions['anymod'] = margs.anymod
      if moptions['anymod'] == None:
         ErrorMessage = ErrorMessage + ("\t There is no parameter for --anymod.")
      moptions['nomod'] = margs.nomod
      if moptions['nomod'] == None:
         ErrorMessage = ErrorMessage + ("\t There is no parameter for --nomod.")
   else:
      ErrorMessage = ErrorMessage + ("\tmotifORPos value (%d) is not supported." % margs.motifORPos)

   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['getfeatures', '--help']);
      sys.exit(1)

   myGetFeatureBasedPos.getFeature_manager(moptions)


'''
def mGetFeaturesE(margs):
   #from scripts import myGetFeatureForRNN_OE
   from scripts import myGetFeatureBasedPos

   moptions, ErrorMessage = mCommonParam(margs)
   moptions['posneg'] = margs.posneg
   moptions['fnum'] = margs.fnum
   ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
   moptions['size_per_batch'] = margs.size_per_batch
   if moptions['size_per_batch'] < 0.001: moptions['size_per_batch'] = 0.001

   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['getfeaturesE', '--help']);
      sys.exit(1)


   #myGetFeatureForRNN_OE.getFeature_manager(moptions)
   myGetFeatureBasedPos.getFeature_manager(moptions)

def mGetNAfeatures(margs):
   from scripts import myNAGetFeature

   moptions, ErrorMessage = mCommonParam(margs)
   #moptions['featureFolder'] = moptions['outFolder']
   moptions['fulmod'] = margs.fulmod
   moptions['anymod'] = margs.anymod
   moptions['nomod'] = margs.nomod
   moptions['fnum'] = margs.fnum
   ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
   moptions['size_per_batch'] = margs.size_per_batch
   if moptions['size_per_batch'] < 0.001: moptions['size_per_batch'] = 0.001

   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['getNAfeatures', '--help']);
      sys.exit(1)


   myNAGetFeature.mGetFeature_manager(moptions)
'''

#####################################################################################

subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

com_group_for_comparison = parent_parser.add_argument_group('Common options.')
com_group_for_comparison.add_argument("--outLevel", type=int, choices=[OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR], default=OUTPUT_WARNING, help=("The level for output: %d for DEBUG, %d for INFO, %d for WARNING, %d for ERROR. Default: %d" % (OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR, OUTPUT_WARNING)))
com_group_for_comparison.add_argument("--wrkBase", help="The base folder for FAST5 files.")
com_group_for_comparison.add_argument("--FileID", default="mod", help="The unique string for intermediate files and final output files. Default: 'mod'")
com_group_for_comparison.add_argument("--outFolder", default='./mod_output', help="The default folder for outputing the results. Default: ./mod_output")
com_group_for_comparison.add_argument("--recursive", type=int, default=1, choices=[0,1], help="Recurise to find fast5 files. Default:1")
com_group_for_comparison.add_argument("--threads", type=int, default=4, help="The number of threads used (not for train). Default:4")
com_group_for_comparison.add_argument("--files_per_thread", type=int, default=1000, help="The number of fast5 files for each thread (not for train). Default:1000")
com_group_for_comparison.add_argument("--windowsize", type=int, default=51, help="The window size to extract features. Default: 51")
com_group_for_comparison.add_argument("--alignStr", type=str, default='bwa', choices=["bwa","minimap2"], help="Alignment tools (bwa or minimap2 is supported). Default: bwa")

parser_detect = subparsers.add_parser('detect', parents=[parent_parser], help="Detect modifications at a genomic scale", description="Detect modifications by integrating all long reads for a genome", epilog="For example, \n \
python %(prog)s --wrkBase ctrl_oligo_SpeI_cut --FileID mod_det --outFolder ./mod_output/detect3 \n \
", formatter_class=RawTextHelpFormatter)
parser_detect.add_argument("--Ref", help="The reference sequence")
parser_detect.add_argument("--predDet", type=int, default=1, choices=[0,1], help="pred first and then detect (1) or only detect (0). Default: 1")
parser_detect.add_argument("--predpath", default=None, help="The file path of predictions for each fast5 file. The file pattern is *_*.detail. Default: './mod_output/pred2/'")
parser_detect.add_argument("--modfile", type=str, default=None, help="The path to load training model. Default: 'mod_output/'")
parser_detect.add_argument("--fnum", type=int, default=53, help="The number of features. Default: 53")
parser_detect.add_argument("--hidden", type=int, default=100, help="The number of hidden node. Default: 100")
parser_detect.add_argument("--basecall_1d", default="Basecall_1D_000", help="Path for basecall_1d. Default: Basecall_1D_000")
#parser_detect.add_argument("--basecall_1d", default="Basecall_1D_000,Basecall_1D_001", help="Path for basecall_1d. Default: Basecall_1D_000")
parser_detect.add_argument("--basecall_2strand", default="BaseCalled_template", help="Path for basecall_2strand. Default: BaseCalled_template")
#parser_detect.add_argument("--region", default=":1000000:2000000", help="The region of interest: for example, chr:1:100000;chr2:10000");
parser_detect.add_argument("--region", default=None, help="The region of interest: for example, chr:1:100000;chr2:10000");
parser_detect.add_argument("--ConUnk", default=True, choices=[False, True], help="Whether contain unknown chromosome"); 
parser_detect.set_defaults(func=mDetect)


parser_predict = subparsers.add_parser('predict', parents=[parent_parser], help="Predict modifications for each long read in FAST5 files", description="Predict modifications for each long read in FAST5 files", epilog="For example, \n \
python %(prog)s --wrkBase ctrl_oligo_SpeI_cut --FileID mod_pred  --outFolder ./mod_output/pred2 \n \
", formatter_class=RawTextHelpFormatter)
parser_predict.add_argument("--modfile", type=str, default="./mod_output/train1/2/mod_train", help="The path to save mod")
parser_predict.add_argument("--fnum", type=int, default=53, help="The number of features. Default: 53")
parser_predict.add_argument("--hidden", type=int, default=100, help="The number of hidden node. Default: 100")
parser_predict.add_argument("--basecall_1d", default="Basecall_1D_000", help="Path for basecall_1d. Default: Basecall_1D_000")
parser_predict.add_argument("--basecall_2strand", default="BaseCalled_template", help="Path for basecall_2strand. Default: BaseCalled_template")
parser_predict.add_argument("--test", help="The number of E Coli genomic position for testing.  Default: 'I,1,2'")
parser_predict.set_defaults(func=mPredict)


parser_training = subparsers.add_parser('train', parents=[parent_parser], help="Training a modification classifier", description="Training a modification classifier", epilog="For example, \n \
python %(prog)s --wrkBase /scr1/users/liuq1/project/deepnanomod/aoe53features/umr --wrkBase2 /scr1/users/liuq1/project/deepnanomod/aoe53features/sss --FileID mod_train --outFolder ./mod_output/train1 \n \
", formatter_class=RawTextHelpFormatter)
parser_training.add_argument("--wrkBase2", help="The base folder for long reads without any modifications.")
parser_training.add_argument("--fnum", type=int, default=53, help="The number of features. Default: 53")
parser_training.add_argument("--hidden", type=int, default=100, help="The number of hidden node. Default: 100")
parser_training.add_argument("--test", help="The number of E Coli genomic position for testing. Default: 'E,1,2'")
parser_training.set_defaults(func=mTrain)

parser_getfeatures = subparsers.add_parser('getfeatures', parents=[parent_parser], help="Get features for all fast5 files", description="Get features for all fast5 files", epilog="For example, \n \
python %(prog)s --wrkBase /mnt/isilon/wang_lab/liuq1/nanopore/nanopolish/umr/160617_ecolilowinput_UMR9/called/pass --threads 48 --recursive 0 --posneg 0 --outFolder /scr1/users/liuq1/project/deepnanomod/aoe53features/umr  \n \
python %(prog)s --wrkBase /mnt/isilon/wang_lab/liuq1/nanopore/nanopolish/sss/160617_ecolilowinput_sssiR9/called/pass --threads 48 --recursive 0 --posneg 1 --outFolder /scr1/users/liuq1/project/deepnanomod/aoe53features/sss \n \
", formatter_class=RawTextHelpFormatter)
parser_getfeatures.add_argument("--posneg", type=int, default=0, choices=[0,1], help="The positive(1) or negative(0) class. Default: 0")
parser_getfeatures.add_argument("--size_per_batch", type=int, default=7, help="The size (unit: 10^7=10M) of a feature file. Default: 7")
parser_getfeatures.add_argument("--fnum", type=int, default=53, help="The number of features. Default: 53")
parser_getfeatures.add_argument("--region", type=str, help="The region of interest. Set to None or empty for all. Format is chr:start_pos:end_pos")
parser_getfeatures.add_argument("--basecall_1d", default="Basecall_1D_000", help="Path for basecall_1d. Default: Basecall_1D_000")
parser_getfeatures.add_argument("--basecall_2strand", default="BaseCalled_template", help="Path for basecall_2strand. Default: BaseCalled_template")

parser_getfeatures.add_argument("--motifORPos", type=int, default=1, help="Use Motif (1) or pos (2) for modified bases. Default: 1")

parser_getfeatures.add_argument("--motif", default='CG', type=str, help="The motif of interest")
parser_getfeatures.add_argument("--ModinMotif", default=0, type=int, help="The motif of interest")
parser_getfeatures.add_argument("--Ref", help="The reference sequence")

parser_getfeatures.add_argument("--fulmod", type=str, help="The file pattern for full modification: bisultfiteseq/chr20_C*_0.95.txt")
parser_getfeatures.add_argument("--anymod", type=str, help="The file pattern for any modification: bisultfiteseq/chr20_any_0.95.txt")
parser_getfeatures.add_argument("--nomod", type=str, help="The file pattern for any modification: bisultfiteseq/chr20_no1_0.95.txt")

parser_getfeatures.set_defaults(func=mGetFeatures)


'''
parser_getfeaturesE = subparsers.add_parser('getfeaturesE', parents=[parent_parser], help="Get features for all fast5 files", description="Get features for all fast5 files", epilog="For example, \n \
python %(prog)s --wrkBase /mnt/isilon/wang_lab/liuq1/nanopore/nanopolish/umr/160617_ecolilowinput_UMR9/called/pass --threads 48 --recursive 0 --posneg 0 --outFolder /scr1/users/liuq1/project/deepnanomod/aoe53features/umr  \n \
python %(prog)s --wrkBase /mnt/isilon/wang_lab/liuq1/nanopore/nanopolish/sss/160617_ecolilowinput_sssiR9/called/pass --threads 48 --recursive 0 --posneg 1 --outFolder /scr1/users/liuq1/project/deepnanomod/aoe53features/sss \n \
", formatter_class=RawTextHelpFormatter)
parser_getfeaturesE.add_argument("--posneg", type=int, default=0, choices=[0,1], help="The positive(1) or negative(0) class. Default: 0")
parser_getfeaturesE.add_argument("--size_per_batch", type=float, default=7, help="The size (unit: 10^8=100M) of a feature file. Default: 7")
# for uni2feature
#parser_getfeaturesE.add_argument("--size_per_batch", type=int, default=3, help="The size (unit: 10^7=10M) of a feature file. Default: 7")
parser_getfeaturesE.add_argument("--fnum", type=int, default=53, help="The number of features. Default: 53")
parser_getfeaturesE.set_defaults(func=mGetFeaturesE)


parser_getNAfeatures= subparsers.add_parser('getNAfeatures', parents=[parent_parser], help="Get features for all fast5 files", description="Get features for all fast5 files", epilog="For example, \n \
python %(prog)s --wrkBase /scr1/users/liuq1/project/nanopore/nanodeepmod/na12878albacore126/chr7 --threads 30 --recursive 1 --outFolder data_NanoDeepMod/na12878/chr7 --fulmod \"bisultfiteseq/chr7_C*_0.95.txt\" --anymod \"bisultfiteseq/chr7_any_0.95.txt\"  \n \
python %(prog)s --wrkBase /scr1/users/liuq1/project/nanopore/nanodeepmod/na12878albacore126/chr20 --threads 30 --recursive 1 --outFolder data_NanoDeepMod/na12878/chr20 --fulmod \"bisultfiteseq/chr20_C*_0.95.txt\" --anymod \"bisultfiteseq/chr20_any_0.95.txt\"  \n \
", formatter_class=RawTextHelpFormatter)
parser_getNAfeatures.add_argument("--size_per_batch", type=float, default=7, help="The size (unit: 10^8=100M) of a feature file. Default: 7")
parser_getNAfeatures.add_argument("--fulmod", default='bisultfiteseq/chr20_C*_0.95.txt', help="The file pattern for full modification")
parser_getNAfeatures.add_argument("--anymod", default="bisultfiteseq/chr20_any_0.95.txt", help="The file pattern for any modification")
parser_getNAfeatures.add_argument("--nomod", default="bisultfiteseq/chr20_no1_0.95.txt", help="The file pattern for no modification")
parser_getNAfeatures.add_argument("--fnum", type=int, default=53, help="The number of features. Default: 53")
parser_getNAfeatures.set_defaults(func=mGetNAfeatures)
'''


if len(sys.argv)<2:
   parser.print_help();
else:
   args = parser.parse_args()
   args.func(args);

