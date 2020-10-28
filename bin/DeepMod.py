#!/usr/bin/env python

import os;
import sys;

import string;

from collections import defaultdict

import argparse;
from argparse import RawTextHelpFormatter

from DeepMod_scripts.myCom import *



# three modules in DeepMod
parser = argparse.ArgumentParser(description="Detect nucleotide modification from nanopore signals data.", epilog="For example, \n \
\tpython %(prog)s train: Training a modification classifier.\n \
\tpython %(prog)s detect: Detect modification by integrating all long reads. \n \
\tpython %(prog)s getfeatures: Get features for training a model.  \n \
", formatter_class=RawTextHelpFormatter);


#
# Return error message when a value<1
# Return an empty string otherwise
#
def non_negative(i, mstr):
   if i<1: return (("\n\tError %d could not be negative(%d)" % (mstr, i)))
   else: return ''

#
# Print all parameters in stdout
#
def printParameters(moptions):
   mpkeys = moptions.keys(); #mpkeys.sort()
   sorted(mpkeys)
   print('%30s: %s' % ('Current directory', os.getcwd()))
   for mpk in mpkeys:
      print ('%30s: %s' % (mpk, str(moptions[mpk])))
   sys.stdout.flush()

#
# Got common argument provided by users or default values.
#
#
def mCommonParam(margs):

   ErrorMessage = ""
   moptions = defaultdict()
   # how to output running message: need more control now.
   moptions['outLevel'] = margs.outLevel
   # the input working base
   moptions["wrkBase"] = margs.wrkBase
   if moptions["wrkBase"]==None:
      ErrorMessage = ErrorMessage + ("\n\tThe input folder is None.")

   # An unique ID for output
   # Usefull for run the program in parallel
   moptions["FileID"] = margs.FileID
   # output folder;
   # make it if the output folder does not exist
   moptions['outFolder'] = margs.outFolder
   moptions['outFolder'] = format_last_letter_of_folder(moptions['outFolder'])
   if moptions['outFolder']==None or (not os.path.isdir(moptions['outFolder'])):
      try:
         os.system('mkdir -p '+moptions['outFolder']);
      except:
         ErrorMessage = ErrorMessage + ("\n\tThe output folder (%s) does not exist and cannot be created." % moptions['outFolder'])

   # check all data in a recurive way
   moptions['recursive'] = margs.recursive
   # the number of threads used and the number of files handled by each thread.
   moptions['files_per_thread'] = margs.files_per_thread
   if moptions['files_per_thread']<2: moptions['files_per_thread'] = 2
   # the number of threads used
   moptions['threads'] = margs.threads
   if moptions['threads']<1: moptions['threads'] = 1

   # windowsize: default=21
   moptions['windowsize'] = margs.windowsize
   ErrorMessage = ErrorMessage + non_negative(moptions['windowsize'], 'windowsize')
   if moptions['windowsize']<1: moptions['windowsize'] = 1

   # aligners: bwa-mem or minimap2
   moptions['alignStr'] = margs.alignStr;

   moptions['SignalGroup'] = margs.SignalGroup;

   moptions['move'] = margs.move

   return [moptions, ErrorMessage]

#
# detect modification for bases of interests
# input is a list of fast5 files, a reference genome and a well-trained model.
#
def mDetect(margs):
   # get common parameters
   moptions, ErrorMessage = mCommonParam(margs)

   # path for basecall information in fast5 files
   moptions['basecall_1d'] = margs.basecall_1d
   moptions['basecall_2strand'] = margs.basecall_2strand
   # Whether consider those chromosome which contain -_:/
   # default: yes;
   moptions['ConUnk'] = margs.ConUnk
   # output layer information for deep learning
   moptions['outputlayer'] = margs.outputlayer
   # base of interest
   moptions['Base'] = margs.Base
   # whether take cluster effect of methylation into consideration
   moptions['mod_cluster'] = margs.mod_cluster
   # base of interest
   if moptions['Base'] in ["", None]:
      ErrorMessage = ErrorMessage + ("\n\t Please provide a base of interest.")

   # predict medification for bases of interest in long reads first
   # only summarize them for each genomic position of interest .
   moptions['predDet'] = margs.predDet
   if moptions['predDet']:
      # path to reference genome
      moptions['Ref'] = margs.Ref
      if moptions['Ref']==None or (not os.path.isfile(moptions['Ref'])):
         ErrorMessage = ErrorMessage + ("\n\t reference file does not exist (%s)" % moptions['Ref'])

      # the number of feature for each event
      moptions['fnum'] = margs.fnum
      ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
      # the number of hidden nodes
      moptions['hidden'] = margs.hidden
      ErrorMessage = ErrorMessage + non_negative(moptions['hidden'], 'hidden')
      # the well-trained model
      moptions['modfile'] = margs.modfile
      if moptions['modfile']==None:
         print("No mod file is provided. The default one is used")
         moptions['modfile'] = ('train_deepmod/rnn_P90wd%d_f53/mod_train_P90wd%d_f53' % (moptions['windowsize'], moptions['windowsize']))
         if (not os.path.isfile(moptions['modfile']+'.meta')):
            moptions['modfile'] = ('{}/lib/python{}.{}/site-packages/DeepMod/train_deepmod/rnn_P90wd{}_f53/mod_train_P90wd{}_f53'.format(sys.prefix,sys.version_info.major,sys.version_info.minor, moptions['windowsize'], moptions['windowsize']))
      if (not os.path.isfile(moptions['modfile']+'.meta')):
         ErrorMessage = ErrorMessage + ("\n\tThe meta file (%s) does not exist" % (moptions['modfile']+'.meta' if not moptions['modfile']==None else ""))
   else:
      # already done the prediction process?
      # Yes: summarize the results only
      moptions['predpath'] = margs.predpath
      if moptions['predpath']==None or (not os.path.isdir(moptions['predpath'])):
         ErrorMessage = ErrorMessage + ("\n\tThe predpath does not exist")

   # specify region of interest
   # not consider bases outside regions in a reference genome
   # None: all bases of interest
   moptions['region'] = [ ]
   if margs.region == None or len(margs.region)==0:
      moptions['region'].append([None, None, None])
   else:
      mregionlist = margs.region.split(';')
      for mr in mregionlist:
         mr_sp = mr.split(':')
         moptions['region'].append([mr_sp[0], int(mr_sp[1]) if len(mr_sp)>1 else None, int(mr_sp[2]) if len(mr_sp)>2 else None ])

   print("\nNanopore sequencing data analysis is resourece-intensive and time consuming. ")
   print("Some potential strong recommendations are below:")
   print("\tIf your reference genome is large as human genome and your Nanopore data is huge,")
   print("\tIt would be faster to run this program parallelly to speed up.")
   print("\tYou might run different input folders of your fast5 files and ")
   print("\tgive different output names (--FileID) or folders (--outFolder)")
   print("\tA good way for this is to run different chromosome individually.\n")

   # print help information if any necessary options are not provided.
   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['detect', '--help']);
      sys.exit(1)

   from DeepMod_scripts import myDetect
   myDetect.mDetect_manager(moptions)

#
# Train a model
# Need to get features first.
#
def mTrain(margs):
   from DeepMod_scripts import myMultiBiRNN

   # gent common options
   moptions, ErrorMessage = mCommonParam(margs)

   # network setting: the number of features and the number of hidden nodes
   moptions['fnum'] = margs.fnum
   ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
   moptions['hidden'] = margs.hidden
   ErrorMessage = ErrorMessage + non_negative(moptions['hidden'], 'hidden')

   # the output function of the deep learning model
   moptions['outputlayer'] = margs.outputlayer
   # whether using different class weights
   moptions['unbalanced'] = margs.unbalanced

   # re-load trained model and continue to train
   moptions['modfile'] = margs.modfile
   if moptions['modfile']==None: pass;
   elif (not os.path.isfile(moptions['modfile']+'.meta')):
      ErrorMessage = ErrorMessage + ("\n\tThe meta file (%s) does not exist" % (moptions['modfile']+'.meta' if not moptions['modfile']==None else ""))

   # read-based or region based independent training
   # E: region-based
   # P: read-based.
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

   # print help document if necessary options are not provided.
   print("Train")
   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['train', '--help']);
      sys.exit(2)

   myMultiBiRNN.mMult_RNN_LSTM_train(moptions)

#
# get features for training
#
#
def mGetFeatures(margs):
   from DeepMod_scripts import myGetFeatureBasedPos

   # get common options
   moptions, ErrorMessage = mCommonParam(margs)
   # motif-based data: positive or negative control data
   moptions['posneg'] = margs.posneg
   # the number of features: 7-description or 57-description
   moptions['fnum'] = margs.fnum
   ErrorMessage = ErrorMessage + non_negative(moptions['fnum'], 'fnum')
   # size of each bacth to store features
   moptions['size_per_batch'] = margs.size_per_batch
   if moptions['size_per_batch'] < 0.001: moptions['size_per_batch'] = 0.001

   # path to basecall inform in fast5 files
   moptions['basecall_1d'] = margs.basecall_1d
   moptions['basecall_2strand'] = margs.basecall_2strand

   # regions of interest
   moptions['region'] = [None, None, None]
   if not (margs.region==None or margs.region.strip()==''):
      rsp = margs.region.split(':')
      for rv_ind in range(len(rsp)):
         rsp[rv_ind] = rsp[rv_ind].strip();
         if not rsp[rv_ind]=='':
            moptions['region'][rv_ind] = rsp[rv_ind]

   # referene genome
   moptions['Ref'] = margs.Ref
   if moptions['Ref']==None or (not os.path.isfile(moptions['Ref'])):
      ErrorMessage = ErrorMessage + ("\n\t reference file does not exist (%s)" % moptions['Ref'])

   # get motif-based modification
   # or specify by --fulmod/--anymod/--nomod
   moptions['motifORPos'] = margs.motifORPos
   if margs.motifORPos==1:
      moptions['motif'] = [margs.motif.upper(), margs.ModinMotif]
   elif margs.motifORPos==2:
      moptions['fulmod'] = margs.fulmod
      if moptions['fulmod']==None: # completely modificated positions
         ErrorMessage = ErrorMessage + ("\t There is no parameter for --fulmod.")
      moptions['anymod'] = margs.anymod
      if moptions['anymod'] == None: # patially modificated positions
         ErrorMessage = ErrorMessage + ("\t There is no parameter for --anymod.")
      moptions['nomod'] = margs.nomod
      if moptions['nomod'] == None: # completely unmodified posisionts
         ErrorMessage = ErrorMessage + ("\t There is no parameter for --nomod.")
   else:
      ErrorMessage = ErrorMessage + ("\tmotifORPos value (%d) is not supported." % margs.motifORPos)

   # print help document if any required options are not provided.
   printParameters(moptions)
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print(ErrorMessage)
      parser.print_help();
      parser.parse_args(['getfeatures', '--help']);
      sys.exit(1)

   myGetFeatureBasedPos.getFeature_manager(moptions)


#####################################################################################

subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

# add common options
com_group_for_comparison = parent_parser.add_argument_group('Common options.')
com_group_for_comparison.add_argument("--outLevel", type=int, choices=[OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR], default=OUTPUT_WARNING, help=("The level for output: %d for DEBUG, %d for INFO, %d for WARNING, %d for ERROR. Default: %d" % (OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR, OUTPUT_WARNING)))
com_group_for_comparison.add_argument("--wrkBase", help="The base folder for FAST5 files.")
com_group_for_comparison.add_argument("--FileID", default="mod", help="The unique string for intermediate files and final output files. Default: 'mod'")
com_group_for_comparison.add_argument("--outFolder", default='./mod_output', help="The default folder for outputing the results. Default: ./mod_output")
com_group_for_comparison.add_argument("--recursive", type=int, default=1, choices=[0,1], help="Recurise to find fast5 files. Default:1")
com_group_for_comparison.add_argument("--threads", type=int, default=4, help="The number of threads used (not for train). Default:4")
com_group_for_comparison.add_argument("--files_per_thread", type=int, default=1000, help="The number of fast5 files for each thread (not for train). Default:500")
com_group_for_comparison.add_argument("--windowsize", type=int, default=21, help="The window size to extract features. Default: 21")
com_group_for_comparison.add_argument("--alignStr", type=str, default='minimap2', choices=["bwa","minimap2"], help="Alignment tools (bwa or minimap2 is supported). Default: minimap2")
com_group_for_comparison.add_argument("--SignalGroup", type=str, default='simple', choices=["simple","rundif"], help="How to associate signals to each called bases. Default: simple")
com_group_for_comparison.add_argument("--move", default=False, action="store_true", help="Whether the basecalled data use move tables instead of event tables. Default: False")

# add detection options
parser_detect = subparsers.add_parser('detect', parents=[parent_parser], help="Detect modifications at a genomic scale", description="Detect modifications by integrating all long reads for a genome", epilog="For example, \n \
python %(prog)s --wrkBase ctrl_oligo_SpeI_cut --FileID mod_det --outFolder ./mod_output/detect3 \n \
", formatter_class=RawTextHelpFormatter)
parser_detect.add_argument("--Ref", help="The reference sequence")
parser_detect.add_argument("--predDet", type=int, default=1, choices=[0,1], help="pred first and then detect (1) or only detect (0). Default: 1")
parser_detect.add_argument("--predpath", default=None, help="The file path of predictions for each fast5 file. The file pattern is *_*.detail. Default: './mod_output/pred2/'")
parser_detect.add_argument("--modfile", type=str, default=None, help="The path to load training model. Default: 'mod_output/'")
parser_detect.add_argument("--fnum", type=int, default=7, help="The number of features. Default: 7")
parser_detect.add_argument("--hidden", type=int, default=100, help="The number of hidden node. Default: 100")
parser_detect.add_argument("--basecall_1d", default="Basecall_1D_000", help="Path for basecall_1d. Default: Basecall_1D_000")
parser_detect.add_argument("--basecall_2strand", default="BaseCalled_template", help="Path for basecall_2strand. Default: BaseCalled_template")
parser_detect.add_argument("--region", default=None, help="The region of interest: for example, chr:1:100000;chr2:10000");
parser_detect.add_argument("--ConUnk", default=True, choices=[False, True], help="Whether contain unknown chromosome");
parser_detect.add_argument("--outputlayer", default="", choices=["", "sigmoid"], help="how to put activation function for output layer")
parser_detect.add_argument("--Base", type=str, default='C', choices=['A', 'C', 'G', 'T'], help="Interest of bases");
parser_detect.add_argument("--mod_cluster", default=0, choices=[0,1], help="1: CpG cluster effect; 0: not");
parser_detect.set_defaults(func=mDetect)

# add training options
parser_training = subparsers.add_parser('train', parents=[parent_parser], help="Training a modification classifier", description="Training a modification classifier", epilog="For example, \n \
python %(prog)s --wrkBase umr --wrkBase2 sss --FileID mod_train --outFolder ./mod_output/train1 \n \
", formatter_class=RawTextHelpFormatter)
parser_training.add_argument("--wrkBase2", help="The base folder for long reads without any modifications.")
parser_training.add_argument("--fnum", type=int, default=7, help="The number of features. Default: 7")
parser_training.add_argument("--hidden", type=int, default=100, help="The number of hidden node. Default: 100")
parser_training.add_argument("--modfile", type=str, default=None, help="The path to load training model. Default: 'mod_output/'")
parser_training.add_argument("--test", help="The number of E Coli genomic position for testing. Default: 'E,1,2'")
parser_training.add_argument("--outputlayer", default="", choices=["", "sigmoid"], help="how to put activation function for output layer")
parser_training.add_argument("--unbalanced", type=int, default=0, choices=[1, 0, None], help="Whether data is unbalanced");
parser_training.set_defaults(func=mTrain)

# add get-feature options
parser_getfeatures = subparsers.add_parser('getfeatures', parents=[parent_parser], help="Get features for all fast5 files", description="Get features for all fast5 files", epilog="For example, \n \
python %(prog)s --wrkBase umr/160617_ecolilowinput_UMR9/called/pass --threads 48 --recursive 0 --posneg 0 --outFolder umr  \n \
python %(prog)s --wrkBase sss/160617_ecolilowinput_sssiR9/called/pass --threads 48 --recursive 0 --posneg 1 --outFolder sss \n \
", formatter_class=RawTextHelpFormatter)
parser_getfeatures.add_argument("--posneg", type=int, default=0, choices=[0,1], help="The positive(1) or negative(0) class. Default: 0")
parser_getfeatures.add_argument("--size_per_batch", type=int, default=1, help="The size (unit: 10^7=10M) of a feature file. Default: 1")
parser_getfeatures.add_argument("--fnum", type=int, default=7, help="The number of features. Default: 7")
parser_getfeatures.add_argument("--region", type=str, help="The region of interest. Set to None or empty for all. Format is chr:start_pos:end_pos")
parser_getfeatures.add_argument("--basecall_1d", default="Basecall_1D_000", help="Path for basecall_1d. Default: Basecall_1D_000")
parser_getfeatures.add_argument("--basecall_2strand", default="BaseCalled_template", help="Path for basecall_2strand. Default: BaseCalled_template")

parser_getfeatures.add_argument("--motifORPos", type=int, default=1, help="Use Motif (1) or pos (2) for modified bases. Default: 1")

parser_getfeatures.add_argument("--motif", default='CG', type=str, help="The motif of interest")
parser_getfeatures.add_argument("--ModinMotif", default=0, type=int, help="The position of modified base in the motif of interest")
parser_getfeatures.add_argument("--Ref", help="The reference sequence")

parser_getfeatures.add_argument("--fulmod", type=str, help="The file pattern for full modification: bisultfiteseq/chr20_C*_0.95.txt")
parser_getfeatures.add_argument("--anymod", type=str, help="The file pattern for any modification: bisultfiteseq/chr20_any_0.95.txt")
parser_getfeatures.add_argument("--nomod", type=str, help="The file pattern for any modification: bisultfiteseq/chr20_no1_0.95.txt")

parser_getfeatures.set_defaults(func=mGetFeatures)

# no provided argument
# print help document
if len(sys.argv)<2:
   parser.print_help();
else:
   args = parser.parse_args()
   args.func(args);
