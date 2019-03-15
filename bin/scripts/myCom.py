
#
# output level
#
OUTPUT_DEBUG=0;
OUTPUT_INFO=1;
OUTPUT_WARNING=2;
OUTPUT_ERROR=3;

#
# base-pairing
# 4 base types with case insensitive
#
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

g_ACGT = ['A', 'C', 'G', 'T']

acgt = na_bp.keys();

#
# get complementary base of a given base
#
def getComplementary(na):
   com_na = []
   for i in range(len(na)):
      com_na.append(na_bp[na[i]])
   ''.join(com_na[::-1])

#
# format the end of a folder
#
def format_last_letter_of_folder(cursub):
   if not cursub==None:
      if cursub[-1]=='/': return cursub;
      elif cursub[-1]=='\\': return cursub[:-1]+'/';
      else: return cursub+'/';

#
# default path of data in fast5 files
#
analyses_base = "Analyses"
basecall_events_base = "Events"
raw_base = 'Raw'
reads_base = "Reads"
signal_base = "Signal"
basecall_fastq_base = "Fastq"




