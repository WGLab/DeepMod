from collections import defaultdict
import csv, time, itertools, copy, h5py, time, re, random

import datetime, os, shutil, argparse, sys, pysam

import multiprocessing as mp
import numpy as np

import numpy.lib.recfunctions as rf

from pathlib import Path





def run(params):
    from src import modDetect
    read_pred_file=modDetect.per_read_predict(params)
    
    site_pred_file=modDetect.per_site_detect(read_pred_file, params)
    
    
    

if __name__ == '__main__':

    t=time.time()

    print('%s: Starting DeepMod.' %str(datetime.datetime.now()))
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--bam", help='Path to bam file if Guppy basecaller is user. BAM file is not needed with Tombo fast5 files.', type=str)
    
    parser.add_argument("-chrom",  "--chrom", nargs='*',  help='A space/whitespace separated list of contigs, e.g. chr3 chr6 chr22. Only applicable with guppy basecaller.')
    
    parser.add_argument("--ref", help='Path to reference file', type=str)
    
    parser.add_argument("--fast5", help='Path to folder containing tombo requiggle Fast5 files. Fast5 files will be recusrviely searched', type=str, required=True)
    
    parser.add_argument("--output", help='Path to folder where features will be stored', type=str, required=True)
    
    parser.add_argument("--threads", help='Number of processors to use',type=int, default=1)
    
    parser.add_argument("--file_name", help='Name of the output file',type=str, default='output')
    
    parser.add_argument("--model", help='Name of the model. Current options are "guppy_na12878" and "tombo_na12878". Use according to the basecaller specified.',type=str, default='guppy_na12878')
        
    parser.add_argument("--guppy_group", help='Name of the guppy basecall group',type=str, default='Basecall_1D_000')
    
    parser.add_argument("--tombo_group", help='Name of the tombo group',type=str, default='RawGenomeCorrected_000')

    parser.add_argument("--basecaller", help='Use Tombo or Guppy output',type=str, default='guppy', choices=['guppy', 'tombo'])
    
    
    parser.add_argument('-wgs_contigs_type','--wgs_contigs_type', \
                        help="""Options are "with_chr", "without_chr" and "all",\
                        "with_chr" option will assume \
                        human genome and run DeepMod on chr1-22 X Y, "without_chr" will \
                        run on chromosomes 1-22 X Y if the BAM and reference genome files \
                        use chromosome names without "chr". "all" option will run \
                        DeepMod on each contig present in reference genome FASTA file. Only applicable with guppy basecaller.""", \
                        type=str, default='all') 
    
    args = parser.parse_args()
    
    if args.basecaller=='guppy' and args.bam==None:
        print('BAM file needed with Guppy basecalled fast5 files')
        sys.exit(1)
    if args.chrom:
        chrom_list= args.chrom
        
    else:
        if args.wgs_contigs_type=='with_chr':
            chrom_list=['chr%d' %d for d in range(1,23)] + ['chrX','chrY']

        elif args.wgs_contigs_type == 'without_chr':
            chrom_list=['%d' %d for d in range(1,23)] + ['X', 'Y']

        elif args.wgs_contigs_type == 'all':
            fastafile=pysam.FastaFile(args.ref)
            chrom_list=fastafile.references
                
    params={'fast5':args.fast5, 'bam_path':args.bam, 'fasta_path':args.ref ,'output':args.output, 'threads':args.threads, 'file_name':args.file_name, 'window':10, 'chrom_list':chrom_list, 'model':args.model, 'guppy_group':args.guppy_group, 'tombo_group': args.tombo_group, 'basecaller':args.basecaller}
    
    os.makedirs(params['output'], exist_ok=True)
    
    print('\n%s: \nCommand: python %s\n' %(str(datetime.datetime.now()), ' '.join(sys.argv)))
    
    with open(os.path.join(args.output,'args'),'w') as file:
        file.write('Command: python %s\n\n\n' %(' '.join(sys.argv)))
        file.write('------Parameters Used For Running DeepMod------\n')
        for k in vars(args):
            file.write('{}: {}\n'.format(k,vars(args)[k]) )
                
    run(params)
    
    print('\n%s: Time elapsed=%.4fs' %(str(datetime.datetime.now()),time.time()-t))
