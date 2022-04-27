from collections import defaultdict, ChainMap

import csv, time, itertools, copy, h5py, time, pysam

import datetime, os, shutil, argparse, sys, random

import multiprocessing as mp
import numpy as np

import numpy.lib.recfunctions as rf

from pathlib import Path

from .utils import *

def get_bam_info(args):
    chrom, bam_path, fasta_path=args
    bam=pysam.Samfile(bam_path,'rb')
    read_info={}
    
    fastafile=pysam.FastaFile(fasta_path)
    
    ref_seq=fastafile.fetch(chrom)

    for pcol in bam.pileup(contig=chrom, flag_filter=0x4|0x100|0x200|0x400|0x800):
        try:
            if ref_seq[pcol.pos]=='C' and ref_seq[pcol.pos+1]=='G':
                for read in pcol.pileups:
                    if read.alignment.is_reverse==False:
                        if not read.is_del:
                            if read.alignment.qname not in read_info:
                                read_info[read.alignment.qname]=[(0, chrom)]
                            read_info[read.alignment.qname].append((pcol.pos+1, read.query_position))
        
            elif ref_seq[pcol.pos]=='G' and ref_seq[pcol.pos-1]=='C':
                for read in pcol.pileups:
                    if read.alignment.is_reverse:
                        if not read.is_del:
                            if read.alignment.qname not in read_info:
                                read_info[read.alignment.qname]=[(1, chrom)]
                            read_info[read.alignment.qname].append((pcol.pos+1, read.alignment.query_length-read.query_position-1))
                            
        except IndexError:
            continue
                            
    return read_info


def process_bam(params, pool):
    fastafile=pysam.FastaFile(params['fasta_path'])
    
    print('%s: Processing BAM File.' %str(datetime.datetime.now()))
    
    bam_info=[x for x in pool.imap_unordered(get_bam_info, zip(params['chrom_list'], itertools.repeat(params['bam_path']), itertools.repeat(params['fasta_path'])))]
    
    read_info=ChainMap(*bam_info)
    
    print('%s: BAM File Processed.' %str(datetime.datetime.now()))
    
    return read_info




def get_dataset_keys(f):
    keys = []
    f.visit(lambda key : keys.append(key) if (isinstance(f[key], h5py.Dataset) and key[:4]=='Raw/' and key[-7:]=='/Signal') else None)
    return keys[0]

def get_read_signal(f5, guppy_group):
    segment=segment=f5['/Analyses/%s/' %guppy_group].attrs['segmentation'].decode('utf-8')
    
    start=int(get_attr(f5['/Analyses/%s/Summary/' %segment],'first_sample_template'))
    stride=int(get_attr(f5['/Analyses/%s/Summary/' %guppy_group],'block_stride'))
    
    signal_key=get_dataset_keys(f5)
    signal=np.array(f5[signal_key])
    
    median=np.median(signal)
    mad=np.median(np.abs(signal-median))
    
    norm_signal=(signal-median)/mad
    
    base_level_data=[]
    base=0
    
    prev=start
    for i,x in enumerate(f5['/Analyses/%s/BaseCalled_template/Move' %guppy_group][()]):
        if x:        
            base_level_data.append([prev,(i+1)*stride+start-prev])
            prev=(i+1)*stride+start
    
    return norm_signal, base_level_data

def getFeatures(f5_chunk, params, read_info):
    base_map={'A':0, 'C':1, 'G':2, 'T':3}
    
    window=params['window']

    features_list=[]
    pos_list=[]
    chr_list=[]
    strand_list=[]
    read_names_list=[]
    
    for filename in f5_chunk:
        f5 = h5py.File(str(filename), 'r')

        try:
            read_name=get_attr(f5['Raw/Reads'], 'read_id').decode('utf-8')
            read_pos_list=read_info[read_name]

        except KeyError:
            continue

        norm_signal, base_level_data = get_read_signal(f5, params['guppy_group'])
        seq_len=len(base_level_data)

        fq=f5['/Analyses/%s/BaseCalled_template/Fastq' %params['guppy_group']][()].decode('utf-8').split('\n')[1]
    
        mapped_strand, mapped_chrom=read_pos_list[0]
        
        for x in read_pos_list[1:]:
            pos, read_pos=x

            if read_pos>window and read_pos<seq_len-window-1:
                mat=[]
                base_seq=[]
                
                for x in range(read_pos-window, read_pos+window+1):
                    mat.append([np.mean(norm_signal[base_level_data[x][0]:base_level_data[x][0]+base_level_data[x][1]]), np.std(norm_signal[base_level_data[x][0]:base_level_data[x][0]+base_level_data[x][1]]), base_level_data[x][1]])
                    base_seq.append(base_map[fq[x]])
                    
                base_seq=np.eye(4)[base_seq]
                mat=np.hstack((np.array(mat), base_seq))

                if np.size(mat)!=21*7:  
                    continue
            
                features_list.append(mat)
                pos_list.append(pos)
                chr_list.append(mapped_chrom)
                strand_list.append(mapped_strand)
                read_names_list.append(read_name)     
            
            
    features_list=np.array(features_list)
    
    return pos_list, chr_list, strand_list, read_names_list, features_list