from collections import defaultdict, ChainMap

import csv, time, itertools, copy, h5py, time, re

import datetime, os, shutil, argparse, sys, random

import multiprocessing as mp
import numpy as np

import numpy.lib.recfunctions as rf

from pathlib import Path
from .utils import *

def get_tombo_events_summary(events):
    rf_events=rf.structured_to_unstructured(events['norm_mean','norm_stdev','length'])
    
    rf_signal_cordinates=rf.structured_to_unstructured(events['length', 'start'])
    
    seq=''.join([x.decode('utf-8') for x in events['base']])
    
    return seq, rf_events, rf_signal_cordinates

def get_tombo_alignment_info(alignment_attrs):
    mapped_chrom = alignment_attrs['mapped_chrom']
    mapped_strand = alignment_attrs['mapped_strand']
    mapped_start = int(alignment_attrs['mapped_start'])
    mapped_end = int(alignment_attrs['mapped_end'])
    
    return mapped_chrom, mapped_strand, mapped_start, mapped_end


def getFeatures(f5_list, params):
    base_map={'A':0, 'C':1, 'G':2, 'T':3}
    strand_map={'+':0, '-':1}
    
    tombo_group=params['tombo_group']

    window=params['window']

    features_list=[]
    pos_list=[]
    chr_list=[]
    label_list=[]
    strand_list=[]
    read_names_list=[]
    
    matcher=re.compile('CG')


    for filename in f5_list:
        f5 = h5py.File(str(filename), 'r')
        
        #get tombo resquiggle data
        try:
            #get event information
            
            read_name=get_attr(f5['Raw/Reads'], 'read_id').decode('utf-8')
            
            events = f5['Analyses/%s/BaseCalled_template/Events' %tombo_group]
            seq, rf_events, rf_signal_cordinates=get_tombo_events_summary(events)
            
            #get alignment attributes
            alignment_attrs = f5['Analyses/%s/BaseCalled_template/Alignment' %tombo_group].attrs
            mapped_chrom, mapped_strand, mapped_start, mapped_end = get_tombo_alignment_info(alignment_attrs)
            
        except KeyError:
            continue            
            
        #for each CG found in sequence extract features
        for t in matcher.finditer(seq,window,len(events)-1-window):
            i=t.start(0)

            pos=i+mapped_start+1 if mapped_strand=='+' else len(events)-i+mapped_start

            m1=rf_events[i-window:i+window+1]
            m2=np.eye(4)[[base_map[x] for x in seq[i-window:i+window+1]]]
            mat=np.hstack([m1,m2])
            
            if np.size(mat)!=21*7:  
                    continue
                    
            #append new data for current position to the lsit of all features 
            features_list.append(mat)
            pos_list.append(pos)
            chr_list.append(mapped_chrom)
            strand_list.append(strand_map[mapped_strand])
            read_names_list.append(read_name) 
    
    
    features_list=np.array(features_list)
    return pos_list, chr_list, strand_list, read_names_list, features_list