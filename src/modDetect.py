from collections import defaultdict, ChainMap

import csv, time, itertools, copy, h5py, time

import datetime, os, shutil, argparse, sys, random

import multiprocessing as mp
import numpy as np

import numpy.lib.recfunctions as rf

from pathlib import Path
from tensorflow import keras
import tensorflow as tf

from .  import guppy
from .  import tombo

from .utils import *

model_dict={'guppy_na12878':'models/guppy/test_model', 'tombo_na12878':'models/tombo/test_model'}

def get_model(model):
    if model in model_dict:
        dirname = os.path.dirname(__file__)
        return os.path.join(dirname, model_dict[model])
        
    elif os.path.exists(model) and os.path.isdir(model_path):
        return model
     
    else:
        return None
    
def detect(args):
    gpus = tf.config.experimental.list_physical_devices('GPU')
    
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
        
    strand_map={0:'+', 1:'-'}
    
    f5files, params, read_info, job_number = args
    
    threshold=0.5
    
    output=os.path.join(params['output'],'intermediate_files', 'part_%d' %job_number)
    
    model=keras.models.load_model(get_model(params['model']))
    
    with open(output, 'w') as outfile: 
    
        for f5_chunk in split_list(f5files, 1000):
            
            if params['basecaller']=='guppy':
                pos_list, chr_list, strand_list, read_names_list, features_list = guppy.getFeatures(f5_chunk, params, read_info)
            
            else:
                pos_list, chr_list, strand_list, read_names_list, features_list = tombo.getFeatures(f5_chunk, params)
                
            if len(features_list)==0:
                continue
            pred_list=model.predict(features_list)
            
            for i in range(len(pos_list)):
                pos, chrom, strand, read_name = pos_list[i], chr_list[i], strand_list[i], read_names_list[i]
                outfile.write('%s\t%s\t%d\t%s\t%.4f\t%d\n' %(read_name, chrom, pos, strand_map[strand], pred_list[i], 1 if pred_list[i]>=threshold else 0))
            outfile.flush()
            os.fsync(outfile.fileno())
    
    return output

def per_read_predict(params):
    
    temp_folder=os.path.join(params['output'],'intermediate_files')
    os.makedirs(temp_folder, exist_ok=True)
    
    
    f5files = list(Path(params['fast5']).rglob("*.fast5"))
    files_per_process=len(f5files)//params['threads'] + 1
    print('Number of files: %d\n' %len(f5files))
    
    pool = mp.Pool(processes=params['threads'])
    
    if params['basecaller']=='guppy':
        read_info=guppy.process_bam(params, pool)
    else:
        read_info=None
    
    job_counter=itertools.count(start=1, step=1)
    
    print('%s: Starting Methylation Detection.' %str(datetime.datetime.now()))
    
    res=pool.imap_unordered(detect, zip(split_list(f5files, files_per_process), itertools.repeat(params), itertools.repeat(read_info), job_counter))
    
    file_list=[file_name for file_name in res]
    
    output=os.path.join(params['output'], '%s.per_read' %params['file_name'])
    with open(output,'wb') as outfile:
        outfile.write(b'read\tchrom\tpos\tstrand\tscore\tmeth\n')
        for f in file_list:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, outfile)
            #os.remove(f)
    
    #shutil.rmtree(temp_folder)
    
    pool.close()
    pool.join()
    
    print('%s: Finishing DeepMod Methylation Detection.' %str(datetime.datetime.now()))
    
    return output


def per_site_detect(read_pred_file, params):
    threshold=0.5
    output_raw=os.path.join(params['output'], '%s.per_site_raw' %params['file_name'])
    
    output=os.path.join(params['output'], '%s.per_site' %params['file_name'])
    
    per_site_pred={}
    
    with open(read_pred_file,'r') as read_file:
        for line in read_file.readlines()[1:]:
            read, chrom, pos, strand, score, meth=line.rstrip('\n').split('\t')
            
            if (chrom, pos, strand) not in per_site_pred:
                per_site_pred[(chrom, pos, strand)]=[0,0]
            
            per_site_pred[(chrom, pos, strand)][int(meth)]+=1
            
    with open(output_raw,'w') as outfile:
        outfile.write('#chrom\tstart\tend\tstrand\tcov\tmeth\tpercent\tprediction\n')
        
        for x,y in per_site_pred.items():
            p=100*y[1]/sum(y)
            outfile.write('%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\n' %(x[0],int(x[1])-1,x[1],x[2], sum(y), y[1], p, 1 if p>=threshold else 0))
    
    run_cmd('bedtools sort -i %s -header > %s' %(output_raw, output))
    
    os.remove(output_raw)
    
    return output