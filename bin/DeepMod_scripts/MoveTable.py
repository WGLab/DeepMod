
import os,sys
import numpy as np
import h5py


def getMove_Info(moptions, sp_param, move_data):
    '''
    sp_param.keys: fq_seq, raw_signals, first_sample_template, duration_template
    '''

    #sp_param['first_sample_template'] = sp_param['f5reader']['/Analyses/Segmentation_001/Summary/segmentation'].attrs['first_sample_template']
    #sp_param['duration_template'] = sp_param['f5reader']['/Analyses/Segmentation_001/Summary/segmentation'].attrs['duration_template']

    seg = "Segmentation_" + moptions['basecall_1d'].split('_')[-1]
    attr_path = '/'.join(['', 'Analyses', seg, 'Summary', 'segmentation'])
    #mv_str = '/'.join(['', 'Analyses', moptions['basecall_1d'], moptions['basecall_2strand'], 'Move'])
    sp_param['first_sample_template'] = sp_param['f5reader'][attr_path].attrs['first_sample_template']
    sp_param['duration_template'] = sp_param['f5reader'][attr_path].attrs['duration_template']
    #move_data = sp_param['f5reader'][mv_str][()]
    nrow = len(sp_param['fq_seq']) # row number of event_info; equals to the base number
    nsig = len(sp_param['raw_signals'])
    first = int(sp_param['first_sample_template'])
    duration = int(sp_param['duration_template'])
    move_info = np.empty([nrow], dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64), ('model_state', 'U5')])
    effect_sig_index = list(range(first, nsig))
    pivot = first
    seg_count = 0 #which segmentation
    for i in range(1, len(move_data)):
        if move_data[i] == 1:
            move_info[seg_count]['mean'] = np.mean(sp_param['raw_signals'][pivot:(2*i + first)])
            move_info[seg_count]['length'] = 2*i + first - pivot
            move_info[seg_count]['stdv'] = np.std(sp_param['raw_signals'][pivot:(2*i + first)])
            move_info[seg_count]['start'] = pivot
            if seg_count == 0:
                move_info[seg_count]['model_state'] = 'N'*2 + sp_param['fq_seq'][seg_count:seg_count+3]
            elif seg_count == 1:
                move_info[seg_count]['model_state'] = 'N' + sp_param['fq_seq'][seg_count-1:seg_count+3]
            elif seg_count == nrow-2:
                move_info[seg_count]['model_state'] = sp_param['fq_seq'][seg_count-2:seg_count+2] + 'N'
            else:
                move_info[seg_count]['model_state'] = sp_param['fq_seq'][seg_count-2 : seg_count+3]
            pivot = 2*i + first
            seg_count += 1
    move_info[seg_count]['mean'] = np.mean(sp_param['raw_signals'][pivot:nsig])
    move_info[seg_count]['length'] = nsig - pivot
    move_info[seg_count]['stdv'] = np.std(sp_param['raw_signals'][pivot:nsig])
    move_info[seg_count]['start'] = pivot
    move_info[seg_count]['model_state'] = sp_param['fq_seq'][seg_count-2:seg_count+1] + 'N'*2
    return move_info


if __name__=='__main__':
   moptions = {}
   sp_param = {}

   exple_data = ['/mnt/isilon/wang_lab/shared/temp_shared/fast5_move/IBDUCAL377261L_20170201_FNfab41074_MN17640_mux_scan_X209_66786_ch12_read102_strand.fast5', \
                 '/mnt/isilon/wang_lab/shared/temp_shared/fast5_move/IBDUCAL377261L_20170201_FNfab41074_MN17640_mux_scan_X209_66786_ch47_read165_strand.fast5', \
                 '/mnt/isilon/wang_lab/shared/temp_shared/fast5_move/IBDUCAL377261L_20170201_FNfab41074_MN17640_mux_scan_X209_66786_ch48_read38_strand.fast5', \
                 '/mnt/isilon/wang_lab/shared/temp_shared/fast5_move/IBDUCAL377261L_20170201_FNfab41074_MN17640_mux_scan_X209_66786_ch52_read12_strand.fast5' \
                 ]

   sp_param['f5reader'] = h5py.File(exple_data[0], 'r');


   fq_str = '/Analyses/Basecall_1D_001/BaseCalled_template/Fastq'
   mv_str = '/Analyses/Basecall_1D_001/BaseCalled_template/Move'
   sg_str = '/Raw/Reads/'

   sp_param['fq_seq'] = sp_param['f5reader'][fq_str][()].splitlines()[1]
   k = list(sp_param['f5reader'][sg_str].keys())[0]
   sp_param['raw_signals'] = sp_param['f5reader'][sg_str][k]['Signal'][()]
   sp_param['first_sample_template'] = sp_param['f5reader']['/Analyses/Segmentation_001/Summary/segmentation'].attrs['first_sample_template']
   sp_param['duration_template'] = sp_param['f5reader']['/Analyses/Segmentation_001/Summary/segmentation'].attrs['duration_template']
   move_data = sp_param['f5reader'][mv_str][()]

   getEvent_Info(moptions, sp_param, move_data)

   sp_param['f5reader'].close();
