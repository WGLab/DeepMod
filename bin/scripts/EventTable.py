
import os,sys
import numpy as np
import h5py


def get_extreme_N(m_signal_dif, n_splits, p_signal_start, p_signal_end, moptions, sp_param):
   cu_region_sort_pos = m_signal_dif[int(p_signal_start-sp_param['min_signal_num']+0.5):int(p_signal_end-sp_param['min_signal_num']+0.5)].argsort()[::-1]+p_signal_start;
   m_nb_pos = set();
   # print n_splits, type(n_splits), p_signal_start, type(p_signal_start), p_signal_end, type(p_signal_end), sp_param['min_signal_num'], type( sp_param['min_signal_num']), type(p_signal_start+sp_param['min_signal_num']-1)
   m_nb_pos.update(range(p_signal_start, int(p_signal_start+sp_param['min_signal_num']-0.5)));
   m_nb_pos.update(range(int(p_signal_end-sp_param['min_signal_num']+1.5), p_signal_end));
   split_points_list = []
   for c_pos in cu_region_sort_pos:
      if c_pos not in m_nb_pos:
         split_points_list.append(c_pos);
         if (len(split_points_list)==n_splits): break;
         m_nb_pos.update(range(c_pos-sp_param['min_signal_num']+1, c_pos+sp_param['min_signal_num']+1));
   return sorted(split_points_list);

def getEvent_Info(moptions, sp_param, events_data):
   event_info = []
   sp_param['min_signal_num'] = 4;

   signal_sum = np.cumsum(np.insert(np.round(sp_param['raw_signals']/50.0,5), 0, 0));
   m_signal_dif = np.abs(signal_sum[sp_param['min_signal_num']:-sp_param['min_signal_num']]*2 - signal_sum[:-2*sp_param['min_signal_num']] - signal_sum[2*sp_param['min_signal_num']:])
   #print (sp_param['raw_signals'][:20]);
   #print (np.round(sp_param['raw_signals']/50.0,5)[:20]);
   #print (signal_sum[:20]);
   #print (m_signal_dif[:20])
   # sp_param['fq_seq'] = fq_data[1]
   last_ev_i = 0;
   last_signal_i = events_data[0]['start'];
   fq_seq_i = 2;
   c_move_num = 1
   incrrt_event_list = []
   for ev_i in range(1, len(events_data)):
      if (events_data['move'][ev_i])==0:
         pass;
      else:
         c_move_num += events_data['move'][ev_i]
         split_points = get_extreme_N(m_signal_dif, c_move_num-1, last_signal_i, events_data[ev_i]['start']+events_data[ev_i]['length'], moptions, sp_param);
         #print c_move_num-1, last_signal_i, ev_i, events_data[ev_i]['start']+events_data[ev_i]['length'], split_points
         #for s_i in range(last_signal_i, events_data[ev_i]['start']+events_data[ev_i]['length']): 
         #    if s_i in split_points:
         #       print '|',
         #    print sp_param['raw_signals'][s_i],
         #print '';
         for c_m_i in range(c_move_num-1):
            if c_m_i < len(split_points): 
               h_m_i = c_m_i;
               c_e_p = split_points[h_m_i]
            else: 
               h_m_i = len(split_points)-1
               c_e_p = last_signal_i + sp_param['min_signal_num']
               incrrt_event_list.append(len(event_info));

            c_mnn = np.mean(sp_param['raw_signals'][last_signal_i:c_e_p]);
            c_std = np.std(sp_param['raw_signals'][last_signal_i:c_e_p]);
            c_start = last_signal_i;
            c_length = c_e_p - last_signal_i;
            c_mode = sp_param['fq_seq'][fq_seq_i-2:fq_seq_i+3];
            event_info.append((c_mnn, c_std, c_start, c_length, c_mode))

            last_signal_i = split_points[h_m_i]
            fq_seq_i += 1;
 
         c_move_num = 1;
   ev_i = len(events_data)-1 
   c_e_p = events_data[ev_i]['start'] + events_data[ev_i]['length']
   c_mnn = np.mean(sp_param['raw_signals'][last_signal_i:c_e_p]);
   c_std = np.std(sp_param['raw_signals'][last_signal_i:c_e_p]);
   c_start = last_signal_i;
   c_length = c_e_p - last_signal_i;
   c_mode = sp_param['fq_seq'][fq_seq_i-2:fq_seq_i+3];
   event_info.append((c_mnn, c_std, c_start, c_length, c_mode))

   event_info = np.array(event_info, dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64), ('model_state', 'U5')])
   #c_seq = ''.join([event_model_state[2] for event_model_state in event_info['model_state'] ] )
 
   #print '\n' 
   for c_ev_i in incrrt_event_list:
      #print c_ev_i, event_info[c_ev_i-1]['start'], event_info[c_ev_i-1]['length'], event_info[c_ev_i]['start'], event_info[c_ev_i]['length'], event_info[c_ev_i+1]['start'], event_info[c_ev_i+1]['length']
      h_2 = int((event_info[c_ev_i+1]['length'] + event_info[c_ev_i+1]['start'] - event_info[c_ev_i]['start'] )/2+0.2)
      event_info[c_ev_i]['length'] = h_2
      event_info[c_ev_i+1]['start'] = event_info[c_ev_i]['start'] + event_info[c_ev_i]['length']
      event_info[c_ev_i+1]['length'] = event_info[c_ev_i+1]['length'] - h_2
      #print '\t', c_ev_i, event_info[c_ev_i-1]['start'], event_info[c_ev_i-1]['length'], event_info[c_ev_i]['start'], event_info[c_ev_i]['length'], event_info[c_ev_i+1]['start'], event_info[c_ev_i+1]['length']

   #for c_ev_i in range(len(event_info)):
   #   print c_ev_i, event_info[c_ev_i]['start'], event_info[c_ev_i]['length'], ':', 
   #   for s_i in range(event_info[c_ev_i]['start'], event_info[c_ev_i]['start']+event_info[c_ev_i]['length']):
   #       pass # print sp_param['raw_signals'][s_i],
   #   print ''
 
   #msi = 50;
   #print (c_seq[:msi])
   #print (sp_param['fq_seq'][2:(msi+2)])
   #print (c_seq[-msi:])
   #print (sp_param['fq_seq'][-(msi+2):-2])
   #print len(events_data), len(event_info), len(sp_param['fq_seq'])
   #ei_i = 0;
   #for ev_i in range(0, len(events_data)):
   #   if (events_data[ev_i]['move']>0):
   #      print ("%d/%s %d-%d vs %d-%d %s=%s%s" % (ev_i, ei_i,events_data[ev_i]['start'], events_data[ev_i]['start']+events_data[ev_i]['length'],  event_info[ei_i]['start'], event_info[ei_i]['start']+event_info[ei_i]['length'], events_data[ev_i]['model_state'][2],event_info[ei_i]['model_state'][2],sp_param['fq_seq'][ei_i+2]))
   #      ei_i += events_data[ev_i]['move']

   return event_info


if __name__=='__main__':
   moptions = {}
   sp_param = {}

   exple_data = ['/home/liuq1/project/DeepNanoRepeat/scripts/fortest/f6343e53-9454-41ae-8398-7be6e1b7557d.fast5', \
                 'data/alb231/S_053119TrainSeq3ctrloligoSpeIcut/workspace/pass/0/000a7916-373c-4cc3-a3f2-6bed205b09cb.fast5', \
                 'data/alb231/S_053119TrainSeq3ctrloligoSpeIcut/workspace/pass/0/00264c38-4945-4263-ae0d-253e6c6a39ba.fast5', \
                 'data/alb231/S_053119TrainSeq3ctrloligoSpeIcut/workspace/pass/0/0039f109-46ac-4a81-883d-b55900924dd4.fast5', \
                 'data/alb231/S_053119TrainSeq3ctrloligoSpeIcut/workspace/pass/0/0045bf1d-d7be-44b1-9b6c-9bb76a634e0f.fast5' \
                ]

   sp_param['f5reader'] = h5py.File(sys.argv[1] if len(sys.argv)>1 else exple_data[0], 'r');   

   fq_str = '/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
   ev_str = '/Analyses/Basecall_1D_000/BaseCalled_template/Events'
   fq_str = '/Analyses/Basecall_1D_001/BaseCalled_template/Fastq'
   ev_str = '/Analyses/Basecall_1D_001/BaseCalled_template/Events'
   sg_str = '/Raw/Reads/'

   sp_param['fq_seq'] = sp_param['f5reader'][fq_str][()].split('\n')[1];
   sp_param['raw_signals'] = sp_param['f5reader'][sg_str].values()[0]['Signal'].value
   events_data = sp_param['f5reader'][ev_str].value;

   getEvent_Info(moptions, sp_param, events_data)

   sp_param['f5reader'].close();   

