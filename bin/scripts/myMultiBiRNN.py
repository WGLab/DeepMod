
#################################
# https://github.com/aymericdamien/TensorFlow-Examples/blob/master/examples/3_NeuralNetworks/bidirectional_rnn.py
#######################################

import tensorflow as tf
from tensorflow.contrib import rnn
import numpy as np

import glob, os, sys, time;

from collections import defaultdict

#batchsize = 4096; # 2048
#batchsize = 3072; 
#
#batchsize = 4096;
#batchsize = 3072;
batchsize = 2048;

class_weights = tf.constant([0.4,0.6])

def mCreateSession(num_input, num_hidden, timesteps):
   num_classes = 2;
   numlayers = 3;

   learning_rate = 0.001

   X = tf.placeholder("float", [None, timesteps, num_input]);
   Y = tf.placeholder("float", [None, num_classes]);

   weights = {'out': tf.Variable(tf.truncated_normal([2*num_hidden, num_classes]))};
   biases = {'out': tf.Variable(tf.truncated_normal([num_classes]))}

   def BiRNN(x, weights, biases):
      x = tf.unstack(x, timesteps, 1);
   
      # https://www.programcreek.com/python/example/102773/tensorflow.contrib.rnn.MultiRNNCell
      # example 14
      #lstm_fw_cell = rnn.BasicLSTMCell(num_hidden, forget_bias=1.0);
      #lstm_bw_cell = rnn.BasicLSTMCell(num_hidden, forget_bias=1.0);
      lstm_fw_cell = rnn.MultiRNNCell([rnn.BasicLSTMCell(num_hidden, forget_bias=1.0) for _ in range(numlayers)]);
      lstm_bw_cell = rnn.MultiRNNCell([rnn.BasicLSTMCell(num_hidden, forget_bias=1.0) for _ in range(numlayers)]);


      try:
         outputs, _, _ = rnn.static_bidirectional_rnn(lstm_fw_cell, lstm_bw_cell, x, dtype=tf.float32);
      except Exception:
         outputs = rnn.static_bidirectional_rnn(lstm_fw_cell, lstm_bw_cell, x, dtype=tf.float32);

      return tf.matmul(outputs[int(timesteps/2)], weights['out']) + biases['out']

   logits = BiRNN(X, weights, biases);
   prediction = tf.nn.softmax(logits)

   mfpred=tf.argmax(prediction,1) 

   ## loss_op = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=Y))
   ###
   loss_op = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=logits, labels=Y))
   #
   #https://stackoverflow.com/questions/40198364/how-can-i-implement-a-weighted-cross-entropy-loss-in-tensorflow-using-sparse-sof/46984951#46984951
   #
   #samples_weights = tf.gather(class_weights, mfpred);
   #loss_op = tf.losses.sparse_softmax_cross_entropy(mfpred, logits, samples_weights) 
   #https://stackoverflow.com/questions/35155655/loss-function-for-class-imbalanced-binary-classifier-in-tensor-flow
   #loss_op = tf.nn.weighted_cross_entropy_with_logits(logits=logits, targets=Y, pos_weight=classes_weights)
   #
   # https://datascience.stackexchange.com/questions/17219/how-to-do-imbalanced-classification-in-deep-learning-tensorflow-rnn
   
   optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate);
   train_op = optimizer.minimize(loss_op);

   correct_pred = tf.equal(tf.argmax(prediction, 1), tf.argmax(Y, 1));
   accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32));

   #auc_op = tf.metrics.auc(Y, prediction)
   #mpre = tf.metrics.precision(Y, prediction)
   #mspf = tf.metrics.recall(Y, prediction)
   #auc_op = tf.metrics.auc(tf.argmax(Y, 1), tf.argmax(prediction, 1))
   auc_op = tf.metrics.auc(Y, prediction)
   mpre = tf.metrics.precision(tf.argmax(Y, 1), tf.argmax(prediction, 1))
   mspf = tf.metrics.recall(tf.argmax(Y, 1), tf.argmax(prediction, 1))
   #mpre = tf.metrics.precision(tf.cast(Y,tf.int32), tf.cast(prediction,tf.int32))
   #mspf = tf.metrics.recall(tf.cast(Y,tf.int32), tf.cast(prediction,tf.int32))
   #auc_op, mpre, mspf = None, None, None

   init = tf.global_variables_initializer();
   init_l = tf.local_variables_initializer()

   saver = tf.train.Saver();
   
   return (init, init_l, loss_op, accuracy, train_op, X, Y, saver, auc_op, mpre, mspf, mfpred)

def train_save_model(filelists, num_input, mhidden, timesteps, moptions):
   #desplay_files = 20
   #desplay_files = 2   # for 700M
   #desplay_files = 30  # for 30M
   training_steps = 4
   training_steps = 10

   init, init_l, loss_op, accuracy, train_op, X, Y, saver, auc_op, mpre, mspf, mfpred = mCreateSession(num_input, mhidden, timesteps)

   desplay_files = len(filelists[0])/100
   if desplay_files<2: desplay_files = 2;
   if desplay_files>10: desplay_files = int(desplay_files/10) * 10; #desplay_files=2
   if desplay_files>100: desplay_files = 100
   file_group_id = [0 for _ in range(len(filelists))];
   sumpsize = 25;

   config = tf.ConfigProto()
   #config.gpu_options.allow_growth = True
   if (timesteps>61 and num_input>50):
      # for 81 of E1m2
      config.gpu_options.per_process_gpu_memory_fraction = 0.6
      # for 81 of P90
      config.gpu_options.per_process_gpu_memory_fraction = 0.5
   else: config.gpu_options.allow_growth = True
   with tf.Session(config=config) as sess:
      if 'modfile' in moptions and (not moptions['modfile']==None) and len(moptions['modfile'])==2:
         new_saver = tf.train.import_meta_graph(moptions['modfile'][0]+'.meta')
         new_saver.restore(sess,tf.train.latest_checkpoint(moptions['modfile'][1])) 
      else:
         sess.run(init);
         sess.run(init_l) 
      start_time = time.time(); start_c_time = time.time();
      io_time = 0;

      for step in range(1, training_steps+1):
         print('===%d=====================step========================%d/%d' % (desplay_files, step, training_steps))
         sys.stdout.flush()
         last_desplay_files_num = -1;
         file_group_id[0] = 0
         while file_group_id[0] < len(filelists[0]):
             io_start_time = time.time();

             featurelist = [[[], []] for _ in range(len(filelists))];
             minsize = None; cur_batch_num = None;
             for ifl in range(len(filelists)):
                if ifl==0:
                   minsize = batchsize * sumpsize
                else: minsize = batchsize * cur_batch_num;
                while len(featurelist[ifl][0])<minsize:
                   if not file_group_id[ifl] < len(filelists[ifl]): 
                      if ifl==0: break;
                      else: file_group_id[ifl] = 0
                   batch_2_x, batch_2_y, _ = getDataFromFile_new(filelists[ifl][file_group_id[ifl]], moptions)
                   if len(batch_2_y)>0:
                      if len(featurelist[ifl][0])==0:
                         featurelist[ifl][0] = batch_2_x
                         featurelist[ifl][1] = batch_2_y
                      else:
                         featurelist[ifl][0] = np.concatenate((featurelist[ifl][0], batch_2_x), axis=0)
                         featurelist[ifl][1] = np.concatenate((featurelist[ifl][1], batch_2_y), axis=0)
                   file_group_id[ifl] += 1;
                if ifl==0:
                   featurelist[ifl][0] = np.array_split(featurelist[ifl][0], int(len(featurelist[ifl][0])/batchsize)) 
                   featurelist[ifl][1] = np.array_split(featurelist[ifl][1], int(len(featurelist[ifl][1])/batchsize))
                   cur_batch_num = len(featurelist[ifl][0])
             if len(featurelist[0][0])<sumpsize*0.8:
                for ifl in range(1, len(filelists)):
                    if len(featurelist[0][0])*batchsize*1.2 < len(featurelist[ifl][0]):
                       featurelist[ifl][0] = featurelist[ifl][0][:int(len(featurelist[0][0])*batchsize*1.2)]
                       featurelist[ifl][1] = featurelist[ifl][1][:int(len(featurelist[0][0])*batchsize*1.2)]
             #
             if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files: msizeprint = ['Read row info', str(len(featurelist[0][0]))]
             for ifl in range(1, len(filelists)):
                if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files: msizeprint.append(str(len(featurelist[ifl][0])))
                featurelist[ifl][0] = np.array_split(featurelist[ifl][0], len(featurelist[0][0]))
                featurelist[ifl][1] = np.array_split(featurelist[ifl][1], len(featurelist[0][0]))
             if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files: print (' '.join(msizeprint))
             io_time += (time.time() - io_start_time)

             ifl=3 if len(featurelist)>3 else len(featurelist)-1
             if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files: 
                sess.run(init_l)
                loss, aucm, acc, p, r = sess.run([loss_op, auc_op[1], accuracy, mpre[1], mspf[1]], feed_dict={X:featurelist[ifl][0][0], Y:featurelist[ifl][1][0]})
                print(">>>Tratin#files "+str(file_group_id[0]+1)+",loss="+"{:.3f}".format(loss)+",AUC="+"{:.3f}".format(aucm)+",acc="+"{:.3f}".format(acc)+",p="+"{:.3f}".format(p)+",r="+"{:.3f}".format(r)+(" Comsuming time: %d(current=%d) IO=%d(%.3f)" % (time.time()-start_time, time.time()-start_c_time, io_time, io_time/float(time.time()-start_time))));

                sys.stdout.flush()
                start_c_time = time.time();

             for subi in range(len(featurelist[0][0])):
                for ifl in range(len(filelists)):
                   sess.run(train_op, feed_dict={X:featurelist[ifl][0][subi], Y:featurelist[ifl][1][subi]})

             ifl=3 if len(featurelist)>3 else len(featurelist)-1
             if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files:
                # sess.run(init_l)
                # loss, aucm, acc, p, r = sess.run([loss_op, auc_op[1], accuracy, mpre[1], mspf[1]], feed_dict={X:featurelist[ifl][0][0], Y:featurelist[ifl][1][0]})
                # sess.run(init_l)
                # if ifl>0: ifl -= 1
                # loss2, aucm2, acc2, p2, r2 = sess.run([loss_op, auc_op[1], accuracy, mpre[1], mspf[1]], feed_dict={X:featurelist[ifl][0][0], Y:featurelist[ifl][1][0]})
                # print("   Train#files "+str(file_group_id[0]+1)+", loss="+"{:.3f}".format(loss)+", AUC="+"{:.3f}".format(aucm)+", acc="+"{:.3f}".format(acc)+", p="+"{:.3f}".format(p)+", r="+"{:.3f}".format(r)+(' %d/%d/%d %d/%d, %d:%d(%d:%d) %d/%d' % (len(featurelist[ifl][0])*batchsize, len(featurelist[ifl][0][0][0]), len(featurelist[ifl][0][0][0][0]), len(featurelist[ifl][1])*batchsize, len(featurelist[ifl][1][0][0]), np.sum(featurelist[ifl][1][0][:,0]), np.sum(featurelist[ifl][1][0][:,1]), np.sum(featurelist[ifl][1][0][:batchsize,1]), np.sum(featurelist[ifl+1][1][0][:batchsize,1]), len(featurelist[ifl][0]), len(featurelist[ifl][1])))+"\n                  loss="+"{:.3f}".format(loss2)+", AUC="+"{:.3f}".format(aucm2)+", acc="+"{:.3f}".format(acc2)+", p="+"{:.3f}".format(p2)+", r="+"{:.3f}".format(r2) + (", %d/%d=%d" % (file_group_id[0], len(filelists[0]), int(file_group_id[0]*100/float(len(filelists[0]))) ) ) );
                # sys.stdout.flush()
                last_desplay_files_num = (file_group_id[0]+1) - ((file_group_id[0]+1) % desplay_files)

             if 49.5<int(file_group_id[0]*100/float(len(filelists[0])))<50.5:
                savp = '.50'
                if (not os.path.isdir(moptions['outFolder']+str(step-1)+savp)):
                  os.system('mkdir -p '+moptions['outFolder']+str(step-1)+savp);
                #if not os.path.isfile(moptions['outFolder']+str(step-1)+savp+'/'+moptions['FileID']+'.meta'):
                #   saver.save(sess, moptions['outFolder']+str(step-1)+savp+'/'+moptions['FileID']);
                saver.save(sess, moptions['outFolder']+str(step-1)+savp+'/'+moptions['FileID']);
         if (not os.path.isdir(moptions['outFolder']+str(step))): 
            os.system('mkdir -p '+moptions['outFolder']+str(step));
         saver.save(sess, moptions['outFolder']+str(step)+'/'+moptions['FileID']);
      print("Training Finished!")

      return (accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred)

def getTFiles1(folder1, moptions):
   t1files = glob.glob(os.path.join(folder1, "*.xy.gz"))
   if moptions['recursive']==1:
      t1files.extend(glob.glob(os.path.join(folder1, "*/*.xy.gz")))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*.xy.gz"))); 
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*.xy.gz"))); 
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*/*.xy.gz")));
   print("Get folder1");
   print(t1files.__sizeof__(), len(t1files))
   if moptions['test'][0] == '0':
      if moptions['test'][1]>0.5:
         t1files = t1files[:int(len(t1files)*moptions['test'][1])]
      else: t1files = t1files[-int(len(t1files)*moptions['test'][1]):]
   print('Sizeinfo: %s sizeof=%d len=%d' % (folder1, t1files.__sizeof__(), len(t1files)))
   sys.stdout.flush();

   return t1files

def getTFiles(folder1, folder2, moptions):
   t1files = glob.glob(os.path.join(folder1, "*.xy.gz")); #print(t1files.__sizeof__(), len(t1files))
   if moptions['recursive']==1:
      t1files.extend(glob.glob(os.path.join(folder1, "*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
   print("Get folder1");
   print(t1files.__sizeof__(), len(t1files))
   if moptions['test'][0] == '0':
      if moptions['test'][1]>0.5:
         t1files = t1files[:int(len(t1files)*moptions['test'][1])] 
      else: t1files = t1files[-int(len(t1files)*moptions['test'][1]):]
   print(t1files.__sizeof__(), len(t1files))
   sys.stdout.flush();

   if folder2==None: t2files = []
   else:
      t2files = glob.glob(os.path.join(folder2, "*.xy.gz"))
      if moptions['recursive']==1:
         t2files.extend(glob.glob(os.path.join(folder2, "*/*.xy.gz")))
         t2files.extend(glob.glob(os.path.join(folder2, "*/*/*.xy.gz")))
         t2files.extend(glob.glob(os.path.join(folder2, "*/*/*/*.xy.gz")))
         t2files.extend(glob.glob(os.path.join(folder2, "*/*/*/*/*.xy.gz")))
      print("Get folder2");
      print(t2files.__sizeof__(), len(t2files))
      if moptions['test'][0] == '0':
         if moptions['test'][1]>0.5:
            t2files = t2files[:int(len(t2files)*moptions['test'][1])]
         else: t2files = t2files[-int(len(t2files)*moptions['test'][1]):]
      print(t2files.__sizeof__(), len(t2files))
      sys.stdout.flush();
   return t1files, t2files
def getDataFromFile(fn, moptions):
   mdata = np.loadtxt(fn, dtype=np.float32)
   t0, ty, tx = np.split(mdata, [0,2], axis=1);
   return (tx, ty, None)
def getDataFromFile_new(fn, moptions, mfind0ld=None):
   mdata = np.loadtxt(fn, dtype=np.float32)
   t0, ty, tx = np.split(mdata, [1,3], axis=1);

   if moptions['test'][0] in ['-', '+']:
      t0 = t0.astype(int)

   m_data = []; m_y = [];
   if not mfind0ld==None:
      pos_to_file_dict = defaultdict();  preind = 0
      mfind0ldkeys = sorted(list(mfind0ld.keys()));
   for mind in range(len(ty)):
      if not mfind0ld==None:
         if preind<len(mfind0ldkeys) and mind == mfind0ldkeys[preind]:
            pos_to_file_dict[len(m_y)] = mfind0ld[ mfind0ldkeys[preind] ]
            preind += 1

      if (ty[mind][0]<0.01 and ty[mind][1]<0.01): continue;
      if (moptions['test'][0]=='-' and moptions['test'][1]<t0[mind]<moptions['test'][2]) or \
         (moptions['test'][0]=='+' and (not moptions['test'][1]<t0[mind]<moptions['test'][2])):
         continue;
      ## for 61; does not work. memory issue;
      m_y.append(ty[mind])
      m_data.append(tx[(mind-int(moptions['windowsize']/2)):(mind+int(moptions['windowsize']/2)+1)])
      ## for 7 to 51;
      #m_y.append(np.copy(ty[mind]))
      #m_data.append(np.copy(tx[(mind-int(moptions['windowsize']/2)):(mind+int(moptions['windowsize']/2)+1)]))
   if not mfind0ld==None:
      file_to_pos_dict = defaultdict();
      ptofkeys = sorted(list(pos_to_file_dict.keys()))
      for npk_ind in range(len(ptofkeys)):
         if (npk_ind+1<len(ptofkeys) and ptofkeys[npk_ind+1]-ptofkeys[npk_ind]<500) or len(m_y)-ptofkeys[npk_ind]<500: continue;
 
         file_to_pos_dict[ pos_to_file_dict[ptofkeys[npk_ind]] ] = [ptofkeys[npk_ind], (ptofkeys[npk_ind+1] if npk_ind+1<len(ptofkeys) else len(m_y))]

   if len(m_data)>0:
      m_data = np.reshape(m_data, (len(m_data), len(m_data[0]), len(m_data[0][0])))
      m_y = np.reshape(m_y, (len(m_y), len(m_y[0]))).astype(int)

   if not mfind0ld==None:
      return (m_data, m_y, file_to_pos_dict);
   else: return (m_data, m_y, None)


def getGZFilePos(gzfile):
   mfind = defaultdict()
   with open(gzfile[:-len('.gz')]+'.ind', 'r') as mr:
      line = mr.readline()
      while line:
          line = line.strip();
          lsp = line.split();
          if len(lsp)>1:
             mfind[int(lsp[0])] = lsp[1]
          line = mr.readline()
   return mfind


def mPred(mfbase, mffolder, accuracy, X, Y, test_gzfile2, pf, num_input, auc_op, mpre, mspf, init_l, mfpred, timesteps, moptions):
   config = tf.ConfigProto()
   config.gpu_options.allow_growth = True
   with tf.Session(config=config) as sess:
      new_saver = tf.train.import_meta_graph(mfbase+'.meta')
      new_saver.restore(sess,tf.train.latest_checkpoint(mffolder))

      pfwriter = open(pf, 'w');
      for test_gzfile in test_gzfile2:
         for test_fn_ind in range(len(test_gzfile)):
            singlefeatfiledict = getGZFilePos(test_gzfile[test_fn_ind])
            test_gzfeature, test_gzlabel, file_to_pos_dict = getDataFromFile(test_gzfile[test_fn_ind], moptions) if is_new==0 else getDataFromFile_new(test_gzfile[test_fn_ind], moptions, singlefeatfiledict)

            if len(test_gzlabel)==0: continue;
            detailwriter = open(test_gzfile[test_fn_ind]+'.pred', 'w')
            file_to_pos_dict_keys = sorted(list(file_to_pos_dict.keys()))
            for test_file in file_to_pos_dict_keys: 
               test_feature = test_gzfeature[file_to_pos_dict[test_file][0]:file_to_pos_dict[test_file][1]]
               test_label   = test_gzlabel[file_to_pos_dict[test_file][0]:file_to_pos_dict[test_file][1]]

               sess.run(init_l)
               testacc, p,r,aucm, mfpred_output = sess.run([accuracy, mpre[1], mspf[1], auc_op[1], mfpred], feed_dict={X:test_feature, Y:test_label})
               print("Testing accuracy=" + "{:.3f}".format(testacc)+", p=" + "{:.3f}".format(p)+", r=" + "{:.3f}".format(r)+", auc=" + "{:.3f}".format(aucm))
               pfwriter.write(('%.3f %.3f %.3f %.3f %s\n' % (testacc, p, r, aucm, test_file)))
               sys.stdout.flush()
               
               detailwriter.write(test_file+' '+str(len(test_feature))+' '+''.join(list(map(str, mfpred_output)))+'\n')
               detailwriter.flush()
            detailwriter.close()
            pfwriter.flush()
      pfwriter.close();

def pred_prepare(moptions, test_file, accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred):
    mPred(moptions['modfile'][0], moptions['modfile'][1], accuracy, X, Y, test_file, moptions['outFolder']+moptions['FileID']+'_mpred.txt', moptions['fnum'], auc_op, mpre, mspf, init_l, mfpred, moptions['windowsize'], moptions)

def mMult_RNN_LSTM_train(moptions):

   filegroups = moptions['wrkBase'].split(';')
   for i in range(len(filegroups)):
      filegroups[i] = filegroups[i].split(',')

   print(filegroups)

   filelists = [[] for _ in range(len(filegroups))]
   for i in range(len(filegroups)):
      for fgj in range(len(filegroups[i])):
         if not len(filegroups[i][fgj])>0: continue
         filelists[i].extend(getTFiles1(filegroups[i][fgj], moptions))
   mostnum, mostid = 0, -1;
   np.random.seed(1)
   for i in range(len(filelists)):
      np.random.shuffle(filelists[i])
      if len(filelists[i])>mostnum:
         mostnum = len(filelists[i])
         mostid = i;

   if 'modfile' in moptions and (not moptions['modfile']==None):
      if moptions['modfile'].rfind('/')==-1:
         moptions['modfile'] = [moptions['modfile'], './']
      else:
         moptions['modfile'] = [moptions['modfile'], moptions['modfile'][:moptions['modfile'].rfind('/')+1]]
 
   if not mostid==0:
      filelists[mostid], filelists[0] = filelists[0], filelists[mostid]

   accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred = train_save_model(filelists, moptions['fnum'], moptions['hidden'], moptions['windowsize'], moptions)

def pred_entry(moptions):

   tfiles = getTFiles(moptions['wrkBase'], None, moptions)

   init, init_l, loss_op, accuracy, train_op, X, Y, saver, auc_op, mpre, mspf, mfpred = mCreateSession(moptions['fnum'], moptions['hidden'], moptions['windowsize'])

   if moptions['modfile'].rfind('/')==-1:
      moptions['modfile'] = [moptions['modfile'], './']
   else:
      moptions['modfile'] = [moptions['modfile'], moptions['modfile'][:moptions['modfile'].rfind('/')+1]]

   pred_prepare(moptions, tfiles, accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred)



