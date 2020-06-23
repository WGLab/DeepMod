

import tensorflow as tf
from tensorflow.contrib import rnn
import numpy as np

import math
import glob, os, sys, time;

from collections import defaultdict

batchsize = 2048;

# different class weights for unbalanced data
class_weights = tf.constant([0.1,0.9])

#
# create a RNN with LSTM
# define performance evaluation operation
#
def mCreateSession(num_input, num_hidden, timesteps, moptions):
   # two classes only
   num_classes = 2;
   # the number of layers
   numlayers = 3;
   # learning rate
   learning_rate = 0.001

   # define input and output
   X = tf.placeholder("float", [None, timesteps, num_input]);
   Y = tf.placeholder("float", [None, num_classes]);

   # for last layers
   weights = {'out': tf.Variable(tf.truncated_normal([2*num_hidden, num_classes]))};
   biases = {'out': tf.Variable(tf.truncated_normal([num_classes]))}

   # define a bidirectional RNN
   def BiRNN(x, weights, biases):
      x = tf.unstack(x, timesteps, 1);

      # define the LSTM cells
      lstm_fw_cell = rnn.MultiRNNCell([rnn.BasicLSTMCell(num_hidden, forget_bias=1.0) for _ in range(numlayers)]);
      lstm_bw_cell = rnn.MultiRNNCell([rnn.BasicLSTMCell(num_hidden, forget_bias=1.0) for _ in range(numlayers)]);

      # define bidirectional RNN
      try:
         outputs, _, _ = rnn.static_bidirectional_rnn(lstm_fw_cell, lstm_bw_cell, x, dtype=tf.float32);
      except Exception:
         outputs = rnn.static_bidirectional_rnn(lstm_fw_cell, lstm_bw_cell, x, dtype=tf.float32);

      # define output layer
      if moptions['outputlayer'] in ['sigmoid']:
         return tf.contrib.layers.fully_connected(outputs[int(timesteps/2)], num_outputs=num_classes, activation_fn=tf.nn.sigmoid);
      else:
         return tf.matmul(outputs[int(timesteps/2)], weights['out']) + biases['out']

   # get prediction
   logits = BiRNN(X, weights, biases);
   prediction = tf.nn.softmax(logits)

   mfpred=tf.argmax(prediction,1)

   ## with different class-weights or not
   if 'unbalanced' in moptions and (not moptions['unbalanced']==None) and moptions['unbalanced']==1:  # class_weights
      loss_op = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=tf.multiply(logits, class_weights), labels=Y))
   else:
      loss_op = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=logits, labels=Y))
   #

   # for optimizer
   optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate);
   train_op = optimizer.minimize(loss_op);

   # get accuracy
   correct_pred = tf.equal(tf.argmax(prediction, 1), tf.argmax(Y, 1));
   accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32));

   # AUC
   auc_op = tf.metrics.auc(Y, prediction)
   # precision
   mpre = tf.metrics.precision(tf.argmax(Y, 1), tf.argmax(prediction, 1))
   # recall
   mspf = tf.metrics.recall(tf.argmax(Y, 1), tf.argmax(prediction, 1))

   # initialization of variables
   init = tf.global_variables_initializer();
   init_l = tf.local_variables_initializer()

   saver = tf.train.Saver();

   return (init, init_l, loss_op, accuracy, train_op, X, Y, saver, auc_op, mpre, mspf, mfpred)

#
# train a model and save it.
#
def train_save_model(filelists, num_input, mhidden, timesteps, moptions):
   training_steps = 4
   #training_steps = 40

   init, init_l, loss_op, accuracy, train_op, X, Y, saver, auc_op, mpre, mspf, mfpred = mCreateSession(num_input, mhidden, timesteps, moptions)

   # display step
   desplay_files = len(filelists[0])/100
   if desplay_files<2: desplay_files = 2;
   if desplay_files>10: desplay_files = int(desplay_files/10) * 10; #desplay_files=2
   if desplay_files>100: desplay_files = 100
   file_group_id = [0 for _ in range(len(filelists))];
   sumpsize = 25;

   # for configuration
   config = tf.ConfigProto()
   if (timesteps>61 and num_input>50):
      config.gpu_options.per_process_gpu_memory_fraction = 0.5
   else: config.gpu_options.allow_growth = True
   with tf.Session(config=config) as sess:
      # initialization
      sess.run(init);
      sess.run(init_l)
      start_time = time.time(); start_c_time = time.time();
      io_time = 0;

      # for each epoch
      for step in range(1, training_steps+1):
         print('===%d=====================step========================%d/%d' % (desplay_files, step, training_steps))
         sys.stdout.flush()
         last_desplay_files_num = -1;
         file_group_id[0] = 0
         while file_group_id[0] < len(filelists[0]):
             io_start_time = time.time();

             # for each input groups.
             # usually two groups: one positive group and one negative group
             # might also one group containing both positive and negative labelling data
             featurelist = [[[], []] for _ in range(len(filelists))];
             minsize = None; cur_batch_num = None;
             # get data from all groups until 'minsize' data is loaded.
             for ifl in range(len(filelists)):
                if ifl==0:
                   minsize = batchsize * sumpsize
                else: minsize = batchsize * cur_batch_num;
                while len(featurelist[ifl][0])<minsize:
                   if not file_group_id[ifl] < len(filelists[ifl]):
                      if ifl==0: break;
                      else: file_group_id[ifl] = 0
                   # get more data
                   batch_2_x, batch_2_y, _ = getDataFromFile_new(filelists[ifl][file_group_id[ifl]], moptions)
                   if len(batch_2_y)>0:
                      if len(featurelist[ifl][0])==0:
                         featurelist[ifl][0] = batch_2_x
                         featurelist[ifl][1] = batch_2_y
                      else:
                         # merge current loading data with previously loading data
                         featurelist[ifl][0] = np.concatenate((featurelist[ifl][0], batch_2_x), axis=0)
                         featurelist[ifl][1] = np.concatenate((featurelist[ifl][1], batch_2_y), axis=0)
                   file_group_id[ifl] += 1;
                # split for small groups for training
                if ifl==0:
                   featurelist[ifl][0] = np.array_split(featurelist[ifl][0], int(len(featurelist[ifl][0])/batchsize))
                   featurelist[ifl][1] = np.array_split(featurelist[ifl][1], int(len(featurelist[ifl][1])/batchsize))
                   cur_batch_num = len(featurelist[ifl][0])
             if len(featurelist[0][0])<sumpsize*0.8:
                for ifl in range(1, len(filelists)):
                    if len(featurelist[0][0])*batchsize*1.2 < len(featurelist[ifl][0]):
                       featurelist[ifl][0] = featurelist[ifl][0][:int(len(featurelist[0][0])*batchsize*1.2)]
                       featurelist[ifl][1] = featurelist[ifl][1][:int(len(featurelist[0][0])*batchsize*1.2)]
                if len(featurelist[0][0])<1: continue
             #
             if len(filelists)>1:
                for ifl in range(1, len(filelists)):
                   #if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files: msizeprint.append(str(len(featurelist[ifl][0])))
                   featurelist[ifl][0] = np.array_split(featurelist[ifl][0], len(featurelist[0][0]))
                   featurelist[ifl][1] = np.array_split(featurelist[ifl][1], len(featurelist[0][0]))
             io_time += (time.time() - io_start_time)

             ifl=3 if len(featurelist)>3 else len(featurelist)-1
             if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files:
                sess.run(init_l)
                try:
                   # print some testing information as progress indicators
                   loss, aucm, acc, p, r = sess.run([loss_op, auc_op[1], accuracy, mpre[1], mspf[1]], feed_dict={X:featurelist[ifl][0][0], Y:featurelist[ifl][1][0]})
                   print(">>>Tratin#files "+str(file_group_id[0]+1)+",loss="+"{:.3f}".format(loss)+",AUC="+"{:.3f}".format(aucm)+",acc="+"{:.3f}".format(acc)+",p="+"{:.3f}".format(p)+",r="+"{:.3f}".format(r)+(" Comsuming time: %d(current=%d) IO=%d(%.3f)" % (time.time()-start_time, time.time()-start_c_time, io_time, io_time/float(time.time()-start_time))));
                except:
                   print(">>>Tratin#filesError "+str(file_group_id[0]+1)+(" Comsuming time: %d(current=%d) IO=%d(%.3f)" % (time.time()-start_time, time.time()-start_c_time, io_time, io_time/float(time.time()-start_time))));
                sys.stdout.flush()
                start_c_time = time.time();

             # using each subgroup of data for training
             for subi in range(len(featurelist[0][0])):
                for ifl in range(len(filelists)):
                   to = sess.run([train_op, loss_op], feed_dict={X:featurelist[ifl][0][subi], Y:featurelist[ifl][1][subi]})
                   if len(featurelist)==1:
                      # print some detail if nan issue happens
                      if math.isnan(to[1]):
                         for toj in range(len(featurelist[ifl][0][subi])):
                            print('{} vs {}'.format(featurelist[ifl][1][subi][toj][0], featurelist[ifl][1][subi][toj][1]))
                            for tok in featurelist[ifl][0][subi][toj]:
                               opstr = []
                               for tol in tok:
                                   opstr.append(str(round(tol, 2)))
                               print("\t\t\t"+','.join(opstr))
                         sys.exit(1)


             # adjust progress output information
             ifl=3 if len(featurelist)>3 else len(featurelist)-1
             if (file_group_id[0]+1) - last_desplay_files_num >= desplay_files:
                last_desplay_files_num = (file_group_id[0]+1) - ((file_group_id[0]+1) % desplay_files)

             # store more models
             if 49.5<int(file_group_id[0]*100/float(len(filelists[0])))<50.5:
                savp = '.50'
                if (not os.path.isdir(moptions['outFolder']+str(step-1)+savp)):
                  os.system('mkdir -p '+moptions['outFolder']+str(step-1)+savp);
                saver.save(sess, moptions['outFolder']+str(step-1)+savp+'/'+moptions['FileID']);
             if len(featurelist)==1:
                cur_per = int(file_group_id[0]*100/float(len(filelists[0])))
                if cur_per in [10, 20, 30, 40, 60, 70, 80, 90]:
                    savp = str(round(cur_per/100.0, 2))
                    if (not os.path.isdir(moptions['outFolder']+str(step-1)+savp)):
                        os.system('mkdir -p '+moptions['outFolder']+str(step-1)+savp);
                    saver.save(sess, moptions['outFolder']+str(step-1)+savp+'/'+moptions['FileID']);
         # for each epoch, store the trained model
         if (not os.path.isdir(moptions['outFolder']+str(step))):
            os.system('mkdir -p '+moptions['outFolder']+str(step));
         saver.save(sess, moptions['outFolder']+str(step)+'/'+moptions['FileID']);
      print("Training Finished!")

      return (accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred)

#
# get all data files in a folder
#
def getTFiles1(folder1, moptions):
   t1files = glob.glob(os.path.join(folder1, "*.xy.gz"))
   # get all data in a recursive way
   if moptions['recursive']==1:
      t1files.extend(glob.glob(os.path.join(folder1, "*/*.xy.gz")))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*.xy.gz")));
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*.xy.gz")));
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*/*.xy.gz")));
   print("Get folder1");
   # for read-based independent testing
   print(t1files.__sizeof__(), len(t1files))
   if moptions['test'][0] == '0':
      if moptions['test'][1]>0.5:
         t1files = t1files[:int(len(t1files)*moptions['test'][1])]
      else: t1files = t1files[-int(len(t1files)*moptions['test'][1]):]
   print('Sizeinfo: %s sizeof=%d len=%d' % (folder1, t1files.__sizeof__(), len(t1files)))
   sys.stdout.flush();

   return t1files

#
# get all data files in two seperate folders
#
def getTFiles(folder1, folder2, moptions):
   t1files = glob.glob(os.path.join(folder1, "*.xy.gz")); #print(t1files.__sizeof__(), len(t1files))
   # get all data in a recursive way
   if moptions['recursive']==1:
      t1files.extend(glob.glob(os.path.join(folder1, "*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
      t1files.extend(glob.glob(os.path.join(folder1, "*/*/*/*/*.xy.gz"))); #print(t1files.__sizeof__(), len(t1files))
   print("Get folder1");
   # for read-based independent testing
   print(t1files.__sizeof__(), len(t1files))
   if moptions['test'][0] == '0':
      if moptions['test'][1]>0.5:
         t1files = t1files[:int(len(t1files)*moptions['test'][1])]
      else: t1files = t1files[-int(len(t1files)*moptions['test'][1]):]
   print(t1files.__sizeof__(), len(t1files))
   sys.stdout.flush();

   if folder2==None: t2files = []
   else:
      # get all data in a recursive way for the second folder
      t2files = glob.glob(os.path.join(folder2, "*.xy.gz"))
      if moptions['recursive']==1:
         t2files.extend(glob.glob(os.path.join(folder2, "*/*.xy.gz")))
         t2files.extend(glob.glob(os.path.join(folder2, "*/*/*.xy.gz")))
         t2files.extend(glob.glob(os.path.join(folder2, "*/*/*/*.xy.gz")))
         t2files.extend(glob.glob(os.path.join(folder2, "*/*/*/*/*.xy.gz")))
      print("Get folder2");
      # for read-based independent testing
      print(t2files.__sizeof__(), len(t2files))
      if moptions['test'][0] == '0':
         if moptions['test'][1]>0.5:
            t2files = t2files[:int(len(t2files)*moptions['test'][1])]
         else: t2files = t2files[-int(len(t2files)*moptions['test'][1]):]
      print(t2files.__sizeof__(), len(t2files))
      sys.stdout.flush();
   return t1files, t2files

#
# get data from a data file
#
def getDataFromFile(fn, moptions):
   mdata = np.loadtxt(fn, dtype=np.float32)
   # get genomic position, label, and feature
   t0, ty, tx = np.split(mdata, [0,2], axis=1);
   return (tx, ty, None)

#
# get X and Y from a data file
#
def getDataFromFile_new(fn, moptions, mfind0ld=None):
   mdata = np.loadtxt(fn, dtype=np.float32)
   # get genomic position, label, and feature
   t0, ty, tx = np.split(mdata, [1,3], axis=1);

   if moptions['test'][0] in ['-', '+']:
      t0 = t0.astype(int)

   nan_file = []
   m_data = []; m_y = [];
   if not mfind0ld==None:
      pos_to_file_dict = defaultdict();  preind = 0
      mfind0ldkeys = sorted(list(mfind0ld.keys()));
   for mind in range(len(ty)):
      if not mfind0ld==None:
         if preind<len(mfind0ldkeys) and mind == mfind0ldkeys[preind]:
            pos_to_file_dict[len(m_y)] = mfind0ld[ mfind0ldkeys[preind] ]
            preind += 1

      # for region-based independent testing
      if (ty[mind][0]<0.01 and ty[mind][1]<0.01): continue;
      if (moptions['test'][0]=='-' and moptions['test'][1]<t0[mind]<moptions['test'][2]) or \
         (moptions['test'][0]=='+' and (not moptions['test'][1]<t0[mind]<moptions['test'][2])):
         continue;
      ## check nan values in feature files
      has_nan_value = False;
      for cur_row in tx[(mind-int(moptions['windowsize']/2)):(mind+int(moptions['windowsize']/2)+1)]:
         if np.isnan(cur_row).any():
            has_nan_value = True;
            break;
      if has_nan_value:
         if fn in nan_file: pass
         else:
             print ("Warning-nan-value {}".format(fn))
             nan_file.append(fn);
      else:
         m_y.append(ty[mind])
         m_data.append(tx[(mind-int(moptions['windowsize']/2)):(mind+int(moptions['windowsize']/2)+1)])
   if not mfind0ld==None:
      # since each feature file contains data from several fast5 files
      # get the staring rows of fast5 files in feature files
      file_to_pos_dict = defaultdict();
      ptofkeys = sorted(list(pos_to_file_dict.keys()))
      for npk_ind in range(len(ptofkeys)):
         if (npk_ind+1<len(ptofkeys) and ptofkeys[npk_ind+1]-ptofkeys[npk_ind]<500) or len(m_y)-ptofkeys[npk_ind]<500: continue;

         file_to_pos_dict[ pos_to_file_dict[ptofkeys[npk_ind]] ] = [ptofkeys[npk_ind], (ptofkeys[npk_ind+1] if npk_ind+1<len(ptofkeys) else len(m_y))]

   # reshape the data
   if len(m_data)>0:
      m_data = np.reshape(m_data, (len(m_data), len(m_data[0]), len(m_data[0][0])))
      m_y = np.reshape(m_y, (len(m_y), len(m_y[0]))).astype(int)

   if not mfind0ld==None:
      return (m_data, m_y, file_to_pos_dict);
   else: return (m_data, m_y, None)

#
# get index information for each fast5 file in a data file
#
def getGZFilePos(gzfile):
   # get the staring row of fast5 files in a feature file
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

#
# make a prediction for data from a data file
#
def mPred(mfbase, mffolder, accuracy, X, Y, test_gzfile2, pf, num_input, auc_op, mpre, mspf, init_l, mfpred, timesteps, moptions):
   config = tf.ConfigProto()
   config.gpu_options.allow_growth = True
   with tf.Session(config=config) as sess:
      # restore the well-trained model
      new_saver = tf.train.import_meta_graph(mfbase+'.meta')
      new_saver.restore(sess,tf.train.latest_checkpoint(mffolder))

      # write the prediction information
      pfwriter = open(pf, 'w');
      for test_gzfile in test_gzfile2:
         for test_fn_ind in range(len(test_gzfile)):
            # get feature files for prediction
            test_gzfeature, test_gzlabel, _ = getDataFromFile_new(test_gzfile[test_fn_ind], moptions)
            if len(test_gzfeature)<1: continue;

            ftlist = np.array_split(test_gzfeature, int(len(test_gzfeature)/batchsize)+1)
            lblist = np.array_split(test_gzlabel, int(len(test_gzlabel)/batchsize)+1)
            for fti in range(len(ftlist)):
               sess.run(init_l)
               mfpred_output = sess.run([ mfpred], feed_dict={X:ftlist[fti], Y:lblist[fti]})
               # get ture positive, false positive, false negative and true negative
               tp, fp, fn, tn = 0, 0, 0, 0
               for pi in range(len(mfpred_output[0])):
                  if mfpred_output[0][pi]==1:
                     if lblist[fti][pi][1]==1: tp += 1
                     else: fp += 1
                  else:
                     if lblist[fti][pi][1]==1: fn += 1
                     else: tn += 1
               pfwriter.write('tp=%d fp=%d fn=%d tn=%d %s\n' % (tp, fp, fn, tn, test_gzfile[test_fn_ind]))
            pfwriter.flush()
      pfwriter.close();

#
# entry for prediction
#
def pred_prepare(moptions, test_file, accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred):
    mPred(moptions['modfile'][0], moptions['modfile'][1], accuracy, X, Y, test_file, moptions['outFolder']+moptions['FileID']+'_mpred.txt', moptions['fnum'], auc_op, mpre, mspf, init_l, mfpred, moptions['windowsize'], moptions)

#
# prepare training process
#
def mMult_RNN_LSTM_train(moptions):
   # get folders of feature files
   filegroups = moptions['wrkBase'].split(';')
   for i in range(len(filegroups)):
      filegroups[i] = filegroups[i].split(',')

   print(filegroups)

   # get feature files
   filelists = [[] for _ in range(len(filegroups))]
   for i in range(len(filegroups)):
      for fgj in range(len(filegroups[i])):
         if not len(filegroups[i][fgj])>0: continue
         filelists[i].extend(getTFiles1(filegroups[i][fgj], moptions))
   mostnum, mostid = 0, -1;
   # set seed so that random setting
   np.random.seed(3)
   for i in range(len(filelists)):
      np.random.shuffle(filelists[i])
      if len(filelists[i])>mostnum:
         mostnum = len(filelists[i])
         mostid = i;

   # set seed so that random setting
   np.random.seed(7)
   # path for model file
   if 'modfile' in moptions and (not moptions['modfile']==None):
      if moptions['modfile'].rfind('/')==-1:
         moptions['modfile'] = [moptions['modfile'], './']
      else:
         moptions['modfile'] = [moptions['modfile'], moptions['modfile'][:moptions['modfile'].rfind('/')+1]]

   if not mostid==0:
      filelists[mostid], filelists[0] = filelists[0], filelists[mostid]

   accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred = train_save_model(filelists, moptions['fnum'], moptions['hidden'], moptions['windowsize'], moptions)

#
# prepare prediction process
#
def pred_entry(moptions):
   # get feature files
   tfiles = [getTFiles1(moptions['wrkBase'], moptions)]

   init, init_l, loss_op, accuracy, train_op, X, Y, saver, auc_op, mpre, mspf, mfpred = mCreateSession(moptions['fnum'], moptions['hidden'], moptions['windowsize'], moptions)

   # path for model file
   if moptions['modfile'].rfind('/')==-1:
      moptions['modfile'] = [moptions['modfile'], './']
   else:
      moptions['modfile'] = [moptions['modfile'], moptions['modfile'][:moptions['modfile'].rfind('/')+1]]

   pred_prepare(moptions, tfiles, accuracy, X, Y, auc_op, mpre, mspf, init_l, mfpred)
