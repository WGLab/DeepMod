This file contains description of well-trained model in `train_mod` directory. The meaning of the short name of the data set can be found in the paper. ***Warnings: the well-trained models below are NEVER retrained with Guppy. You will get unpredictable results if you use them to make the prediction with the data basecalled with Guppy.***

# 1. Modificaiton prediction model
These are several trained models of LSTM(Long short-term memory)-based RNN (Recurrent neural network), and their setting is given below.

## 1.1 `rnn_sinmodC_P100wd21_f7ne1u0_4`
This model is trained on E. Coli data with synthetically introduced 5mC.
```
Window size: 21
#Feature per event: 7
Base of interest: C
#Training epoch of negative data: 4
Training data: SSS and UMR
```
Please refer to our paper for the meanings of `SSS` and `UMR`.

## 1.2 `rnn_conmodC_P100wd21_f7ne1u0_4`
This model is trained on E. Coli data with synthetically introduced 5mC.
```
Window size: 21
#Feature per event: 7
Base of interest: C
#Training epoch of negative data: 4
Training data: positive control (SSS, Cg_sssl, Cg_mpel and gCgC), and negative control: (UMR, con1 and con2)
```
Please refer to our paper for the meanings of `SSS`, `UMR`, `con1`, `con2`, `Cg_sssl`, `Cg_mpel` and `gCgC`.

## 1.3 `rnn_conmodA_P100wd21_f7ne1u0_4`
This model is trained on E. Coli data with synthetically introduced 6mA.
```
Window size: 21
#Feature per event: 7
Base of interest: A
#Training epoch of negative data: 4
Training data: Three positive control of 6mA and the negative control (con1, con2)
```

## 1.4 `rnn_conmodA_E1m2wd21_f7ne1u0_4`
This is a region-based training model. Reads and bases mapped to 1,000,000 on E. Coli is used for testing and others for training.
```
Window size: 21
#Feature per event: 7
Base of interest: A
#Training epoch of negative data: 4
Training data: Three positive control of 6mA and the negative control (con1, con2)
```

## 1.5 `rnn_f7_wd21_chr1to10_4`
This model is trained on Chr 1 to 10 of NA12878 with completely methylated positions (>90% methylation percentage in both replicates of bisulfite sequences) and completely un-methylated positions (<=0% methylation percentage in both replicates of bisulfite sequences). Bases in long reads mapped to the two groups of reference positions are used in training process. This model is used to make 5mC prediction on HX1, and both HX1 and NA12878 are basecalled using Albacore v2.3.1.
```
Window size: 21
#Feature per event: 7
Base of interest: C
#Training epoch of all data: 4
Training data: Chr 1 to 10 of NA12878
```

# 2. Cluster-effect model (the second neural network)
This is the model of the second neural network to consider modificatoin cluster of 5mC. 
## 2.1 `na12878_cluster_train_mod-keep_prob0.7-nb25-chr1`
The model is only trained on Chr1 of NA12878.
```
Window size: 25
#Feature per event: 11
Base of interest: C
#Training epoch of negative data: 1
Training data: Chr 1 of NA12878
```
