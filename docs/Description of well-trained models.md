This file contains description of well-trained model in `train_mod` directory.  ***This file is under construction now. Will be done soon.***

# Modificaiton prediction model (LSTM(Long short-term memory)-based RNN (Recurrent neural network))
## `rnn_sinmodC_P100wd21_f7ne1u0_4`
This model is trained on E. Coli data with synthetically introduced 5mC.
```
Window size: 21
Feature/event: 7
Base of interest: C
Training epoch of negative data: 4
Training data: SSS and UMR
```

## `rnn_conmodC_P100wd21_f7ne1u0_4`
This model is trained on E. Coli data with synthetically introduced 5mC.
```
Window size: 21
Feature/event: 7
Base of interest: C
Training epoch of negative data: 4
Training data: positive control (SSS, Cg_sssl, Cg_mpel and gCgC), and negative control: (UMR, con1 and con2)
```

## `rnn_conmodA_P100wd21_f7ne1u0_4`
This model is trained on E. Coli data with synthetically introduced 6mA.
```
Window size: 21
Feature/event: 7
Base of interest: A
Training epoch of negative data: 4
Training data: Three positive control of 6mA and the negative control (con1, con2)
```

## `rnn_conmodA_E1m2wd21_f7ne1u0_4`
This is a region-based training model. Reads and bases mapped to 1,000,000 on E. Coli is used for testing and others for training.
```
Window size: 21
Feature/event: 7
Base of interest: A
Training epoch of negative data: 4
Training data: Three positive control of 6mA and the negative control (con1, con2)
```

## `rnn_f7_wd21_chr1to10_4`
This model is trained on Chr 1 to 10 of NA12878 with completely methylated positions (>90% methylation percentage in both replicates of bisulfite sequences) and completely un-methylated positions (<=0% methylation percentage in both replicates of bisulfite sequences). Bases in long reads mapped to the two groups of reference positions are used in training process.
```
Window size: 21
Feature/event: 7
Base of interest: C
Training epoch of all data: 4
Training data: Chr 1 to 10 of NA12878
```

# Cluster-effect model (the second neural network)
## `na12878_cluster_train_mod-keep_prob0.7-nb25-chr1`
```
Window size: 25
Feature/event: 11
Base of interest: C
Training epoch of negative data: 1
Training data: Chr 1 of NA12878
```
