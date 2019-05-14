This document describes how to use module in DeepMod to build a well-trained model which can be used in later prediction. Now, 2-class classification is supported. 

There are two steps to build a model: one is to get features from FAST5 files and user-defined modification labels, and the other step is to train a model.

# 1. Feature extraction
There are two ways to extract features for next training process. One is motif-based, and the other is position-based. For large genome, position-based process is recommended.

## 1.1 Motif-based process
For control data, the command below is used for motif-based modification of CpG sites (C is modified) and 10 threads are used.
```python bin/DeepMod getfeatures --posneg 0 --motif CG --ModinMotif 0 --recursive 1 --threads 10 --Ref m_ref_file --wrkBase fast5-folder_control --FileID m_unique_id-control --outFolder m_output/ctrl_features/```

For modification data, the command is
```python bin/DeepMod getfeatures --posneg 1 --motif CG --ModinMotif 0 --recursive 1 --threads 10 --Ref m_ref_file --wrkBase fast5-folder_modification --FileID m_unique_id-mod --outFolder m_output/mod_features/```

## 1.2 Position-based process
For control data, the command below is used for position-based modification and 10 threads are used.
```python bin/DeepMod getfeatures --posneg 0 --fulmod file-complete-modified-position --anymod file-partil-modified-position --nomod file-complete-unmodified-position --recursive 1 --threads 10 --Ref m_ref_file --wrkBase fast5-folder_control --FileID m_unique_id-control --outFolder m_output/ctrl_features/```

For modification data, the command is
```python bin/DeepMod getfeatures --posneg 1 --fulmod file-complete-modified-position --anymod file-partil-modified-position --nomod file-complete-unmodified-position --recursive 1 --threads 10 --Ref m_ref_file --wrkBase fast5-folder_modification --FileID m_unique_id-mod --outFolder m_output/mod_features/```

***For large genome, it is recommended that modified position and un-modified position are generated first and then feed into position-based process. If this is the case, `--anymod` is not necessary, and `--nomod` is necessary and must be given an empty file.***



# 2. Model training
`python bin/DeepMod train `

