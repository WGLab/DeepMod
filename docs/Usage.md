The inputs of DeepMod is a group of FAST5 files and a reference genome. FAST5 files need to be basecalled already, and `Events` data must be availabe in FAST5 files. 

DeepMod has a functional module called "detect" which will detect a specific modification in a single run. However, if the dataset and genome size is very larger or one wants to have the results soon, it would be better to run "detect" in a separate process simultaneously and then merge them together. For some special cases, if there is cluster effect between modifications (such as 5mC in CpG cluster), the third process would be used for additional prediction after "detect". How to use the three functional modules is described below.

# 1. How to detect modifications from FAST5 files.


# 2. How to merge different run of modification detection


# 3. How to consider modification cluster effect.



