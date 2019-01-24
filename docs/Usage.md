This is an explanation of how to run DeepMod without examples. If you want to run some example, please refer to [demo](https://github.com/WGLab/DeepMod/blob/master/docs/Reproducibility.md).


The inputs of DeepMod is a group of FAST5 files and a reference genome. FAST5 files need to be basecalled already, and `Events` data must be availabe in FAST5 files. 

DeepMod has a functional module called "detect" which will detect a specific modification in a single run. However, if the dataset and genome size is very larger or one wants to have the results soon, it would be better to run "detect" in a separate process simultaneously and then merge them together. For some special cases, if there is cluster effect between modifications (such as 5mC in CpG cluster), the third process would be used for additional prediction after "detect". How to use the three functional modules is described below.

# 1. How to detect modifications from FAST5 files.
The command for modification detection is to run `python DeepMod.py detect`. Without any other parameters, the help document will be shown. An example of how to use it is given below.

```
python DeepMod.py detect --wrkBase FAST5-Folder --Ref Ref_genome_path --outFolder out_folder --Base C --modfile train_mod/rnn_f7_wd21_chr1to10_4/mod_train_f7_wd21_chr1to10 --FileID User_Uniq_name --threads 4
```
where users need to provide where is the FAST5 files (`--wrkBase`), where is the reference genome (`--Ref`), where is the output folder (`--outFolder`), and also the base of interest and the mod file. Users can optionally specify unique string for the results file names (`--FileID`) and how many threads are used (`--threads`).

If you want to make the prediction for base `A`, the following command could be used.
```
python DeepMod.py detect --wrkBase FAST5-Folder --Ref Ref_genome_path --outFolder out_folder --Base A --modfile train_mod/rnn_conmodA_P100wd21_f7ne1u0_4/mod_train_conmodA_P100wd21_f3ne1u0 --FileID User_Uniq_name --threads 4
```


# 2. How to merge different runs of modification detection


# 3. How to consider modification cluster effect.



