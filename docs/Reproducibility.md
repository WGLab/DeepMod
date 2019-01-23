# Demo
This section contains examples and also reproducibility process of our analysis. It is assumed that DeepMod has been successfully installed. If not, please install it and its required packages first according to the [instruction ]https://github.com/WGLab/DeepMod/blob/master/docs/Install.md). It is also assumed that the virtual environment `mdeepmod` is used. Please remove `source activate mdeepmod` if virtual environment is not used in your system. 

To prepare to run DeepMod, It is assumed that there is a directory which has a *DeepMod* directory for our DeepMod tool, a *data* directory for storing Nanopore sequencing data, and a *ref* directory for reference data.

## Reference genomes
Please used `bwa` to index the genome before runnning `DeepMod`
### E. coli reference genome
The E. coli reference fasta could be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/556503834. We assumed the reference file name under the *ref* directory is 'Ecoli_k12_mg1655.fasta' and indexed by `bwa`

### Human reference genome
Hg38 was used for human Nanopore sequencing data. We assumed the fasta human file name under the *ref* directory is 'hg38' and indexed by `bwa` 

## Run on E. coli data
The Nanopore data set are large, and thus users need to contact the original authors<sup>1</sup> to get the downloading URL. 
### Example 1
#### Step 1. datasets
Please download Nanopore sequencing data for *con1*, *con2* and *CG* motif with SSSl, and untar them into seperate sub-folder under the *data* directory, and assumed their sub-folder names are *Control_lib1*， *Control_lib3* and *meth10_lib3* respectively.

#### Step 2.
```
mkdir ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/meth10_lib3/ --Ref ref/Ecoli_k12_mg1655.fasta --FileID cgsss --modfile DeepMod/train_mod/rnn_sinmodC_P100wd21_f7ne1u0_4/mod_train_sinmodC_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/Control_lib1/ --Ref ref/Ecoli_k12_mg1655.fasta --FileID con1 --modfile DeepMod/train_mod/rnn_sinmodC_P100wd21_f7ne1u0_4/mod_train_sinmodC_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/Control_lib3/ --Ref ref/Ecoli_k12_mg1655.fasta --FileID con2 --modfile DeepMod/train_mod/rnn_sinmodC_P100wd21_f7ne1u0_4/mod_train_sinmodC_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
```
The three commands might take ~150 minutes, 350 minutes and 500 minutes to be done. After that, you will found the results under the directory of *ecoli_pred* with *--FileID* as the sub-folder name. The results are in the *bed* format with the file names as *mod_pos.NC_000913.3-.C.bed* and *mod_pos.NC_000913.3+.C.bed*. The detail of the methylation prediction for each long read is also provided. 

One can run DeepMod on other datasets with the commands below.
```
time python DeepMod/bin/DeepMod.py detect --wrkBase data/meth8_lib2/ --Ref ref/Ecoli_k12_mg1655.fasta --FileID gCgc --modfile DeepMod/train_mod/rnn_sinmodC_P100wd21_f7ne1u0_4/mod_train_sinmodC_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/meth9_lib2/ --Ref ref/Ecoli_k12_mg1655.fasta --FileID Cgmpe --modfile DeepMod/train_mod/rnn_sinmodC_P100wd21_f7ne1u0_4/mod_train_sinmodC_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
```

### Example 2
#### Step 1. datasets
You need contact the original authors<sup>1</sup> to get the downloading URL.

Please download Nanopore sequencing data for *con1*, *con2* and *gAtc* motif, *tcgA* and *gaAttc*, and untar them into seperate sub-folder under the *data* directory, and assumed their sub-folder names are *Control_lib1*， *Control_lib3*, *meth11_lib3*, *meth1_lib1* and *meth4_lib1* respectively. If you have run Example 1, you might already have the dataset for *con1* and *con2*.

#### Step 2.
```
mkdir ecoli_pred/ (if not exist)
time python DeepMod/bin/DeepMod.py detect --wrkBase data/meth11_lib3/ --Ref ref/Ecoli_k12_mg1655.fasta --Base A --FileID gatc --modfile DeepMod/train_mod/rnn_conmodA_P100wd21_f7ne1u0_4/mod_train_conmodA_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/meth1_lib1/ --Ref ref/Ecoli_k12_mg1655.fasta --Base A --FileID tcga --modfile DeepMod/train_mod/rnn_conmodA_P100wd21_f7ne1u0_4/mod_train_conmodA_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/meth4_lib1/ --Ref ref/Ecoli_k12_mg1655.fasta --Base A --FileID gaattc --modfile DeepMod/train_mod/rnn_conmodA_P100wd21_f7ne1u0_4/mod_train_conmodA_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/Control_lib1/ --Ref ref/Ecoli_k12_mg1655.fasta --Base A --FileID con1a --modfile DeepMod/train_mod/rnn_conmodA_P100wd21_f7ne1u0_4/mod_train_conmodA_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
time python DeepMod/bin/DeepMod.py detect --wrkBase data/Control_lib3/ --Ref ref/Ecoli_k12_mg1655.fasta --Base A --FileID con2a --modfile DeepMod/train_mod/rnn_conmodA_P100wd21_f7ne1u0_4/mod_train_conmodA_P100wd21_f3ne1u0 --threads 15 --outFolder ecoli_pred/
```


## Reference
1. Stoiber MH, et al. De novo Identification of DNA Modifications Enabled by Genome-Guided Nanopore Signal Processing. bioRxiv 10.1101/094672,  (2017).
