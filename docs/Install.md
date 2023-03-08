
# Installation Guide

## Prerequisites:
 ### The required packages for running DeepMod are listed below:
	* Python 3.6
	* python packages:
		+ h5py
		+ numpy
		+ scipy
	* BWA MEM/minimap2
	* SAMTOOLS
	* HDF5
	* Tensorflow
	* If running performance evaluation, you might need
		+ python packages: rpy2, matplotlib, scikit-learn
		+ R packages: ggplot2, gridExtra, plyr, RColorBrewer

The packages I used are 
```
h5py                      2.7.1
hdf5                      1.10.1
numpy                     1.14.0
scikit-learn              0.19.1
scipy                     1.0.0
tensorflow                1.7.0
samtools                  1.5
minimap2                  2.12
bwa                       0.7.15
matplotlib                2.1.2
rpy2                      2.8.6
r                         3.4.2
```

  ### Package installation		
Users have their own choice of how to install required package above. But to efficiently manage the installation packages, Anaconda is recommended. After installing Annocoda, it would also be an good option to use virtual environment in annocoda. `conda create -n mdeepmod python=3.6` can be used to create a virtual environment, and `source activate mdeepmod` can be used to activate the virtual environment, and then install required packages or run DeepMod. If users want to exit the virtual environment, simply type `source deactivate`. 

After creating `mdeepmod` virtual environment using annocoda, the following commands can install majority of necessary packages:

```
source activate mdeepmod
conda install -c anaconda h5py hdf5 numpy scipy scikit-learn
conda install -c bioconda minimap2 bwa samtools
conda install -c r rpy2 r-ggplot2 r-gridextra r-plyr r-rcolorbrewer
```
Users might need to install `tensorflow` with version 1.13 by yourself or by using `conda install -c conda-forge tensorflow matplotlib` (if there are conflicts when using this command, uers need to give up and install `tensorflow` and `matplotlib` by yourself.).

### Additional notes

Some users reported that the compression format of their fast5 is vbz instead of gzip. Install `ont-vbz-hdf-plugin` solves the problem. Otherwise, an error message such as "Error!!! No Raw_reads/Signal data" will be shown.


## Install DeepMod:
	* git clone https://github.com/WGLab/DeepMod
	* The script to be run is in bin/DeepMod.py: 
		+ Run `python bin/DeepMod.py` for help document or 
		+ Run `python bin/DeepMod.py detect` for the detection help document.
Users might need to run `source activate mdeepmod` if the virtual environment and required packages are installed with the commands above.

## Installation time
Without GPU-version tensorflow, it would take ~30 minutes to install required packages and DeepMod. 

## Usage:
 For how to use them, please refer to [Usage](https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md)

