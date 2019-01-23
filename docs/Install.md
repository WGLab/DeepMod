
## Prerequisites:
The following packages are needed to run DeepMod. To efficiently manage the installation packages, Anaconda is recommended but users feel free to have their own choices. If Annocoda is used, it is also an good option to use virtual environment in annocoda (such as typing `conda create -n mdeepmod python=2.7 anaconda` to create a virtual environment, and typing `source activate mdeepmod` to activate the virtual environment. After activating `mdeepmod`, and you can install packages and then use DeepMod). `source deactivate` can be used if you want to exit the virtual environment.

 ### The required packages are listed below:
	* Python 2.7
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
If you use annoconda, the following three commands can install majority of necessary packages (`conda install -c anaconda h5py hdf5 numpy scipy scikit-learn`, `conda install -c bioconda minimap2 bwa samtools`, `conda install -c r rpy2 r-ggplot2 r-gridextra r-plyr r-rcolorbrewer`). But you might need to install `tensorflow` by yourself.

## Install:
	* git clone https://github.com/WGLab/DeepMod
	* The script to be run is in bin/DeepMod.py: `python bin/DeepMod.py`
	
## Usage:
 For how to use them, please refer to [Usage](https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md)

