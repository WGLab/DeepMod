
# Installation Guide
## Install DeepMod:

First, install Miniconda, a minimal installation of Anaconda, which is much smaller and has a faster installation:

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Go through all the prompts (installation in `$HOME` is recommended). After Anaconda is installed successfully, simply run:

```
git clone https://github.com/WGLab/DeepMod
cd DeepMod
conda env create -f environment.yml
conda activate deepmod
```

Afterwards you can run `python <PATH_TO_REPOSITORY>/DeepMod.py --help` to see usage help.

## Installation time
Without GPU-version tensorflow, it would take ~30 minutes to install required packages and DeepMod. 

## Usage:
 For how to use them, please refer to [Usage](https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md)

