# DeepMod: a deep-learning tool for genomic-scale, strand-sensitive and single-nucleotide based detection of DNA modifications

## Methodology of DeepMod

DeepMod is a computational tool which takes long-read signals as input and outputs modification summary for each genomic position in a reference genome together with modification prediction for each base in a long read. The modification prediction model in DeepMod is a well-trained bidirectional recurrent neural network (RNN) with long short-term memory (LSTM) units. LSTM RNN is a class of artificial neural network for modeling sequential behaviors with LSTM to preclude vanishing gradient problem.  To detect DNA modifications, normalized signals of events in a long read were rescaled from -5 and 5, and signal mean, standard deviation and the number of signals together with base information (denoted by 7-feature description) were obtained for each event as input of a LSTM unit with 100 hidden nodes. In DeepMod with 3 hidden layers in RNN. Predicted modification summary for each position would be generated in a BED format, suggesting how many reads cover genomic positions, how many mapped bases in long reads were predicted to be modified and the coverage percentage of prediction modifications. This modification prediction by DeepMod is strand-sensitive and single-nucleotide based. 

### Inputs of DeepMod

The input of DeepMod is Nanopore long read data together a refrence genome. 

## System Requirements
### Hardware requirements
DeepMod is based on deep learning framework, and needs to access raw data of Nanopore sequencing. Thus, it needs enough RAM to support deep learning framework and enough hard drive for raw data of Nanopore sequencing. GPU can substantially speedup the detection process. For optimal performance, we recommend a computer with:
 * RAM: 20+ GB per thread
 * GPU or CPU with 8+ cores
 * HDD or better with SSD. Dependent on how large raw data is (for 30X E coli data, it might need 10+GB, while for 30X human data, it might need 10+TB)
 
### Software requirements
The developmental version of DeepMod has been tested on Linux operating system: CentOS 7.0 with both CPU and GPU machines.

## Installation
Please refer to [Installation](https://github.com/WGLab/DeepMod/blob/master/docs/Install.md) for how to install DeepMod.

## Usage

Please refer to [Usage](https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md) for how to use DeepMod.

## Examples and Reproducibility of our analysis.

Please refer to [Examples and Reproducibility](https://github.com/WGLab/DeepMod/blob/master/docs/Reproducibility.md) for examples of how to run DeepMod.

## Revision History

For release history, please visit [here](https://github.com/WGLab/NanoDeepMod/releases). For details, please go [here](https://github.com/WGLab/DeepMod/blob/master/README.md).

## Contact

If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/WGLab/DeepMod/issues). They would also be helpful to other users. 

## Reference

