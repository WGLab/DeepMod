# NanoDeepMod: a deep-learning tool for genomic-scale, strand-sensitive and single-nucleotide based detection of DNA modifications

## Methodology of NanoDeepMod

NanoDeepMod is a computational tool which takes long-read signals as input and outputs modification summary for each genomic position in a reference genome together with modification prediction for each base in a long read. The modification prediction model in NanoDeepMod is a well-trained bidirectional recurrent neural network (RNN) with long short-term memory (LSTM) units. LSTM RNN is a class of artificial neural network for modeling sequential behaviors with LSTM to preclude vanishing gradient problem.  To detect DNA modifications, normalized signals of events in a long read were rescaled from -5 and 5, and discretized into 50 bins, and the count in 50 bins of normalized signals together with signal mean, standard deviation and the number of signals (denoted by 53-feature description) were obtained for each event as input of a LSTM unit with 100 hidden nodes. In NanoDeepMod with 3 hidden layers in RNN. Predicted modification summary for each position would be generated in a BED format, suggesting how many reads cover genomic positions, how many mapped bases in long reads were predicted to be modified and the coverage percentage of prediction modifications. This modification prediction by NanoDeepMod is strand-sensitive and single-nucleotide based. 

## Inputs of NanoDeepMod

The input of NanoDeepMod is Nanopore long read data together a refrence genome. 

## Usage

Please refer to [Usage](https://github.com/WGLab/NanoDeepMod/blob/master/docs/Usage.md) for how to use NanoDeepMod.

## Revision History

For release history, please visit [here](https://github.com/WGLab/NanoDeepMod/releases). For details, please go [here](https://github.com/WGLab/NanoDeepMod/blob/master/README.md).

## Contact

If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/WGLab/NanoDeepMod/issues). They would also be helpful to other users. 

## Reference

