```
Usage: python DeepMod.py [-h] [--bam BAM] [-chrom [CHROM [CHROM ...]]] [--ref REF] --fast5 FAST5 --output OUTPUT [--threads THREADS] [--file_name FILE_NAME]
                  [--model MODEL] [--guppy_group GUPPY_GROUP] [--tombo_group TOMBO_GROUP] [--basecaller {guppy,tombo}] [-wgs_contigs_type WGS_CONTIGS_TYPE]

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             Path to bam file if Guppy basecaller is user. BAM file is not needed with Tombo fast5 files.
  -chrom [CHROM [CHROM ...]], --chrom [CHROM [CHROM ...]]
                        A space/whitespace separated list of contigs, e.g. chr3 chr6 chr22. Only applicable with guppy basecaller.
  --ref REF             Path to reference file
  --fast5 FAST5         Path to folder containing tombo requiggle Fast5 files. Fast5 files will be recusrviely searched
  --output OUTPUT       Path to folder where features will be stored
  --threads THREADS     Number of processors to use
  --file_name FILE_NAME
                        Name of the output file
  --model MODEL         Name of the model. Current options are "guppy_na12878" and "tombo_na12878". Use according to the basecaller specified.
  --guppy_group GUPPY_GROUP
                        Name of the guppy basecall group
  --tombo_group TOMBO_GROUP
                        Name of the tombo group
  --basecaller {guppy,tombo}
                        Use Tombo or Guppy output
  -wgs_contigs_type WGS_CONTIGS_TYPE, --wgs_contigs_type WGS_CONTIGS_TYPE
                        Options are "with_chr", "without_chr" and "all", "with_chr" option will assume human genome and run DeepMod on chr1-22 X Y,
                        "without_chr" will run on chromosomes 1-22 X Y if the BAM and reference genome files use chromosome names without "chr". "all" option
                        will run DeepMod on each contig present in reference genome FASTA file. Only applicable with guppy basecaller.
```
