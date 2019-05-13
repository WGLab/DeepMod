This document describe the fomrat of the output of DeepMod.

# 1. Check whether the running of DeepMod is successful.
After finishing the running of DeepMod, usually you will find a "\*.done" file under the output folder. If you cannot see that file, some errors might happen and the running of DeepMod is not successful. One common cause behind this is the memory issue, that is, DeepMod needs more memory than what you requested or have. Increasing memory generally solve the issue.

# 2 Format of output
The output is in a BED format like below. The first six columns are `Chr`, `Start pos`, `End pos`, `Base`, `Capped coverage`, and `Strand`, and the last three columns are `Real coverage`, `Mehylation percentage` and `Methylation coverage`.

```
chr6 148655 148656 C 10 -  148655 148656 0,0,0 10 10 1
chr6 148657 148658 C 12 +  148657 148658 0,0,0 12 8 1
chr6 148674 148675 C 14 -  148674 148675 0,0,0 14 7 1
chr6 148675 148676 C 15 -  148675 148676 0,0,0 15 6 1
chr6 148676 148677 C 14 -  148676 148677 0,0,0 14 7 1
chr6 148684 148685 C 12 -  148684 148685 0,0,0 12 25 3
chr6 148685 148686 C 16 -  148685 148686 0,0,0 16 6 1
chr6 148689 148690 C 11 +  148689 148690 0,0,0 11 72 8
chr6 148691 148692 C 10 +  148691 148692 0,0,0 10 50 5
chr6 148693 148694 C 8 +  148693 148694 0,0,0 8 100 8
chr6 148694 148695 C 11 -  148694 148695 0,0,0 11 54 6
chr6 148695 148696 C 10 +  148695 148696 0,0,0 10 90 9
chr6 148697 148698 C 12 +  148697 148698 0,0,0 12 50 6
chr6 148699 148700 C 9 +  148699 148700 0,0,0 9 22 2
chr6 148701 148702 C 13 -  148701 148702 0,0,0 13 7 1
chr6 148703 148704 C 13 -  148703 148704 0,0,0 13 15 2
chr6 148706 148707 C 9 -  148706 148707 0,0,0 9 22 2
```
