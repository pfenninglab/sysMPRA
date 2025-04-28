# sysMPRA
This repository is for **sys**temic **MPRA** code.

## Scripts
designMPRA1Sequences.sh: code for designing sequences related to MEF2C binding for the gMPRA

getEnsemblCoordinatesFromGff3.py: code for obtaining gene coordinates from gff3 files

convertFIMOToMotifHitBed.py: code for converting motif hits from FIMO from the MEME suite to a bed file

chooseBestScramble.py: code for choosing the shuffled sequence with the least significant motif match

filterFastaFile.py: filter a fasta file based on the fasta headers

predictMPRASequences.sh: make open chromatin predictions for gMPRA sequences

arrayProc: directory with scripts for processing MPRA data

## Dependencies
bedtools version 2.26.0

meme suite version 4.12.0

liftOver (UCSC Genome Browser utility)

python version 2.7.17 (chooseBestScramble.py also works with python version 3.7.6)

numpy version 1.13.3

biopython version 1.70

twoBitToFa (UCSC Genome Browser utility)

OCROrthologPrediction (https://github.com/pfenninglab/OCROrthologPrediction)

## Contacts
Andreas Pfenning (apfenning@cmu.edu)

Irene Kaplow (ikaplow@cs.cmu.edu)

Ashley Brown (arb1@andrew.cmu.edu)

## Reference
Brown AR, Fox GA*, Kaplow IM*, Lawler AJ, Phan BN, Gadey L, Wirthlin ME, Ramamurthy E, May GE, Chen Z, Su Q, McManus CJ, van der Weerd R, Pfenning AR.  An in vivo systemic massively parallel platform for deciphering animal tissue-specific regulatory function.  Frontiers in Genetics, 16:1533900, 2025.  https://pubmed.ncbi.nlm.nih.gov/40270544/.
