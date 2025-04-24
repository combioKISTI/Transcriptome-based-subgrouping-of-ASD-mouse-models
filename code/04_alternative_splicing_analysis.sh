################################################################################
##### 04. Alternative splicing analysis
##### April 2025, Yukyung Jun
################################################################################

#!/bin/bash

################################################################################
##### rMATS: Replicate Multivariate Analysis of Transcript Splicing
### https://rnaseq-mats.sourceforge.net/
### Version: 4.0.2

################################################################################
##### Define base working directory
wkdir="/your/project/path"

################################################################################
##### Define tool and reference paths
RMATS="$wkdir/tools/rMATS/rmats.py"
GENOME_FASTA="$wkdir/ref/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
GTF="$wkdir/ref/GRCm38/gtf/genes.gtf"
STAR_INDEX="$wkdir/ref/GRCm38/STAR"

################################################################################
##### Experimental setup
MODEL="Adnp" # model name (e.g., Adnp, Chd8-N, etc.)
SEX="M"      # M = Male, F = Female

##### Input FASTQ list (comma-separated or file with full paths)
MUT_FASTQ="$wkdir/fastq_list/${MODEL}_${SEX}_Mut.txt"    # Mutant group
WT_FASTQ="$wkdir/fastq_list/${MODEL}_${SEX}_WT.txt"     # Wildtype group

##### Output directories
OUTDIR="$wkdir/splicing/$MODEL/${MODEL}_${SEX}"
TMPDIR="$wkdir/splicing/tmp/${MODEL}_${SEX}"

################################################################################
##### Run rMATS
python $RMATS  --nthread 8 --readLength 101 --bi $STAR_INDEX \
  --s1 $MUT_FASTQ --s2 $WT_FASTQ --gtf $GTF --od $OUTPUT_DIR --tmp $TMP_DIR

