############################################################
##### 01. bulk RNA-seq raw data(FASTQ) mapping using salmon
##### April 2025, Hyojin Kang
############################################################

#!/bin/sh

############################################################
##### salmon version 1.1.0
### https://combine-lab.github.io/salmon/
SALMON=/path/salmon/salmon-1.1.0/bin/salmon

############################################################
##### salmon index file
GENCODE=/path/db/gencode
SALMON_INDEX=$GENCODE/Gencode_mouse/salmon_index
GTF=$GENCODE/Gencode_mouse/gencode.vM23.annotation.gtf

############################################################
##### mouse reference genome file
### https://support.illumina.com/sequencing/sequencing_software/igenome.html
VERSION=/path/iGenomes/Mus_musculus/Ensembl/GRCm38
GENOME_FASTA=$VERSION/Sequence/WholeGenomeFasta/genome.fa

############################################################
##### salmon run script
FASTQ1=$WORKDIR'/fastq/'$TAG'_1.fastq.gz'
FASTQ2=$WORKDIR'/fastq/'$TAG'_2.fastq.gz'
SALMON_OUT=$WORKDIR'/salmon/'$TAG

$SALMON quant -i $SALMON_INDEX --libType ISR -p 8 \
 --seqBias --validateMappings --mimicBT2 \
 --gcBias -1 $FASTQ1 -2 $FASTQ2 -o $SALMON_OUT -g $GTF

############################################################
### --libType ISR: Strand-specific RNA-seq using dUTP method
# I(Inward): Paired-end reads are aligned facing each other
# S(Stranded): strand-specific library
# R: Second read is in the same direction as the transcript (first read is in the opposite direction)
