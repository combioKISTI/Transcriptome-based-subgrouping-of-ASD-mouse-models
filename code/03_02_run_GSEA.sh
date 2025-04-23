############################################################
##### 03.02 Gene Set Enrichment Analysis (GSEA)
##### April 2025, Hyojin Kang

############################################################
### Gene Set Enrichment Analysis (GSEA) 
### homepage: http://software.broadinstitute.org/gsea

############################################################
### GSEA was used to capture if the expressions of genes in a specific gene set are changed in a consistent direction, 
### although each might not be significant enough to be counted as differentially expressed genes. 
### Enrichment analysis was performed using GSEAPreranked (gsea-4.3.3.jar) module on gene set collections downloaded 
### from Molecular Signature Database (MSigDB) v2023.2.Hs (https://www.gsea-msigdb.org/gsea/. 
### GSEAPreranked was applied using the list of all genes expressed, ranked by the fold change and multiplied 
### by the inverse of the P value with recommended default settings (1,000 permutations and a classic scoring scheme). 
### The gene sets with an FDR of less than 0.05 were considered as significantly enriched. 


#!/bin/bash

############################################################
### GSEA version
GSEA_HOME='/path/GSEA'
GSEA_RUN=$GSEA_HOME'/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked'

### MsigDB version
VERSION='v2023.2.Hs'
HUMAN_CHIP_DIR='msigdb_'$VERSION'_chip_files_to_download_locally/'
HUMAN_GMT_DIR='msigdb_'$VERSION'_GMTs/'

############################################################
### input rnk file
RNK=$WORKDIR'/'$GENE'_deseq2_'$COMPTAG'_all_human.txt.rnk'
GMT=$GSEA_HOME'/'$HUMAN_GMT_DIR'c5.go.bp.'$VERSION'.symbols.gmt'
HUMAN_CHIP='Human_HGNC_ID_MSigDB.'$VERSION'.chip'

### output file
OUTPUT=$WORKDIR'/gsea/'$SAMPLE

############################################################
### GSEA analysis with c5.go.bp category

LABEL=$SAMPLE'_c5_bp'
$GSEA_RUN \
 -gmx $GMT -rnk $RNK -chip $HUMAN_CHIP -collapse false \
 -mode Max_probe -norm meandiv -nperm 1000 \
 -scoring_scheme classic -rpt_label $LABEL \
 -include_only_symbols true -make_sets true -plot_top_x 30 \
 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false \
 -out $OUTPUT

############################################################
### GSEA analysis with c5.go.cc category
GMT=$GSEA_HOME'/'$HUMAN_GMT_DIR'c5.go.cc.'$VERSION'.symbols.gmt'
LABEL=$SAMPLE'_c5_cc'
$GSEA_RUN \
 -gmx $GMT -rnk $RNK -chip $HUMAN_CHIP -collapse false \
 -mode Max_probe -norm meandiv -nperm 1000 \
 -scoring_scheme classic -rpt_label $LABEL \
 -include_only_symbols true -make_sets true -plot_top_x 30 \
 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false \
 -out $OUTPUT

############################################################
### GSEA analysis with c5.go.mf category
GMT=$GSEA_HOME'/'$HUMAN_GMT_DIR'c5.go.mf.'$VERSION'.symbols.gmt'
LABEL=$SAMPLE'_c5_mf'
$GSEA_RUN \
 -gmx $GMT -rnk $RNK -chip $HUMAN_CHIP -collapse false \
 -mode Max_probe -norm meandiv -nperm 1000 \
 -scoring_scheme classic -rpt_label $LABEL \
 -include_only_symbols true -make_sets true -plot_top_x 30 \
 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false \
 -out $OUTPUT
