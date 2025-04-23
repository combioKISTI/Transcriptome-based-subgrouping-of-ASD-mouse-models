############################################################
##### 03.01 rank file (RNK) generation 
##### April 2025, Hyojin Kang
############################################################

#!/bin/sh
############################################################
### rnkgen.sh converts a differential gene expression spreadsheet (XLS) into a rank file (RNK)
### generate a ranked list of genes from most up-expressed to most down-expressed based on the p-value
### reference: https://genomespot.blogspot.com/2015/01/how-to-generate-rank-file-from-gene.html

############################################################
### Specify the input file
XLS=$1
### Specify the gene ID column
ID=$2
### Specify the fold change value column
FC=$3
### Specify the raw p-value column
P=$4

############################################################
sed 1d $XLS | tr -d '"' \
| awk -v I=$ID -v F=$FC -v P=$P '{FS="\t"} $I!="" {print $I, $F, $P}' \
| awk '$1!="NA" && $2!="NA" && $3!="NA"' \
| awk '{s=1} $2<0{s=-1} {print $1"\t"s*-1*log($3)/log(10)}' \
| awk '{print $1,$2,$2*$2}' | sort -k3gr \
| awk '{OFS="\t"} !arr[$1]++ {print $1,$2}' \
| sort -k2gr > ${XLS}.rnk

############################################################
### The usage of the script is as follows:
#./rnkgen.sh /path/to/spreadsheet.xls GeneID_column FoldChange_column P-value_column
