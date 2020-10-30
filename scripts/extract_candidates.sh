#!/bin/bash


# Generate BED file with 220nt up/down stream of TSS
if [ $1 == "CRISPRa" ]; then
	cat $2 | awk 'BEGIN {OFS = "\t"} {if ($6=="+") {print $1,$2-320, $3, $4, $5, $6} else {print $1, $2, $3+320, $4, $5, $6} }' > tmp/FANTOM_increased.bed
fi

if [ $1 == "CRISPRi" ]; then
	cat $2 | awk 'BEGIN {OFS = "\t"} {if ($6=="+") {print $1,$2, $3+320, $4, $5, $6} else {print $1, $2-320, $3, $4, $5, $6} }' > tmp/FANTOM_increased.bed
fi

# Extract nucleotide sequence
bedtools getfasta -fi $3 -bed tmp/FANTOM_increased.bed -fo tmp/seqs.out -name





