#!/bin/bash

###################### CRISPRia v1 ######################
#							#
# New sorting method:					#
# - Based on Doench et al. 2018 (Up, Down and Out...)	#
#							#
# Command example:					#
# - bash ./CRISPRia.sh myfile.config			#
#							#
# Requeriments:						#
# - Python-2.7						#
#	- pickle					#
#	- pandas					#
#	- numpy						#
#	- scipy						#
#	- Bio.Seq					#
#	- _mysql					#
# - MySQL (with off-targets database)			#
# - bedtools						#
# - scikit-learn 0.16.1					#
#							#
#########################################################

# Load config file
. $1

# Create tmp folder
cd "${0%/*}"
mkdir -p tmp
mkdir -p results
mkdir -p plots

# Gather scripts dir
scripts_dir=$( echo $0 | sed 's/CRISPRia\.sh/scripts/' )

# Change to FANTOM TSS when distance < dist(nt)
if [ -e "results/input_to_FANTOM.bed" ]; then
	echo "results/input_to_FANTOM.bed already exist. Skipping: Change to FANTOM TSS"
else
	bash $scripts_dir/bed_to_FANTOM.sh -i $INFILE -f $FANTOM -r $RANGE > results/input_to_FANTOM.bed 2>/dev/null
fi

# Extract genomic sequences up/downstream of TSS (for CRISPRa and CRISPRi respectively)
if [ -e "tmp/seqs.out" ]; then
	echo "tmp/seqs.out already exist. Skipping: Extract genomic sequences"
else
	bash $scripts_dir/extract_candidates.sh $MODE results/input_to_FANTOM.bed $FASTA
fi

# Extract sgRNAs for each genomic region
if [ -e "results/sgRNAs.txt" ]; then
	echo "results/sgRNAs.txt already exist. Skipping: Extract sgRNAs"
else
	python $scripts_dir/extract_seqs.py tmp/seqs.out results/sgRNAs.txt 
fi

# Obtain off-target information for each sgRNA
if [ -e "results/sgRNAs_w_offtargets.txt" ]; then
	echo "results/sgRNAs_w_offtargets.txt already exist. Skipping: Obtaining off-target information"
else
	python $scripts_dir/off_targets.py results/sgRNAs.txt results/sgRNAs_w_offtargets.txt
fi

# Score each sgRNA with Rule_Set_2
if [ -e "results/sgRNAs_avaluated.txt" ]; then
	echo "results/sgRNAs_avaluated.txt already exist. Skipping: Scoring sgRNAs"
else
	python $scripts_dir/score.py results/sgRNAs_w_offtargets.txt results/sgRNAs_avaluated.txt $scripts_dir
fi

# Sort results
python $scripts_dir/sortDoenchCRISPRa.py results/sgRNAs_avaluated.txt tmp/FANTOM_increased.bed results/$MODE"_"$SUBMODE"_top_"$n_gRNAs"_selected.txt" $MODE $SUBMODE $n_gRNAs


# Result to BED format
# For $SUBMODE == paired
if [ $SUBMODE == "paired" ]; then
	echo "track name='paired_CRISPRi_designs' description=' ' visibility=3 itemRgb='On'" > results/$MODE"_"$SUBMODE"_designs.bed";
	awk 'NR>1 {{
			if ( $3 > $17 ){
				MIN=$17
			}
			else {
				MIN=$3
			}
			
			if ( $4 > $18 ){
				MAX=$4
			}
			else {
				MAX=$18
			}
		  }
		  {OFS="\t"; print $2, MIN, MAX, $1, $NF, ".", MIN, MAX, "150.0,0,0", 2, "20,20", "0,"MAX-MIN-20}}' results/$MODE"_"$SUBMODE"_top_"$n_gRNAs"_selected.txt" >> results/$MODE"_"$SUBMODE"_designs.bed";
fi  ##Changed "23,23" to "20,20" and "MAX-MIN-23 to MAX-MIN-20"

# For $SUBMODE == single
if [ $SUBMODE == "single" ]; then
	echo "track name='single_CRISPR_designs' description=' ' visibility=3 itemRgb='On'" > results/$MODE"_"$SUBMODE"_designs.bed";
	awk 'NR>1 {OFS="\t"; print $2, $3, $4, $1, $NF, $6}' results/$MODE"_"$SUBMODE"_top_"$n_gRNAs"_selected.txt" >> results/$MODE"_"$SUBMODE"_designs.bed";
fi



