# Original code written by carlospq. Use the bug-fix branch as it contains slight modifications to ensure correct output sequences for guide RNAs on the negative / anti-sense strand. 


# How to run CRISPETa-i/a

Update the file paths in "myfile.config"

bash ./CRISPRia.sh myfile.config			
							
# Requeriments:						
  - Python-2.7						
  - pickle					
  - pandas					
  - numpy						
  - scipy						
  - Bio.Seq					
  - _mysql					
  - MySQL (with off-targets database)			
  - bedtools (tested with bedtools v2.26.0)						
  - scikit-learn 0.16.1

# Other files required:

1. Fasta file for the genome of interest. 
2. Fantom TSS CAGE peaks (https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/)
   
# IMPORTANT: CRISPETa needs off-target scores to run! Find instructions here for how to set up the MySQL database for off-targets:

https://github.com/Carlospq/CRISPETa

# To test CRISPETa-i_a:

Output file containing top 10 CRISPRi designs along with the input "test.bed" has been provided in test_files. 

A fantom range of 10 was used with masked fasta and the fantom TSS file for the hg19 human genome (see link above). 
    
