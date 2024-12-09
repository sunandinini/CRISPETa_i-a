# Original code written by carlospq. Slight modifications have been made to it to ensure correct output sequences for guide RNAs on the negative / anti-sense strand. 


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

# IMPORTANT: CRISPETa needs off-target scores to run! 
# Find instructions here for how to set up the MySQL database for off-targets:

https://github.com/Carlospq/CRISPETa

# To test output:

Output file containing top 10 CRISPRi designs for test.bed in has been provided in test_files. A fantom range of 10 was used with hg19_masked.fa
    
