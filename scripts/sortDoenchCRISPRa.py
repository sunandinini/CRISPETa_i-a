import re, sys
import subprocess
from operator import itemgetter
import operator

# read input and output files from command line
input_file = open(sys.argv[1], 'r')
strand_file = open(sys.argv[2], 'r')
output_file = open(sys.argv[3], 'w')

# options from command line
mode = sys.argv[4]
submode = sys.argv[5]

#number of top sgRNAs returned
n_gRNA = int(sys.argv[6])


## CRISPRa ##
#Picking_Round	Position_relative_to_TSS	Off_target_criteria	On_target_criteria-Rule_Set2_score
#	1				-150 to -75				<=1 perfect match			>=0.2
#	2				-200 to -25				<=1 perfect match			>=0.2
#	3				-250 to 0				<=1 perfect match			>=0.2
#	4				-150 to -75				<=5 perfect match			>=0.2
#	5				-200 to -25				<=5 perfect match			>=0.2
#	6				-250 to 0				<=5 perfect match			>=0.2
#	7				-150 to -75				<=5 perfect match			>0
#	8				-200 to -25				<=5 perfect match			>0
#	9				-250 to 0				<=5 perfect match			>0
#	10				-250 to 0				<=1 perfect match			>0
#	11				-300 to 0				<=1 perfect match			>0
#	12				-300 to +300			Any 						Any

## CRISPRi ##
#Picking_Round	Position_relative_to_TSS	Off_target_criteria	On_target_criteria-Rule_Set2_score
#	1				+25 to +75				<=1 perfect match			>=0.2
#	2				0 to +100				<=1 perfect match			>=0.2
#	3				+175 to +250			<=1 perfect match			>=0.2
#	4				+25 to +75				<=5 perfect match			>=0.2
#	5				0 to +100				<=5 perfect match			>=0.2
#	6				+175 to +250			<=5 perfect match			>=0.2
#	7				+25 to +75				<=5 perfect match			>0
#	8				0 to +100				<=5 perfect match			>0
#	9				+175 to +250			<=5 perfect match			>0
#	10				-50 to 300				<=1 perfect match			>0
#	11				-50 to +300				<=1 perfect match			>0
#	12				-300 to +300			Any 						Any

def dist_to_TSS(sgRNA, targetStrand):

	StartEndTSS = int(sgRNA[0].split(":")[-1].split("-")[1])

	#Cutting possition of the sgRNA
	if sgRNA[5] == '+':
		sgRNAPosition = int(sgRNA[2])+17
	else:
		sgRNAPosition = int(sgRNA[2])+5

	#distance in absolute values
	dist = max(StartEndTSS, sgRNAPosition) - min(StartEndTSS, sgRNAPosition)

	return dist

def Picking_Round(method, dist, offTargets, score):

	if method == 'CRISPRa':
		if dist >= 75 and dist <= 150 and offTargets <= 1 and score >= 0.2:
			ranking = 1
		elif dist >= 25 and dist <= 200 and offTargets <= 1 and score >= 0.2:
			ranking = 2
		elif dist >= 0 and dist <= 250 and offTargets <= 1 and score >= 0.2:
			ranking = 3
		elif dist >= 75 and dist <= 150 and offTargets <= 5 and score >= 0.2:
			ranking = 4
		elif dist >= 25 and dist <= 200 and offTargets <= 5 and score >= 0.2:
			ranking = 5
		elif dist >= 0 and dist <= 250 and offTargets <= 5 and score >= 0.2:
			ranking = 6
		elif dist >= 75 and dist <= 150 and offTargets <= 5 and score >= 0:
			ranking = 7
		elif dist >= 25 and dist <= 200 and offTargets <= 5 and score >= 0:
			ranking = 8
		elif dist >= 0 and dist <= 250 and offTargets <= 5 and score >= 0:
			ranking = 9
		elif dist >= 0 and dist <= 250 and offTargets <= 1 and score >= 0:
			ranking = 10
		elif dist >= 0 and dist <= 300 and offTargets <= 1 and score >= 0:
			ranking = 11
		elif dist >= 0 and dist <= 320:
			ranking = 12


		return ranking

	if method == 'CRISPRi':
		if dist >= 25 and dist <= 75 and offTargets <= 1 and score >= 0.2:
			ranking = 1
		elif dist >= 0 and dist <= 100 and offTargets <= 1 and score >= 0.2:
			ranking = 2
		elif dist >= 175 and dist <= 250 and offTargets <= 1 and score >= 0.2:
			ranking = 3
		elif dist >= 25 and dist <= 75 and offTargets <= 5 and score >= 0.2:
			ranking = 4
		elif dist >= 0 and dist <= 100 and offTargets <= 5 and score >= 0.2:
			ranking = 5
		elif dist >= 175 and dist <= 250 and offTargets <= 5 and score >= 0.2:
			ranking = 6
		elif dist >= 25 and dist <= 75 and offTargets <= 5 and score >= 0:
			ranking = 7
		elif dist >= 0 and dist <= 100 and offTargets <= 5 and score >= 0:
			ranking = 8
		elif dist >= 175 and dist <= 250 and offTargets <= 5 and score >= 0:
			ranking = 9
		elif dist >= -50 and dist <= 300 and offTargets <= 1 and score >= 0:
			ranking = 10
		elif dist >= -50 and dist <= 300 and offTargets <= 1 and score >= 0:
			ranking = 11
		elif dist >= -300 and dist <= 320:
			ranking = 12

		return ranking

def sgRNAs_single(dic, output_file=output_file, n_gRNA=n_gRNA):
	
	for target in dic:
		n = 0

		for sgRNA in dic[target]:
		 	if n == n_gRNA:
		 		continue

		 	for i in range(0,len(sgRNA)):
		 	 	sgRNA[i] = str(sgRNA[i])
			
		 	output_file.write("\t".join(sgRNA+[str(n)]+["\n"]))
		 	n += 1

def sgRNAs_paired(dic, output_file=output_file, n_gRNA=n_gRNA, mode=mode):

	if mode == "CRISPRa":
		limit = 150
	elif mode == "CRISPRi":
		limit = 150

	for target in dic:
		n = 0
		RNA1 = []
		RNA2 = []
		pairs = []
		counter = {}

		for sgRNA in dic[target]:
			if int(sgRNA[12]) > limit:
				RNA2.append(sgRNA)
			else:
				RNA1.append(sgRNA)

		for g1 in RNA1:
			for g2 in RNA2:
				pr1=int(g1[-1])
				pr2=int(g2[-1])
				counter[g1[4]] = 0
				counter[g2[4]] = 0
				pairs.append(g1+g2+[pr1+pr2])
		
		## MODIFY this sorting; do it by its individual picking round and not by the combination of both
		#pairs.sort(key = lambda x: x[28])
		pairs.sort(key = lambda x: (x[13],x[27]))
		#pairs.sort(key = operator.itemgetter(1, 2))

		for pair in pairs:
			if n >= n_gRNA:
				continue
			
			#Number of times the same sgRNA can appear in results
			if counter[pair[4]] >= 1 or counter[pair[18]] >= 1:
				continue
			counter[pair[4]] += 1
			counter[pair[18]] += 1
		
			for i in range(0,len(pair)):
			 	pair[i] = str(pair[i])
			output_file.write("\t".join(pair+[str(n)]+["\n"]))
			
			n += 1
			
# Print header in output
output_file.write("\t".join(["target_id", "chr", "sgRNA_start", "sgRNA_end", "sgRNA_context", "sgRNA_strand", "off0", "off1", "off2", "off3", "off4", "score_rule_set2", "dist_to_TSS", "Picking_Round", "\n"]))

# Get strand of targets
GeneStrands = {}
for line in strand_file:
	line = line.strip().split("\t")
	GeneStrands[line[3]] = line[5]
strand_file.close()

# Get sgRNAs for each target
dic = {}
for line in input_file:

	if line.startswith("target_id"): continue
	line = line.strip().split("\t")

	target = line[0]
	if target in dic:
		dic[target].append(line)
	else:
		dic[target] = [line]
input_file.close()

# Rank sgRNAs [ and pair results ]
for target in dic:

	x = [2,3,6,7,8,9,10,11]

	targetID = target.split(":")[0].replace(">","")
	targetStrand = GeneStrands[targetID]

	for sgRNA in dic[target]:
	 	for i in x:
	 		if i == 11:
	 			sgRNA[i] = float(sgRNA[i])
	 		else:
	 			sgRNA[i] = int(sgRNA[i])
	
		dist = dist_to_TSS(sgRNA, targetStrand)
		offTargets = int(sgRNA[6])
		score = float(sgRNA[-1])

		ranking = Picking_Round(mode, dist, offTargets, score)
		
		sgRNA.append(dist)
		sgRNA.append(ranking)

	dic[target] = [item for item in dic[target] if item[-1] < 12]
	dic[target].sort(key=lambda x: x[13])

if submode == "single":
	sgRNAs_single(dic)

if submode == "paired":
	sgRNAs_paired(dic)


output_file.close()
