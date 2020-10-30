import re, sys, _mysql
from Bio.Seq import Seq
from operator import itemgetter
import subprocess

# Read filenames from command line
input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[2], 'w')

#############################   Functions   ####################################
def use_database():
	conn = _mysql.connect(user="crispeta",passwd="crispeta",host="",db="crispeta")
	return conn

def off_targets(sgrna):
	seq = sgrna[-2][4:-6]
	if sgrna[-1].strip() == "-":
		seq = str(Seq(sgrna[-2]).reverse_complement())
		seq = seq[4:-6]
	if 'N' in seq:
		return ['100', '100', '100', '100', '100']

	conn = use_database()
	conn.query("SELECT * FROM hg19 where grna ='"+seq+"';")
	result = conn.store_result()
	off = result.fetch_row()

	#if int(off[0][1]) <= 1 and int(off[0][2]) <= 0 and int(off[0][3]) <= 0 and int(off[0][4]) <= t[3] and int(off[0][5]) <= t[4]:
	if len(off) == 0:
		return ['100', '100', '100', '100', '100']
	else:
		return off[0][1:]

################################################################################

# Print header in output
output_file.write("\t".join(["target_id", "chr", "sgRNA_start", "sgRNA_end", "sgRNA_context", "sgRNA_strand", "off0", "off1", "off2", "off3", "\n"]))

# Read file and generate pairs
sgRNAs = []
for line in input_file:

	if line.startswith(">"):

		if len(sgRNAs) > 0:

			for sgRNA in sgRNAs:
				off = off_targets(sgRNA)
				sgRNA[-1] = sgRNA[-1].strip()
			
				output_file.write("\t".join([target]+sgRNA+list(off)+["\n"]))
			

		target = line.strip()
		sgRNAs = []
		continue

	sgRNA = line.split(",")
	sgRNAs.append(sgRNA)

for sgRNA in sgRNAs:
	off = off_targets(sgRNA)
	sgRNA[-1] = sgRNA[-1].strip()
	output_file.write("\t".join([target]+sgRNA+list(off)+["\n"]))

input_file.close()
output_file.close()


