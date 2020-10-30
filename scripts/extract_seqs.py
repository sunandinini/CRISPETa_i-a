import re, sys

#############################   Functions   ####################################
def myfindall(seq, candidates):
    resultlist = list()
    pos=0
    while pos + len(seq) <= len(candidates):
        result = re.search(seq, candidates[pos:])
        if result is None:
            pos = len(candidates)
        else:
            start = pos + result.start()
            end = pos + result.end()
            pos = pos + result.start() + 1
            resultlist.append(candidates[start:end])
    return resultlist

def get_genome_coordinates(header):
	h = header.split(":")[-1]
	chrom = header.split(":")[2]
	strand = h.split("(")[-1][0]
	start = h.split("(")[0].split("-")[0]
	end = h.split("(")[0].split("-")[1]

	return chrom, start, end, strand

def get_sgrna_info(m, header, strand):
	genome_coordinates = get_genome_coordinates(header)
	sgrna_chrom = genome_coordinates[0]
	sgrna_start = str(int(genome_coordinates[1])+m.start()+4)
	sgrna_end = str(int(genome_coordinates[1])+m.start()+30-3)
	sgrna_strand = strand
	if sgrna_strand == "-":
		sgrna_start = str(int(sgrna_start)-1)
		sgrna_end = str(int(sgrna_end)-1)
	sgrna_seq = m.group(1)

	return genome_coordinates, [sgrna_chrom, sgrna_start, sgrna_end, sgrna_seq, sgrna_strand]
################################################################################

seq1 = '.{24}.GG...'
seq2 = '...CC.{25}'

# Choose sequences from epifactor file
input_file = open(sys.argv[1], "r")
output_file = open(sys.argv[2], "w")

i = 0
for candidates in input_file:
	# Write the transcript ID
	if i%2 == 0:

		header = candidates
		output_file.write(candidates)
	# Write all sequences with a PAM
	else:
		for m in re.finditer(r'(?=(.{24}.GG...))',str(candidates.upper())):
			features = get_sgrna_info(m, header, "+")
			sgrna_features = features[1]
			output_file.write(",".join(sgrna_features)+ "\n")

		for m in re.finditer(r'(?=(...CC.{25}))',str(candidates.upper())):
			features = get_sgrna_info(m, header, "-")
			sgrna_features = features[1]
			output_file.write(",".join(sgrna_features)+ "\n")
	i = i+1
input_file.close()
output_file.close()


