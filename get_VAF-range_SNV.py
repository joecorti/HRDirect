import sys

if len(sys.argv) != 3+1:
	sys.exit("""
	Usage: %s start end SNV.file.vcf|-
	""" % sys.argv[0])

DEBUG = False

start = float(sys.argv[1])
end = float(sys.argv[2])
SNV_fd = (sys.argv[3] != '-') and open(sys.argv[3]) or sys.stdin

def get_FORMAT_ID(VCF_row, ID):
	FORMAT_string = VCF_row[8]
	TUMOR = VCF_row[10]
	idx = FORMAT_string.split(':').index(ID)
	return TUMOR.split(':')[idx]

def parse_VCF(fd):
# filter VCF and return frequency and the whole row
	for row in fd:
		if row[0] == '#':
			continue
		row = row.strip().split('\t')
		#tumor_info = row[-1].split(':')
		#support = int(tumor_info[-1])
		#depth =int( tumor_info[-2])
		VAF = 100*float(get_FORMAT_ID(row, "PM"))

		yield (VAF, row)

def read_ref_distribution(fd):
# read the reference distribution file in the format as the GSL-histogram output
	distr = []
	for row in fd:
		row = row.strip().split()
		start = float(row[0])
		end = float(row[1])
		value = float(row[2])
		distr.append((start,end, value))

	distr.sort() # sort by coordinate
	return distr

#############################

# read VCF and assign indels to bins
SNV = {(start, end): []}

for f, row in parse_VCF(SNV_fd):
	if start <= f <= end:
		SNV[(start, end)].append(row)

for k, v in SNV.items():
	for i in v:
		print('\t'.join(i))
