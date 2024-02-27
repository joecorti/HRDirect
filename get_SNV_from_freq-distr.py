import sys
import random
import numpy as np

if len(sys.argv) != 2+1:
	sys.exit("""
	Usage: %s reference_freq-distr SNV.file.vcf|-
	""" % sys.argv[0])

DEBUG = False

REF_fd = open(sys.argv[1])
SNV_fd = (sys.argv[2] != '-') and open(sys.argv[2]) or sys.stdin

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

# read reference distribution
distr = read_ref_distribution(REF_fd)

distr_dict = {}
for i in distr:
	distr_dict[(i[0], i[1])] = i[2]
bins = sorted(distr_dict.keys())

if DEBUG: print("Bins:", bins, file=sys.stderr)

# read VCF and assign indels to bins
SNV = {}
for i in bins:
	SNV[i] = []

for f, row in parse_VCF(SNV_fd):
	for b in bins:
		if f <= b[1]:
			SNV[b].append(row)
			break

# clever "support_unit"
ratios = [len(SNV[(start, end)])/value for (start, end, value) in distr if value != 0]
non_zero_ratios = [i for i in ratios if i!=0]
if non_zero_ratios == []:
	support_unit = 0
	print("WARNING: empty ratios list.... setting support_unit to 0.", file=sys.stderr)
else:
	support_unit = min(non_zero_ratios)

if DEBUG: print("support_unit:", support_unit, file=sys.stderr)

for b in bins:
	N = int(support_unit*distr_dict[b])
	if DEBUG:
		if N<=len(SNV[b]):
			print("From bin {} {} extract {} calls, instead of total {}.".format(
				b,
				distr_dict[b],
				N,
				len(SNV[b])
				), file=sys.stderr)
		else:
			print("From bin {} {} extract {} calls, instead of total {} <<<< we'd like {}.".format(
				b,
				distr_dict[b],
				len(SNV[b]),
				len(SNV[b]),
				N
				), file=sys.stderr)

	if N<=len(SNV[b]):
		sampled_SNV = random.sample(SNV[b], N)
	else:
		sampled_SNV = SNV[b]

	for i in sampled_SNV:
		print('\t'.join(i))
