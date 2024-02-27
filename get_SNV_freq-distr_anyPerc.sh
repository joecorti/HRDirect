#!/bin/bash
set -e
set -u

if [ $# -ne 2 ]; then
	echo "Usage: $0 SNV.vcf[.gz] bin_length" >&2
	echo "Output the GSL-histogram output" >&2
	exit 78
fi

SNV=$1
L=$2

eval $(echo $L | awk '{
			L=$1
			N=int(100/$1)+1;
			start=-1*L/2;
			end=100+L/2;
			print "N="N";", "start="start";", "end="end";"
		}')

zcat -f $SNV | \
	grep -v '^#' | \
	awk 'BEGIN{OFS="\t"}
		{
			n1=split($(NF-1),n,":")
			n2=split($(NF),t,":")
			print 100*t[n2]
		}' | \
	gsl-histogram -- $start $end $N
