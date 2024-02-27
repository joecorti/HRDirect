#!/bin/bash
set -e
set -u

if [ $# -ne 2 ]; then
	echo "Usage: $0 distribution bin_length" >&2
	echo "Output the normalized GSL-histogram output" >&2
	exit 78
fi

DISTR=$1
L=$2

eval $(echo $L | awk '{
			L=$1
			N=int(100/$1)+1;
			start=-1*L/2;
			end=100+L/2;
			print "N="N";", "start="start";", "end="end";"
		}')

cat $DISTR | \
	awk 'BEGIN{OFS="\t"}
		{
			x = ($1+$2)/2
			y = $3
			for(i=1;i<=y;i++) print x
		}' | \
	gsl-histogram -u -- $start $end $N
