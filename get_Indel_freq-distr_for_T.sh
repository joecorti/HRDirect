set -e
set -u

if [ $# -ne 1 ]; then
	echo "Usage: $0 Indel.vcf[.gz]" >&2
	echo "Output the normalized GSL-histogram output" >&2
	exit 78
fi

INDEL=$1

zcat -f $INDEL | \
	grep -v '^#' | \
	awk 'BEGIN{OFS="\t"}
		{
			n1=split($(NF-1),n,":")
			n2=split($(NF),t,":")
			print 100*t[n2]/t[n2-1]
		}' | \
	gsl-histogram -u -- -0.5 100.5
