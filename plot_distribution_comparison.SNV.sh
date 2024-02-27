set -e
set -u
set -o pipefail

if [ $# -lt 4 ]; then
	echo "Usage: $0 output_base metaNormal_SNV_distribution subtracted_SNV_distribution mirror_SNV_distribution [title]" >&2
	echo "Create the plot <output_base>.freq-distr_comparison_after-subtract.pdf." >&2
	exit 90
fi

OUT=$1.SNV.freq-distr_comparison_after-subtract
META_SNV=$2
SUBTRACT_SNV=$3
MIRROR_SNV=$4
TITLE=${5:-}

# crea la figura di confronto
cat <<EOF | gnuplot
set term post color solid noenh land
set out "$OUT.ps"
set grid
set title "$TITLE"
plot [1:] "$META_SNV" u (\$1+\$2)/2:3 w lp lw 6 t "metaNormal SNV", \
"$SUBTRACT_SNV" u (\$1+\$2)/2:3 w lp lw 6 t "subtract SNV", \
"$MIRROR_SNV" u (\$1+\$2)/2:3 w lp lw 3 dt 2 t "mirror SNV"

EOF
# convert
ps2pdf $OUT.ps $OUT.pdf
rm $OUT.ps
