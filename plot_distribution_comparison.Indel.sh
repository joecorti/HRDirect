set -e
set -u
set -o pipefail

if [ $# -lt 4 ]; then
	echo "Usage: $0 output_base metaNormal_Indel_distribution subtracted_Indel_distribution mirror_Indel_distribution [title]" >&2
	echo "Create the plot <output_base>.freq-distr_comparison_after-subtract.pdf." >&2
	exit 90
fi

OUT=$1.Indel.freq-distr_comparison_after-subtract
META_Indel=$2
SUBTRACT_Indel=$3
MIRROR_Indel=$4
TITLE=${5:-}

# crea la figura di confronto
cat <<EOF | gnuplot
set term post color solid noenh land
set out "$OUT.ps"
set grid
set title "$TITLE"
plot [1:] "$META_Indel" u (\$1+\$2)/2:3 w lp lw 6 t "metaNormal Indel", \
"$SUBTRACT_Indel" u (\$1+\$2)/2:3 w lp lw 6 t "subtract Indel", \
"$MIRROR_Indel" u (\$1+\$2)/2:3 w lp lw 3 dt 2 t "mirror Indel"

EOF
# convert
ps2pdf $OUT.ps $OUT.pdf
rm $OUT.ps
