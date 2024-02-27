set -e
set -u

if [ $# -ne 3 ]; then
	echo "Usage: $0 metaNormal.SNV.vcf[.gz] outputBase L_bin" >&2
	exit 6
fi

META=$1
OUT=$2
L_BIN=$3
BIN_DIR=$(readlink -e $(dirname $0))

# get normalized metaNormal distribution
bash $BIN_DIR/get_SNV_freq-distr_anyPerc.NORM.sh $META $L_BIN > $OUT.meta.norm
# get metaNormal distribution
bash $BIN_DIR/get_SNV_freq-distr_anyPerc.sh      $META $L_BIN > $OUT.meta

# get mirror distribution
python3 $BIN_DIR/_get_mirror_distribution.py $OUT.meta > $OUT.mirror
bash $BIN_DIR/get_freq-distr_anyPerc_from_distribution.NORM.sh $OUT.mirror $L_BIN > $OUT.mirror.norm

# subtract mirror distribution from the metaNormal one
python3 $BIN_DIR/_subtract_distribution.py $OUT.meta $OUT.mirror > $OUT.subtracted

# get normalized subtracted distribution
bash $BIN_DIR/get_freq-distr_anyPerc_from_distribution.NORM.sh $OUT.subtracted $L_BIN > $OUT.subtracted.norm

#	# plot
bash $BIN_DIR/plot_distribution_comparison.SNV.sh $OUT.norm    $OUT.meta.norm $OUT.subtracted.norm $OUT.mirror.norm "$OUT - Normalized"
bash $BIN_DIR/plot_distribution_comparison.SNV.sh $OUT.raw     $OUT.meta $OUT.subtracted $OUT.mirror "$OUT - NON normalized"
