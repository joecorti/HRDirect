set -e
set -u
set -o pipefail

export BIN_DIR=$(readlink -e $(dirname $0))

if [ $# -ne 3 ]; then
	echo "Usage: $0 variant_dir tumor_sample out_dir" >&2
	echo "Recalibrate SNV and Indel files and recalculate HRD." >&2
	exit 78
fi

FILES_DIR=$(readlink -e $1)
TUMOR=$2
OUT_DIR=$3
L=30 # VAF range
SHIFT=5 # sliding window shift
PARALLEL=10 # parallel
N=11 # sampling number
export TUMOR
mkdir -p $OUT_DIR

################ functions #############
function sliding_indel {
# $1: start
# $2: end
# $3: vcf.gz
# $4: out ($OUT_I/$(basename $1))
	(
		zcat -f $3 | \
			grep '^#'
		zcat -f $3 | \
			python3 $BIN_DIR/get_VAF-range_indel.py $1 $2 - | \
			sort -k1,1 -k2,2n
	) | \
	bgzip > $4
	tabix $4
}

function sliding_SNV {
# $1: start
# $2: end
# $3: vcf.gz
# $4: out ($OUT_S/$(basename $1))
	(
		zcat -f $3 | \
			grep '^#'
		zcat -f $3 | \
			python3 $BIN_DIR/get_VAF-range_SNV.py $1 $2 - | \
			sort -k1,1 -k2,2n
	) | \
	bgzip > $4
	tabix $4
}

function run_HRDetect {
# $1: tumor sample
# $2: SNV
# $3: SV
# $4: Indel
# $5: CNV
	local SMP=$1
	local SNV=$(readlink -e $2)
	local SV=$(readlink -e $3)
	local INDEL=$(readlink -e $4)
	local CNV=$(readlink -e $5)
	# activate Conda environment
	#set +eu
	#source .../etc/profile.d/conda.sh
	#conda activate signature_tools_lib

	Rscript $BIN_DIR/HRDetect_fullPipeline-hg38.R $SMP $SNV $SV $INDEL $CNV

	# deactivate Conda env
	#conda deactivate # signature_tools_lib
	#conda deactivate # base
	#set -eu
}

export -f sliding_indel sliding_SNV run_HRDetect
#######################################################

INDEL_SLIDING=Indel-sliding
SNV_SLIDING=SNV-sliding
OUT_HRD=out_HRD
cd $OUT_DIR
mkdir $INDEL_SLIDING $SNV_SLIDING $OUT_HRD

#
# 1. sliding
#
# Indel
echo Strategy B: SNV and Indel results.............
INDEL=$(ls $FILES_DIR/Indel.*.vcf.gz)
for start in $(seq 0 $SHIFT $((100-L)))
do
	end=$(($start+$L))
	echo $start $end $INDEL $INDEL_SLIDING/$(basename $INDEL .gz).SLIDE-$start-$end.gz
done | \
	parallel -j $PARALLEL --delay 0.1 --tmpdir . --colsep ' ' "sliding_indel {1} {2} {3} {4}" :::

# SNV
SNV=$(ls $FILES_DIR/SNV.*.vcf.gz)
for start in $(seq 0 $SHIFT $((100-L)))
do
	end=$(($start+$L))
	echo $start $end $SNV $SNV_SLIDING/$(basename $SNV .gz).SLIDE-$start-$end.gz
done | \
	parallel -j $PARALLEL --delay 0.1 --tmpdir . --colsep ' ' "sliding_SNV {1} {2} {3} {4}" :::


#
# run HRDetect
#
echo Running HRDetect................
SV=$(readlink -e $FILES_DIR/SV.*.gz)
CNV=$(readlink -e $FILES_DIR/CNV.*.OK)

for indel in $(readlink -e $INDEL_SLIDING/Indel.*.gz)
do
	range=$(basename $indel | awk '{n=split($1,a,"."); print a[n-1]}')
	snv=$(readlink -e $SNV_SLIDING/SNV.*$range.gz)
	OUTBASE=INDEL-$range.SNV-$range
	OUT=$(readlink -f $OUT_HRD/$OUTBASE)
	for try_n in $(seq $N)
	do
		echo tmp-sliding.$OUTBASE.$try_n $TUMOR $snv $SV $indel $CNV $OUT $try_n
	done
done | \
	parallel -j $PARALLEL --delay 0.1 --tmpdir . --colsep ' ' "mkdir {1} && cd {1} && run_HRDetect {2} {3} {4} {5} {6} > {7}.HRDetect_fullPipeline.out.{8} 2> {7}.HRDetect_fullPipeline.err{8}" :::

# median
for indel in $(readlink -e $INDEL_SLIDING/Indel.*.gz)
do
	range=$(basename $indel | awk '{n=split($1,a,"."); print a[n-1]}')
	OUTBASE=INDEL-$range.SNV-$range
	OUT=$(readlink -f $OUT_HRD/$OUTBASE)
	for i in $OUT.HRDetect_fullPipeline.out.*
	do
		awk '{score=$NF}END{print score+0}' $i
	done | \
		datamash median 1 mean 1 min 1 max 1 > $OUT.HRDetect_fullPipeline.out
done

rm -r tmp-sliding.*

exit
