set -e
set -u
set -o pipefail

BIN_DIR=$(readlink -e $(dirname $0))
LOCAL_BIN_DIR=$(readlink -e $(dirname $0))
export BIN_DIR

if [ $# -ne 3 ]; then
	echo "Usage: $0 variant_dir tumor_sample L_bin" >&2
	echo "Recalibrate SNV and Indel files and recalculate HRD." >&2
	exit 78
fi

FILES_DIR=$(readlink -e $1)
TUMOR=$2
L_BIN=$3
PARALLEL=10 # parallel
N=50 # sampling number
export TUMOR

################ functions #############
function sampling_indel {
# $1: vcf.gz
# $2: out ($OUT_I/$(basename $1))
	(
		zcat -f $1 | \
			grep '^#'
		zcat -f $1 | \
			python3 $BIN_DIR/get_indel_from_freq-distr.py $INDEL_NT_DISTRIBUTION - | \
			sort -k1,1 -k2,2n
	) | \
	bgzip > $2
	tabix $2
}

function sampling_SNV {
# $1: vcf.gz
# $2: out ($OUT_S/$(basename $1))
	(
		zcat -f $1 | \
			grep '^#'
		zcat -f $1 | \
			python3 $BIN_DIR/get_SNV_from_freq-distr.py $SNV_NT_DISTRIBUTION - | \
			sort -k1,1 -k2,2n
	) | \
	bgzip > $2
	tabix $2
}
export -f sampling_indel sampling_SNV
#######################################################

RECALIBRATE=$FILES_DIR/recalibrateD
INDEL_SMPL=Indel-sampling
SNV_SMPL=SNV-sampling
mkdir $RECALIBRATE
cd $RECALIBRATE
mkdir $INDEL_SMPL $SNV_SMPL

# 0.
# get SNV_NT_DISTRIBUTION & INDEL_NT_DISTRIBUTION
#
SNV=$(ls ../SNV.*.vcf.gz)
bash $LOCAL_BIN_DIR/subtract_germline_distribution.SNV.sh $SNV SNV_distribution $L_BIN 
export SNV_NT_DISTRIBUTION=$(readlink -e SNV_distribution.subtracted)
INDEL=$(ls ../Indel.*.vcf.gz)
bash $LOCAL_BIN_DIR/subtract_germline_distribution.Indel.sh $INDEL Indel_distribution $L_BIN 
export INDEL_NT_DISTRIBUTION=$(readlink -e Indel_distribution.subtracted)

#
# 1.
# sampling on N-T distribution
#
# Indel
echo Sampling SNV and Indel results based on N-T frequency distribution
INDEL=$(ls ../Indel.*.vcf.gz)
for n in $(seq $N)
do
	echo $INDEL $INDEL_SMPL/$(basename $INDEL .gz).RANDOM-$n.gz
done | \
	parallel -j $PARALLEL --delay 0.1 --tmpdir . --colsep ' ' "sampling_indel {1} {2}" :::

# SNV
SNV=$(ls ../SNV.*.vcf.gz)

for n in $(seq $N)
do
	echo $SNV $SNV_SMPL/$(basename $SNV .gz).RANDOM-$n.gz
done | \
	parallel -j $PARALLEL --delay 0.1 --tmpdir . --colsep ' ' "sampling_SNV {1} {2}" :::

#
# 2.
# recreate an updated '.files' dir
#

echo Recreate an unpdated '.files' directory
indel_file=$(basename $(ls ../Indel.*gz))
indel_file_noExt=$(basename $indel_file .gz)

for rnd in $INDEL_SMPL/$indel_file_noExt.*.gz # all random
do
	SMP=$(basename $rnd | cut -d. -f2)
	X=$(basename $rnd .gz)
	RND="${X##*.}" # random number
	RND_DIR=$RND.files
	mkdir $RND_DIR
	# copy original files
	for i in CNV Indel SNV SV
	do
		cp -i ../$i.* $RND_DIR
	done
	# replace indels
	cp $rnd $RND_DIR/$indel_file
	cp $rnd.tbi $RND_DIR/$indel_file.tbi
	# replace SNV
	SNV_file=$(basename $(ls ../SNV.*gz))
	SNV_file_noExt=$(basename $SNV_file .gz)
	cp $SNV_SMPL/SNV.$SMP.annot.muts.vcf.$RND.gz $RND_DIR/$SNV_file
	cp $SNV_SMPL/SNV.$SMP.annot.muts.vcf.$RND.gz.tbi $RND_DIR/$SNV_file.tbi
done

#
# 3.
# calculate HRD for all RANDOM directory
#
echo Calculate HRD for all RANDOM directory
parallel -j $PARALLEL --delay 0.1 --tmpdir . --colsep ' ' "cd {1} && $BIN_DIR/HRDetect_fullPipeline-hg38.AUTO.sh $TUMOR ." ::: RANDOM-*.files

#
# 4.
# recalculate HRD based on median values from datamatrix
#
echo Finally, recalculate HRD based on median values from datamatrix coming from RANDOM directories

FIRST_RANDOM=RANDOM-1.files
DATA_MATRIX=$(basename $FILES_DIR .files).data-matrix
for out in RANDOM-*.files/HRDetect_fullPipeline.out
do
	bash $BIN_DIR/get_ordered_data_matrix.sh $out | \
	awk 'BEGIN{OFS="\t"}
		NR==1{
			for(i=1;i<=NF;i++){
				if($i=="del.mh.prop") MH_i=i+1
				if($i=="SNV3") SNV3_i=i+1
				if($i=="SNV8") SNV8_i=i+1
			}
		}
		NR==2{
			MH = $MH_i
			SNV3 = $SNV3_i
			SNV8 = $SNV8_i
			print MH, SNV3, SNV8
		}'
done | \
	datamash median 1 median 2 median 3 | \
	while read MH SNV3 SNV8
	do
		bash $BIN_DIR/get_ordered_data_matrix.sh $FIRST_RANDOM/HRDetect_fullPipeline.out | \
		awk -v MH=$MH -v SNV3=$SNV3 -v SNV8=$SNV8 'BEGIN{OFS="\t"}
			NR==1 {
				$1=$1
				print "", $0
				# del.mh.prop SNV3 SV3      SV5 hrd SNV8
				for(i=1;i<=NF;i++) {
						if($i=="del.mh.prop") MH_i=i+1
						if($i=="SNV3") SNV3_i=i+1
						if($i=="SNV8") SNV8_i=i+1
						if($i=="SV3") SV3_i=i+1
						if($i=="SV5") SV5_i=i+1
						if($i=="hrd") hrd_i=i+1
					}
			}
			NR==2 {  # sostituisco
				#print $1, $MH_i, $SNV3_i, $SNV8_i, SV3_i, SV5_i, hrd_i
				#print $1, MH, SNV3, SV3_i, SV5_i, hrd_i, SNV8 # XXX bug!!!!!!!!! Cancella SV3, SV5, hrd
				print $1, MH, SNV3, $SV3_i, $SV5_i, $hrd_i, SNV8
			}' > $DATA_MATRIX
	done

# activate Conda environment
#set +eu
#source .../etc/profile.d/conda.sh
#conda activate signature_tools_lib

Rscript $BIN_DIR/_HRDetect-score_from_data-matrix.R $DATA_MATRIX > HRDetect_fullPipeline.out 2> HRDetect_fullPipeline.err

# deactivate Conda env
#conda deactivate # signature_tools_lib
#conda deactivate # base
#set -eu

echo Done.
