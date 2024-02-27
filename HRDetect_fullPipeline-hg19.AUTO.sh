#!/bin/bash
#set -e
#set -u
if [ $# -ne 2 ]; then
	echo "Usage: $0 TUMOR_sample_name variant_dir" >&2
	exit 8
fi

SMP=$1
DIR=$2
BIN_DIR=$(readlink -e $(dirname $0))


SNV=$(readlink -e $DIR/SNV.*_vs_*.annot.muts.vcf.gz)
SV=$(readlink -e $DIR/SV.*_vs_*.annot.bedpe.OK.gz)
INDEL=$(readlink -e $DIR/Indel.*_vs_*.annot.vcf.gz)
CNV=$(readlink -e $DIR/CNV.*.copynumber.caveman.OK)

if [ -e "$DIR/HRDetect_fullPipeline.out" ]; then
	echo "ERROR: output file exists: $DIR/HRDetect_fullPipeline.out. Exit." >&2
	exit 8
fi

cd $DIR

# activate Conda environment
#set +eu
#source ..../etc/profile.d/conda.sh
#conda activate signature_tools_lib

Rscript $BIN_DIR/HRDetect_fullPipeline-hg19.R $SMP $SNV $SV $INDEL $CNV > HRDetect_fullPipeline.out 2> HRDetect_fullPipeline.err

# deactivate Conda env
#conda deactivate # signature_tools_lib
#conda deactivate # base
#set -eu
