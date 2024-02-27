#!/bin/bash
set -e
set -u

if [ $# -ne 2 ]; then
	echo "Usage: $0 variant_dir tumor_sample" >&2
	exit 7
fi

DIR=$1
SMP=$2

BINDIR=$(readlink -e $(dirname $0))

# HRDetect
if [ ! -e "$DIR/HRDetect_fullPipeline.out" ]; then
	echo "...................... Running classic HRDetect_fullPipeline ................" >&2
	bash $BINDIR/HRDetect_fullPipeline-hg38.AUTO.sh $SMP $DIR 
else
	echo "WARNING: file exists: $DIR/HRDetect_fullPipeline.out. Skipping." >&2
fi

# strategy A
if [ ! -e "$DIR/recalibrateA" ]; then
	echo "" >&2
	echo "...................... Running strategy A ................" >&2
	bash $BINDIR/strategy_A.sh $DIR $SMP
else
	echo "WARNING: directory exists: $DIR/recalibrateA. Skipping." >&2
fi

# strategy B
if [ ! -e "$DIR/recalibrateB" ]; then
	echo "" >&2
	echo "...................... Running strategy B ................" >&2
	bash $BINDIR/strategy_B.sh $DIR $SMP $DIR/recalibrateB
else
	echo "WARNING: directory exists: $DIR/recalibrateB. Skipping." >&2
fi

# strategy C
if [ ! -e "$DIR/recalibrateC" ]; then
	echo "" >&2
	echo "...................... Running strategy C ................" >&2
	bash $BINDIR/strategy_C.sh $DIR $SMP
else
	echo "WARNING: directory exists: $DIR/recalibrateC. Skipping." >&2
fi

# strategy D
if [ ! -e "$DIR/recalibrateD" ]; then
	echo "" >&2
	echo "...................... Running strategy D ................" >&2
	bash $BINDIR/strategy_D.sh $DIR $SMP 5
else
	echo "WARNING: directory exists: $DIR/recalibrateD. Skipping." >&2
fi

# Ensemble
echo "" >&2
echo "...................... Running Ensemble byC................" >&2
bash $BINDIR/Ensemble_strategy-byC.sh $DIR
echo "" >&2
echo "...................... Running Ensemble byB................" >&2
bash $BINDIR/Ensemble_strategy-byB.sh $DIR
