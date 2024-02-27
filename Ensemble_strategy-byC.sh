#!/bin/bash
set -e
set -u

if [ $# -lt 1 ]; then
	echo "Usage: $0 variant_dir [variant_dir_2 ....]" >&2
	echo "Collect score for Ensemble Strategy" >&2
	echo ">> Based on strategy C <<" >&2
	exit 6
fi

OUT=recalibrate_Ensemble-byC

for d in $*
do
	mkdir -p $d/$OUT
	StrategyC=$d/recalibrateC
        if [ -e "$StrategyC" ]; then
               ensemble=$(awk '{if(NF>0) score = $NF}END{ensemble=0; if(score<0.25) ensemble=1; if(score>0.75) ensemble=2; print ensemble}' $StrategyC/HRDetect_fullPipeline.out)
	fi
	if [ $ensemble -ne 2 ]; then #NON HRD (HRP oppure NA)
		# strategy B
		StrategyB=$d/recalibrateB
		if [ -e "$StrategyB" ]; then
			R1=$(awk '{if(NF>0) score = $NF}END{result=0; if(score<0.25) result=1; if(score>0.75) result=-1; print result}' $StrategyB/out_HRD/INDEL-SLIDE-5-35.SNV-SLIDE-5-35.HRDetect_fullPipeline.out)
		fi
		# strategy A
		StrategyA=$d/recalibrateA
		if [ -e "$StrategyA" ]; then
			R3=$(awk '{if(NF>0) score = $NF}END{result=0; if(score<0.25) result=1; if(score>0.75) result=-1; print result}' $StrategyA/HRDetect_fullPipeline.out)
		fi
		# strategy D
        	StrategyD=$d/recalibrateD
        	if [ -e "$StrategyD" ]; then
                	R4=$(awk '{if(NF>0) score = $NF}END{result=0; if(score<0.25) result=1; if(score>0.75) result=-1; print result}' $StrategyD/HRDetect_fullPipeline.out)
        	fi

		echo $R1 $R3 $R4 | awk 'BEGIN {FS = " "};{sum=$1+$2+$3+$4}END{if(sum>0) print "0.1";if(sum<0) print "0.9";if(sum==0) print "0.5"}'> $d/$OUT/HRDetect_fullPipeline.out 
	else 
		echo "0.9"> $d/$OUT/HRDetect_fullPipeline.out  #HRD
	fi
done


