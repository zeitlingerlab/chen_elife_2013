#!/bin/bash

BASE_PATH="/data/rdata"
ETRACK_PATH="/home/ubuntu/analysis_code/shared_code"

# pre-MBT samples

# Pol II 
for i in 1 2 3 4
do
  echo Rscript $ETRACK_PATH/enrichment_track.r -t $BASE_PATH/preMBT_pol_${i}.cov.RData -c $BASE_PATH/preMBT_wce_1.cov.RData -n preMBT_pol_${i}_enrichment -s -l log2
done

# TBP
for i in 1 2
do
  echo Rscript $ETRACK_PATH/enrichment_track.r -t $BASE_PATH/preMBT_tbp_${i}.cov.RData -c $BASE_PATH/preMBT_wce_1.cov.RData -n preMBT_tbp_${i}_enrichment -s -l log2
done

# K4
for i in 1 2
do
  echo Rscript $ETRACK_PATH/enrichment_track.r -t $BASE_PATH/preMBT_k4_${i}.cov.RData -c $BASE_PATH/preMBT_wce_1.cov.RData -n preMBT_k4_${i}_enrichment -s -l log2
done

# MBT samples

# Pol II
for i in 1 2 3 
do
  echo Rscript $ETRACK_PATH/enrichment_track.r -t $BASE_PATH/MBT_pol_${i}.cov.RData -c $BASE_PATH/MBT_wce_1.cov.RData -n MBT_pol_${i}_enrichment -s -l log2
done

# TBP
for i in 1 2
do
  echo Rscript $ETRACK_PATH/enrichment_track.r -t $BASE_PATH/MBT_tbp_${i}.cov.RData -c $BASE_PATH/MBT_wce_1.cov.RData -n MBT_tbp_${i}_enrichment -s -l log2
done

# K4
for i in 1 2
do
  echo Rscript $ETRACK_PATH/enrichment_track.r -t $BASE_PATH/MBT_k4_${i}.cov.RData -c $BASE_PATH/MBT_wce_1.cov.RData -n MBT_k4_${i}_enrichment -s -l log2
done

# H3Ac
for i in 1 2
do
  echo Rscript $ETRACK_PATH/enrichment_track.r -t $BASE_PATH/MBT_h3ac_${i}.cov.RData -c $BASE_PATH/MBT_wce_1.cov.RData -n MBT_h3ac_${i}_enrichment -s -l log2
done

