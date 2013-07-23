#!/bin/bash

E_PATH="/data/rdata/enrichment"
RSCRIPT="/home/ubuntu/analysis_code/shared_code/build_minimum_enrichment_coverage.r"

Rscript $RSCRIPT $E_PATH/preMBT_pol_?_enrichment.cov.RData preMBT_pol_minimum_enrichment.cov.RData
Rscript $RSCRIPT $E_PATH/preMBT_tbp_?_enrichment.cov.RData preMBT_tbp_minimum_enrichment.cov.RData
Rscript $RSCRIPT $E_PATH/preMBT_k4_?_enrichment.cov.RData preMBT_k4_minimum_enrichment.cov.RData

Rscript $RSCRIPT $E_PATH/MBT_pol_?_enrichment.cov.RData MBT_pol_minimum_enrichment.cov.RData
Rscript $RSCRIPT $E_PATH/MBT_tbp_?_enrichment.cov.RData MBT_tbp_minimum_enrichment.cov.RData
Rscript $RSCRIPT $E_PATH/MBT_k4_?_enrichment.cov.RData MBT_k4_minimum_enrichment.cov.RData
Rscript $RSCRIPT $E_PATH/MBT_h3ac_?_enrichment.cov.RData MBT_h3ac_minimum_enrichment.cov.RData

