#!/bin/bash

# log of test file generation
export windows="../../hg19/hg19_lite.250kb_windows"
../simulations/sconce_sim infileA.txt paramfile.txt 0 > ref_sim.log 2>&1
mkdir -p diploid
mkdir -p cancer
mv *healthy* diploid
find . -maxdepth 1 -name "*_cancer_cell_[1-9]*" -exec mv {} cancer \;
find . -name "*_cell_*[0-9]" | while read simu ; do
  depth="$simu"".hg19_lite.bed"
  paste "$windows" <(awk '{print $4}' $simu) | awk 'BEGIN{OFS="\t"} {if($4 == "") {$4=0} print}' > $depth
done

find diploid -name "simu_healthy*hg19_lite.bed" | sort -V > diploidFileList

Rscript --vanilla ../scripts/avgDiploid.R diploidFileList ref_healthy_avg.bed
Rscript --vanilla ../scripts/fitMeanVarRlnshp.R diploidFileList ref.meanVar

time ../sconce -d ref_healthy_avg.bed -k 5 -s simu_cancer_cell_0.hg19_lite.bed --meanVarCoefFile ref.meanVar -o ref_k5_simu_cancer_cell_0 > ref_k5_simu_cancer_cell_0.log 2> ref_k5_simu_cancer_cell_0.err

Rscript --vanilla ../scripts/plotGenomeTrace.R ref_healthy_avg.bed simu_cancer_cell_0.hg19_lite.bed ref_k5_simu_cancer_cell_0.bed ref_k5_simu_cancer_cell_0.png "Genome Trace for SCONCE (k=5)" true_cancer_cell_0.hg19_lite.bed

