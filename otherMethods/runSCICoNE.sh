#!/bin/bash

# Wed 03 Nov 2021 01:59:05 PM PDT
# script to run SCICoNE

export windows="/space/s2/sandra/hg19/hg19_lite.250kb_windows"
export chrBoundaries="$windows"".chrBoundaries"
export sciconeBin="/space/s2/sandra/methodComparisons/SCICoNE/build"

sciconeDir="$1"
depthFiles="$2" # needed to combine into an input csv file
tumorIDs="$3" # needed to split back into indv files
outputBase="$4"

mkdir -p $sciconeDir
# takes a cell x bins matrix, comma separated
if [ ! -f "$sciconeDir"/"$outputBase".csv ] ; then
  Rscript --vanilla "$scripts"/combineDepthsIntoCsv.R "$depthFiles" "$sciconeDir"/"$outputBase".csv
fi

# first run breakpoint detection
# input_breakpoints_file should include the the index of the last chr, but not the beginning index of the first
n_bins=$(wc -l < "$windows")
n_cells=$(wc -l < "$depthFiles")
time $sciconeBin/breakpoint_detection --d_matrix_file "$sciconeDir"/"$outputBase".csv --n_bins $n_bins --n_cells $n_cells --postfix "$sciconeDir"/scicone --input_breakpoints_file $chrBoundaries &> "$sciconeDir"/scicone.log
retVal=$?

# if breakpoint detection fails with brekapoints file, try again without it
if [ $retVal -ne 0 ] ; then
  echo "retrying without breakpoint file" >> "$sciconeDir"/scicone.log
  time $sciconeBin/breakpoint_detection --d_matrix_file "$sciconeDir"/"$outputBase".csv --n_bins $n_bins --n_cells $n_cells --postfix "$sciconeDir"/scicone &>> "$sciconeDir"/scicone.log
fi

retVal=$?
if [ $retVal -ne 0 ] ; then
  echo "SCICoNE failure, quitting" >> "$sciconeDir"/scicone.log
  exit 1
fi

# then call segment_counts.py
time python3 $sciconeBin/../scripts/segment_counts.py "$sciconeDir"/"$outputBase".csv "$sciconeDir"/scicone_segmented_region_sizes.txt &>> "$sciconeDir"/scicone.log

# then infer copy number tree
n_regions=$(wc -l < "$sciconeDir"/scicone_segmented_region_sizes.txt)
time $sciconeBin/inference --n_cells $n_cells --d_matrix_file "$sciconeDir"/"$outputBase"_segmented_counts.txt --region_sizes_file "$sciconeDir"/scicone_segmented_region_sizes.txt --n_regions $n_regions --postfix "$sciconeDir"/scicone &>> "$sciconeDir"/scicone.log

# then split scicone_inferred_cnvs.csv back into cell specific bed files
while read cnvs <&3 && read cellID <&4 ; do
  outFile="$sciconeDir""/""$cellID"".scicone.bed"
  paste $windows <(echo $cnvs | tr ',' '\n') > $outFile
done 3< "$sciconeDir"/scicone_inferred_cnvs.csv 4< "$tumorIDs"

