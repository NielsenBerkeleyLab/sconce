# SCONCE

This directory is for the program SCONCE (Single Cell cOpy Numbers in CancEr).

## Dependencies
SCONCE is written in C++11, and requires
- GNU make (tested on v4.1)
- g++ (tested on 7.5.0)
- BOOST (tested on v1.65.1)
- GSL (tested on v2.4)

Additional [R scripts](scripts/) require
- R (tested on v3.6.3)
- ggplot2 (tested on v3.2.1)
- reshape2 (tested on v1.4.3)
- cowplot (tested on v1.1.1)
- scales (tested on v1.1.0)

SCONCE was developed and tested on Ubuntu 18.04.

## Installation instructions
1. Clone this repo:
```
git clone git@github.com:NielsenBerkeleyLab/sconce.git
```
2. Run `make`. This will build intermediates into the `build/` directory and create an executable named `sconce`.


## Example test run
To ensure `sconce` was built correctly, we include some test files. Run the following:
```
time ./sconce --diploid test/test_healthy_avg.bed --tumor test/test_cancer_cell.bed --meanVarCoefFile test/test.meanVar --outputBase test/output_k5 --maxKploid 5 > test/output_k5.log 2> test/output_k5.err
```
Your output (with the exception of timing information) should match the provided `test/ref_output_k5*` files.


## Brief parameter descriptions
- `--diploid` This file should be the averaged read depth across all diploid cells. It can be generated using [scripts/avgDiploid.R](scripts/avgDiploid.R).
- `--tumor` This file should be the per window read depth of a tumor cell to analyze.
- `--meanVarCoefFile` This file should define the coefficients for the relationship between the mean and variance of the negative binomial distribution used for emission probabilities. It should be generated using [scripts/fitMeanVarRlnshp.R](scripts/fitMeanVarRlnshp.R).
- `--outputBase` This gives the basename for all [output files](#output-files).
- `--maxKploid` This gives the maximum allowed ploidy (recommended `k=10`).


## Input files
Averaged diploid read depth files should be tab separated, with columns `<chr>\t<start>\t<end>\t<meanReadDepth>\t<varianceOfReadDepth>`. They should be generated using [scripts/avgDiploid.R](scripts/avgDiploid.R). For example:
```
chr1	0	250000	325.87	2014.1131
chr1	250000	500000	313.32	2017.6376
chr1	500000	750000	314.41	2165.7219
chr1	750000	1000000	330.23	1751.4371
chr1	1000000	1250000	313.69	2156.5539
chr1	1250000	1500000	321.75	2743.3475
chr1	1500000	1750000	327.08	2228.4736
chr1	1750000	2000000	318.79	2606.3459
chr1	2000000	2250000	326.76	2733.8824
chr1	2250000	2500000	319.27	3082.6371
```

Tumor read depth files should be tab separated, with columns `<chr>\t<start>\t<end>\t<readDepth>`. See [simulations/README.md](simulations/README.md) for how to generate simulations with this format. For real data, a tool like [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) can be used to create this file from a bam file. For example:
```
chr1	0	250000	699
chr1	250000	500000	804
chr1	500000	750000	627
chr1	750000	1000000	701
chr1	1000000	1250000	521
chr1	1250000	1500000	685
chr1	1500000	1750000	616
chr1	1750000	2000000	583
chr1	2000000	2250000	736
chr1	2250000	2500000	634
```

Mean and variance coefficient files should have one parameter and value pair per line. They should be generated using [scripts/fitMeanVarRlnshp.R](scripts/fitMeanVarRlnshp.R). For example:
```
intercept=10.80433928393
slope=1.17232224213
poly2=0.01918299582
```


## Output files
SCONCE will create the following files automatically:
- `<output>.hmm` This file contains the state of the HMM after the Baum Welch step and the state of the HMM after the BFGS step.
- `<output>.bed` This file contains the copy number calls in tab separated bed format (`<chr>\t<start>\t<end>\t<ploidy>`)
- `<output>.viterbiDecoded`. This tab separated files contains the copy number calls with more detail. The columns are `<coord>` (chr:start-end), `<diploid_mean>`, `<diploid_variance>`, `<tumor>` (observed tumor read count), `<ploidy_0-<k>>` (ploidy call in the range of 0 to `<k>`)

SCONCE also prints log messages to stdout and error messages to stderr.
If the `--verbose` flag is used, debugging statements will be printed to stderr.

Genome traces can be plotted using the included [scripts/plotGenomeTrace.R](scripts/plotGenomeTrace.R) script. The arguments for this script are:
- Arg 1: /path/to/healthy/average/bed/file
- Arg 2: /path/to/observed/tumor/read/depths/bed/file
- Arg 3: /path/to/sconce/output/bed/file
- Arg 4: /path/to/output/plot
- Arg 5: quoted text for the plot title (ex sample name)
- Arg 6: [optional] /path/to/ground/truth/ploidy/bed/file (ie for simulations)

Using the test files from before:
```
Rscript scripts/plotGenomeTrace.R test/test_healthy_avg.bed test/test_cancer_cell.bed test/ref_output_k5.bed test/ref_plot_k5.png "Reference Genome Trace for SCONCE (k=5)" test/true_cancer_cell.bed
```
[The output plot](test/ref_plot_k5.png) is included for reference.


## Simulations
To compile and run the simulation program, see [simulations/README.md](simulations/README.md).

