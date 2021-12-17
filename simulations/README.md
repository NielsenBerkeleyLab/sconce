# sconce_sim
This file describes the simulation program for SCONCE.

## Compiling
To compile the program, run `make`. All dependencies are included here.

## Generating simulated datasets
To run it, run `./sconce_sim <infile> <paramfile> [optional random seed]`. This will create the following files in the current working directory:
- `simu_cancer_cell_*` (observed read depth for cancer cells)
- `simu_healthy_cell_*` (observed read depth for healthy cells)
- `true_cancer_cell_*` (true ploidy for cancer cells)
- `true_healthy_cell_*` (true ploidy for healthy cells)

To convert the `simu_*_cell*` files (which have dummy genomic coordinates) into hg19 (or similar) coordinates, first create a bin file of your reference genome (format should be tab separated `<chr>\t<start>\t<end>`). Make sure the paramfile specifies the same number of bins in your reference genome, then run
```
paste /path/to/reference/bin/file <(awk '{print $4}' /path/to/simu/file) | awk 'BEGIN{OFS="\t"} {if($4 == "") {$4=0} print}' > /path/to/output
```

## Parameter file formats
Any lines in the infile and paramfile after the arguments are ignored (ie good for comments).

The format of the `<infile>` varies between the line segment and binned models.

The line segment model `<infile>` format is:
```
<flag for linesegment (0) model> <genome length> <# tumor cells> <flag for coalescence simulations (0)> <# healthy cells>
<rate of deletion> <rate of amplification> <mean deletion length> <mean amplification length>
<length of edge leading to root of tree = length of simulation time for single cell>
<the parameter of the growth rate model from Hudson and Slatkin, alpha is 0.0. A value of 0.0 corresponds to neutral coalescence simulations>
```

The binned model `<infile>` format is:
```
<flag for binned (1) model> <genome length> <# tumor cells> <flag for coalescence simulations (0)> <# healthy cells>
<maximum ploidy, k> <geometric distribution parameter, probability of stopping>
<transition matrix, must be k rows by k cols, and each row must sum to 0>
<...>
<length of edge leading to root of tree = length of simulation time for single cell>
<the parameter of the growth rate model from Hudson and Slatkin, alpha is 0.0. A value of 0.0 corresponds to neutral coalescence simulations>
```

All published simulation sets use the same `<paramfile>` format. The `<paramfile>` format is:
```
<# bins to divide the genome into> <parameter r of the negative binomial (must be integer, converges to Poisson as r goes to infinity> <total expected number of reads (note that the observed number is random)> <flag for even coverage in expectation before accounting for CNAs>
```

## Example parameter files
Parameter files for the SCONCE paper dataset are included in the inputFiles directory for reproducibility. All of the simulation sets use the same paramfile.txt. Sets A-D are under the line segment model, and E-K are under the binned model.

Note that set B has a whole genome duplication before any CNAs. This requires calling `makeDiploidGenomeNoGapsDuplicated` instead of `makeDiploidGenomeNoGaps` in `sconce_sim.c`. To do this, uncomment `makeDiploidGenomeNoGapsDuplicated` (line 849) and comment out `makeDiploidGenomeNoGaps` (line 848) and recompile.

The included test files are from parameter set A, with seed 0. The tumor cell is cell 0.

