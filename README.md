# BB Bottleneck Estimator
## Description
This is the implementation for [transmission bottleneck estimation based on beta-binomial sampling](https://www.biorxiv.org/content/10.1101/101790v1).  

## Requirements
- R 3.5.2+
- rmutil
- argparse
- mdatools
- Rscript

Be sure to add the Rscript path to your environment variables or include the path when calling Rscript from the command line.  


## Examples

For the approximate code run

  Rscript --vanilla Bottleneck_size_estimation_approx.r "example_data/donor_and_recipient_freqs.txt" TRUE 0.03 1 200 .95

For the exact code run

  Rscript --vanilla Bottleneck_size_estimation_exact.r "example_data/donor_freqs_recip_freqs_and_reads.txt" TRUE 0.03 1 200 .95

The six command line arguments for the codes are:

- file with lists of donor frequencies and recipient frequencies and reads*
- TRUE or FALSE (determines if pdf plot is produced)
- variant calling threshold
- minimum bottleneck size
- maximum bottleneck size
- confidence level

For the approximate code the first argument is a two column file, with the columns containing variant frequencies at each locus for the donor and recipient respectively.  For the exact code the first argument is a file with donor and recipient frequencies in the first two columns, total recipient reads in the third column, and variant recipient reads in the fourth column.  Files structured to run for the exact code can also be used with the approximate code.

## Output
- Bottleneck Estimate
- Bounds of Confidence Interval
- Plots in pdf format (optional)
