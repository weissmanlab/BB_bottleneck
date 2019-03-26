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

## Usage

usage: Bottleneck_size_estimation_approx.r [-h] [--file FILE]
                                           [--plot_bool PLOT_BOOL]
                                           [--var_calling_threshold VAR_CALLING_THRESHOLD]
                                           [--Nb_min NB_MIN] [--Nb_max NB_MAX]
                                           [--confidence_level CONFIDENCE_LEVEL]

optional arguments:
  -h, --help            show this help message and exit
  --file FILE           file containing variant frequencies
  --plot_bool PLOT_BOOL
                        determines whether pdf plot approx_plot.pdf is
                        produced or not
  --var_calling_threshold VAR_CALLING_THRESHOLD
                        variant calling threshold
  --Nb_min NB_MIN       Minimum bottleneck value considered
  --Nb_max NB_MAX       Maximum bottleneck value considered
  --confidence_level CONFIDENCE_LEVEL
                        Confidence level (determines bounds of confidence
                        interval)







usage: Bottleneck_size_estimation_exact.r [-h] [--file FILE]
                                          [--plot_bool PLOT_BOOL]
                                          [--var_calling_threshold VAR_CALLING_THRESHOLD]
                                          [--Nb_min NB_MIN] [--Nb_max NB_MAX]
                                          [--confidence_level CONFIDENCE_LEVEL]

optional arguments:
  -h, --help            show this help message and exit
  --file FILE           file containing variant frequencies and reads
  --plot_bool PLOT_BOOL
                        determines whether pdf plot exact_plot.pdf is produced
                        or not
  --var_calling_threshold VAR_CALLING_THRESHOLD
                        variant calling threshold
  --Nb_min NB_MIN       Minimum bottleneck value considered
  --Nb_max NB_MAX       Maximum bottleneck value considered
  --confidence_level CONFIDENCE_LEVEL
                        Confidence level (determines bounds of confidence
                        interval)






## Examples

For the approximate code run

  Rscript  Bottleneck_size_estimation_approx.r --file "example_data/donor_and_recipient_freqs.txt"  --plot_bool TRUE  --var_calling_threshold 0.03  --Nb_min 1 --Nb_max 200 --confidence_level .95

For the exact code run

  Rscript  Bottleneck_size_estimation_exact.r --file "example_data/donor_freqs_recip_freqs_and_reads.txt"  --plot_bool TRUE  --var_calling_threshold 0.03  --Nb_min 1 --Nb_max 200 --confidence_level .95

The six optional command line arguments for the codes are:

- file with lists of donor frequencies and recipient frequencies and reads*
- TRUE or FALSE (determines if pdf plot is produced)
- variant calling threshold
- minimum bottleneck size
- maximum bottleneck size
- confidence level

For the approximate code the first argument is a two column file, with the columns containing variant frequencies at each locus for the donor and recipient respectively.  For the exact code the first argument is a file with donor and recipient frequencies in the first two columns, total recipient reads in the third column, and variant recipient reads in the fourth column.  Files structured to run for the exact code can also be used with the approximate code.

If an argument is not given the default value will be assigned.

## Output
- Bottleneck Estimate
- Bounds of Confidence Interval
- Plots in pdf format (optional)
