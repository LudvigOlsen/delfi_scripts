
<!-- README.md is generated from README.Rmd. Please edit that file -->

# delfi_scripts

This fork of `delfi_scripts` is for benchmarking against the DELFI
method (short/long fragment ratios in 5mb bins). While we have attempted
to stay as true to the original approach as possible, the following
changes have been made:

1)  We have made it runnable with the HG38 assembly.

2)  We have removed the summary statistics, as we only want the feature
    creation. We run the cross-validation in python to have the same
    setup for all benchmarks.

3)  We have ignored the z-scores (a meaasure of chromosome-arm-specific
    copy number changes). In the paper, these do not add significantly
    to the model AUC.

4)  The features are saved as a csv file, so we can read them into
    python more easily. The feature creation is performed directly in
    `04-5mb_bins.r` instead of `06-gbm_full.r`.

5)  Added/replaced command line arguments to use `optparse` and handle
    in- and output files in parent workflow.

6)  GC correction failed when the GC content (`bias` argument) had
    `NA`s. So we set them to the median value.

7)  Extraction of GC contents per fragment required insane amounts of
    RAM (\>512gb available for some files), so we replaced the original
    approach with the one used in DELFI Lucas (later paper).

8)  When summing across the 50 100kbp bins at a time, we ignore `NA`s
    instead of making the full 5mbp bin `NA`. This should reduce the
    number of 5mbp bins that are `NA` leading to removal of entire bins
    (which is problematic in a cross-dataset context when the test set
    has new NAs - would require imputation or similar).
