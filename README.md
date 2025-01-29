# ONT Processing

ONT data were first processed with a modified version of nanoranger available from https://github.com/ForrestCKoch/nanoranger. The pipeline was modified so that barcodes must be an exact match to the whitelist in order to be considered.

The relevant job submission script for performing nanoranger preprocessing is `scripts/run_nanoranger_v2.sh`

Mutation calling is then performed by a custom script `scripts/get_mutations_cb_only_with_quality_threshold_v2.py` with relative job submission script `scripts/run_get_mutations_cb_only_with_quality_threshold.py`. Versions without quality thresholds are also available (see `scripts/get_mutations_cb_only_with_qualities.py`). Note that the other variants of `get_mutations...py` have a bug relating to read direction and should not be used.

In brief, the mutation calling code assigns a label of wild-type (WT), mutant-type (MT), or unknown (UN), to each read which covers a target mutation. UN indicates that the base pair sequence matches neither the MT or WT sequence (note that this is needed due to the high error rate of ONT). Total counts for each mutation category are calcualted per-barcode and written to `results/mutations-called-cb-only-with-quality-threshold_v2/mutations-called_*_threshodl-*.csv` 



