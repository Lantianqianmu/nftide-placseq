# Pre-edit pipeline backup: 1.0.1

Git branch `1.0.1` (commit `4f110f680b76a1b4cd0c6cb800a27bb25f123f4c`) preserves the pipeline code and analysis configuration immediately before the samplesheet peak-file redesign. The historical files recorded in `BACKUP_MANIFEST.json` at that commit are byte-for-byte copies of the working local pipeline, including the matched 1D bigWig implementation, broad/narrow MACS2 selection, replicate-name fix, and 10 kb FitHiChIP configuration.

The original machine-specific input samplesheet and reference data are preserved separately in a local snapshot. This public branch keeps the repository's existing example samplesheet. To reproduce an analysis, supply the original samplesheet, FASTQs, matching reference assets, and external HiC-Pro/FitHiChIP/tool installations, and adapt installation/data paths if running on another machine. The conda export and R lockfile record the software environment.

The manifest describes the historical backup only; files in a later working tree or release can differ.
