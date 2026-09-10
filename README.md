# nftide-placseq #

Nextflow pipeline for paired-end PLAC-seq / HiChIP data using HiC-Pro and FitHiChIP.

## Introduction ##

The pipeline merges FASTQs with the same sample name, uses cutadapt to remove Nextera adapters, and processes the trimmed reads with HiC-Pro. FitHiChIP then uses the HiC-Pro valid pairs, raw contact matrix, and bin intervals to call significant interactions.

If peak files are not provided, narrow peaks are inferred from HiC-Pro results using `PeakInferHiChIP.sh` and MACS2. Each sample is processed with both FitHiChIP backgrounds:

- **FitHiChIP_L**: `UseP2PBackgrnd=0`, using peak-to-peak and peak-to-nonpeak interactions for background modeling.
- **FitHiChIP_S**: `UseP2PBackgrnd=1`, using peak-to-peak interactions for background modeling.

The default reference is **hg38_HBV_A**, with MboI restriction fragments. HiC-Pro produces raw and ICE-normalized contact maps at **5 kb, 10 kb, 100 kb, and 1 Mb**. FitHiChIP uses the **5 kb raw matrix**.

## Software dependencies ##

Dependencies | Version / configuration
------------- | -------------
Nextflow | 25.10.2
Java | Compatible with the installed Nextflow; existing launcher uses the `hic` environment
HiC-Pro | 3.1.0
FitHiChIP | Manually installed by user
MACS2 | 2.2.9.1; resolved from the active conda environment
samtools | 1.16; resolved from the active conda environment for FitHiChIP
Bowtie2 | Required by HiC-Pro
bedtools, bgzip, tabix | Required by FitHiChIP
Python | With HiC-Pro dependencies and FitHiChIP's networkx dependency
R | With packages required by HiC-Pro and FitHiChIP, including optparse, data.table, ggplot2, fdrtool and GenomicRanges
bedToBigBed | Browser-output helper available in the pipeline's `bin/` directory

HiC-Pro uses the dependency paths in its installed `config-system.txt`, so install HiC-Pro **after** the installation of the conda enviroment, and provide dependency paths of the placseq conda environment. FitHiChIP and peak calling use tools from the active environment. MACS2 and samtools have no explicit executable overrides in this workflow.

## Installation ##

(1) Clone the repository and enter the pipeline directory. Replace `<repository-url>` with the repository URL.

```bash
git clone <repository-url> nftide-placseq
cd nftide-placseq
```

(2) Create and activate the `placseq` mamba environment from the pinned package export:

```bash
mamba env create --file placseq.yml
mamba activate placseq
```

(3) Restore the recorded R package library from `renv.lock`. Run these commands from the `nftide-placseq` directory after activating `placseq`.

```bash
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::restore(lockfile = "renv.lock", prompt = FALSE)'
```

`placseq.yml` recreates the pinned mamba/conda packages. `renv.lock` restores the R package versions into the repository-local `renv` library. The restore needs network access to CRAN and Bioconductor the first time it is run.

This workflow also requires HiC-Pro and FitHiChIP. Users needed to install them manually. Set installation paths of `hicpro` and `fithichip_dir` in `nextflow.config`.

(4) Prepare a Bowtie2 index, chromosome-size file, and restriction-fragment BED for the same reference. Modify `nextflow.config` accordingly. Current defaults for demonstration are:

```text
bowtie2_index: /data/xrz/ref/hg38_HBV_A_bowtie2/hg38_HBV_A
chrom_sizes: hg38_original_chrom_sizes.tsv
genome_fragment: hg38_HBV_A_MboI.bed
ligation_site: GATCGATC
```

The chromosome-size file must include the reference sequences to be represented in the contact maps. The Bowtie2 parameter is an **index prefix**, not just the containing folder. Prepare restriction enzyme cut sites following the instruction of HiC-Pro, and provide the correct logation sequences for HiC-Pro.

## Usage ##

(1) Prepare `samplesheet.csv`. The CSV file **must** contain three columns with these names:

`sample`: Name of the sequenced library and its output folder. FASTQs with the same sample name are merged before processing. Use letters, numbers, underscores, periods, or hyphens; start the name with a letter or number.  
`fastq_1`: Absolute path to the gzipped read 1 FASTQ.  
`fastq_2`: Absolute path to the matching gzipped read 2 FASTQ.

```csv
sample,fastq_1,fastq_2
sample1,/path/sample1_lane1_R1.fastq.gz,/path/sample1_lane1_R2.fastq.gz
sample1,/path/sample1_lane2_R1.fastq.gz,/path/sample1_lane2_R2.fastq.gz
sample2,/path/sample2_R1.fastq.gz,/path/sample2_R2.fastq.gz
```

Rows sharing a sample name are concatenated in CSV order. Input filenames must be distinct within each merge task because the workflow stages them by filename.

(2) Run the Nextflow pipeline:

```bash
nextflow run main.nf \
  --input_csv samplesheet.csv \
  -output-dir /data/xrz/PLAC/output \
  -with-report nf_plac_report.html \
  -with-timeline nf_plac_timeline.html \
  -with-trace nf_plac_trace.tsv \
  -resume -bg
```

If Nextflow and Java are not on PATH in `placseq`, use the existing installation explicitly:

```bash
mamba run -n placseq env JAVA_HOME=/home/xrz/miniforge3/envs/hic \
  /home/xrz/miniforge3/envs/hic/bin/nextflow run main.nf \
  --input_csv samplesheet.csv \
  -output-dir /data/xrz/PLAC/output \
  -resume
```

`-output-dir`: Output directory; default is `/data/xrz/PLAC/output`.  
`--input_csv`: Samplesheet; default is `samplesheet.csv`. Use `samplesheet_test.csv` for the included small inputs.  
`-resume`: Reuse completed tasks whose inputs and settings match the cache. Keep both `work/` and `.nextflow/`. Interrupted tasks may need to restart.  
`-bg`: Run Nextflow in the background.

All analysis settings are exposed in `nextflow.config` and can be overridden with `--parameter value`. The workflow generates per-task HiC-Pro and FitHiChIP configuration files; it does **not** read the legacy `config-hicpro.txt`.

Parameter | Default | Description
------------- | ------------- | -------------
`bowtie2_index` | hg38_HBV_A index prefix shown above | Shared Bowtie2 index
`chrom_sizes` | `hg38_original_chrom_sizes.tsv` | Chromosome names and lengths
`genome_fragment` | `hg38_HBV_A_MboI.bed` | Restriction fragments
`ligation_site` | `GATCGATC` | HiC-Pro ligation junction
`bin_size` | `5000 10000 100000 1000000` | HiC-Pro matrix resolutions in bp
`hicpro_cpus` | `16` | HiC-Pro task CPUs and generated N_CPU
`min_read_length` | `2` | Minimum cutadapt read length for each mate
`adapter_overlap` | `1` | Minimum adapter overlap
`min_mapq` | `10` | HiC-Pro mapping-quality threshold
`fithichip_bin_size` | `5000` | FitHiChIP BINSIZE; must occur in `bin_size`
`fithichip_circular_genome` | `0` | CircularGenome
`fithichip_int_type` | `3` | IntType: peak-to-all interactions
`fithichip_low_dist` | `20000` | Minimum loop distance in bp
`fithichip_upp_dist` | `2000000` | Maximum loop distance in bp
`fithichip_bias_type` | `1` | Coverage-bias correction
`fithichip_merge_int` | `1` | Merge filtering enabled
`fithichip_qvalue` | `0.05` | Loop-significance QVALUE threshold
`fithichip_cpus` | `16` | Nextflow CPU allocation and OMP_NUM_THREADS
`fithichip_peak_file` | Empty | Shared peak file or sample-to-path map
`fithichip_read_length` | `150` | Read length supplied to PeakInferHiChIP.sh
`fithichip_peak_genome` | `hs` | MACS2 effective-genome-size argument

Both backgrounds are always run: `UseP2PBackgrnd` is assigned to 0 and 1 by the workflow. `ValidPairs`, `Interval`, and `Matrix` come directly from the corresponding sample's HiC-Pro output. `Interval` and `Matrix` use the raw matrix resolution selected by `fithichip_bin_size`, not the ICE-normalized matrix.

By default, each process has `maxForks=2` and the executor is local. Different pipeline stages can run concurrently. Adjust `process.cpus`, `hicpro_cpus`, `fithichip_cpus`, and `maxForks` in `nextflow.config` for the available resources.

(3) Configure peak input.

Leave `fithichip_peak_file` empty to infer narrow peaks from HiC-Pro dangling-end, self-circle, re-ligation, and short-range valid pairs. The required interaction-class files are retained by the default `get_all_interaction_classes=1`.

The MACS2 options are exposed as `fithichip_macs2_options`:

```text
-q 0.05 --nomodel --shift 0 --extsize 200 --keep-dup all -B --SPMR --nolambda
```

Alternatively, provide a shared BED/narrowPeak file, optionally gzipped:

```bash
nextflow run main.nf --fithichip_peak_file /path/reference.narrowPeak -resume
```

For different peak files per sample, set a map in a Nextflow config:

```groovy
params.fithichip_peak_file = [
    GM12878_1: '/path/GM12878_1.narrowPeak',
    GM12878_2: '/path/GM12878_2.narrowPeak'
]
```

Samples missing from the map use automatic peak inference. Supplied BED files are copied to the task's standard filename `peaks.narrowPeak`; their contents are not converted into the ten-column narrowPeak format.

## Expected output ##

The pipeline creates a separate folder for each sample under `-output-dir`:

```text
output/<sample>/
├── fastqs/
│   ├── <sample>_merged_R1.fq.gz
│   ├── <sample>_merged_R2.fq.gz
│   └── <sample>_cutadapt.log
├── hicpro/
│   ├── config-hicpro.txt
│   ├── rawdata/
│   ├── bowtie_results/
│   ├── hic_results/
│   │   ├── data/<sample>/
│   │   ├── matrix/<sample>/
│   │   │   ├── raw/
│   │   │   └── iced/
│   │   ├── stats/<sample>/
│   │   └── pic/<sample>/
│   └── logs/<sample>/
├── FitHiChIP_peaks/
│   └── peaks/
│       ├── peaks.narrowPeak
│       ├── peak_calling.log                 # Automatic peak inference only
│       └── MACS2_ExtSize/                  # Automatic peak inference only
├── FitHiChIP_L/
│   ├── configuration_file
│   ├── run.log
│   ├── Summary_results_FitHiChIP.html
│   └── ...                                # Native FitHiChIP results
└── FitHiChIP_S/
    ├── configuration_file
    ├── run.log
    ├── Summary_results_FitHiChIP.html
    └── ...
```

The `fastqs/` folder publishes **merged FASTQs and the cutadapt log**. Trimmed reads are passed to HiC-Pro and are represented in its `rawdata/` folder; they are not separately published as cutadapt FASTQs in `fastqs/`.

HiC-Pro's main contact-level result is `hic_results/data/<sample>/<sample>.allValidPairs`, containing valid contacts after duplicate removal. Matrix results include the raw sparse `.matrix`, bin-coordinate `_abs.bed`, ICE-normalized matrix, and bias files at each configured resolution.

For FitHiChIP, open `Summary_results_FitHiChIP.html` in each L/S folder to locate the native output files. With the default settings, scored and significant interaction files are under:

```text
FitHiChIP_Peak2ALL_b5000_L20000_U2000000/
  P2PBckgr_<0-or-1>/Coverage_Bias/FitHiC_BiasCorr/
```

`FitHiChIP.interactions_FitHiC.bed` contains scored interactions.  
`FitHiChIP.interactions_FitHiC_Q0.05.bed` contains interactions passing the significance threshold.  
Merge-filtered and browser-compatible outputs are produced when applicable.

A completed analysis may find no significant loops; this is distinct from a runtime failure. Failed tasks may not publish their results: inspect `run.log` in the work directory printed by Nextflow, along with `.command.err` and `.nextflow.log`.

`main_hic_pro.nf` is the HiC-Pro-only backup. It uses the current `nextflow.config` but does not execute FitHiChIP.
