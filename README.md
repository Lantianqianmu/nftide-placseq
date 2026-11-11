# nftide-placseq

Nextflow pipeline for paired-end PLAC-seq / HiChIP data using HiC-Pro and FitHiChIP, with configurable loop-anchor peaks and per-sample 1D bigWig tracks.

## Workflow

1. Merge FASTQ rows with the same **sample** ID.
2. Trim adapters with cutadapt and process each sample independently with HiC-Pro.
3. Select immunoprecipitation-enriched reads from each sample, run MACS2, and generate its 1D track.
4. Assign loop-anchor peaks using the samplesheet's **peak_file** column. Missing shared paths trigger pooled peak discovery; blank fields use individual sample peaks.
5. Run FitHiChIP for each sample with both backgrounds:
   - **L:** `UseP2PBackgrnd=0`, peak-to-peak and peak-to-nonpeak background.
   - **S:** `UseP2PBackgrnd=1`, peak-to-peak background.

The default FitHiChIP resolution is **10 kb**, with a **20 kb–2 Mb** distance window. HiC-Pro generates raw and ICE-normalized matrices at the configured resolutions. The HiC-Pro-only workflow is maintained separately in [nftide-hicpro](https://github.com/Lantianqianmu/nftide-hicpro).

## Dependencies and installation

| Software | Configuration |
| --- | --- |
| Nextflow | 25.10.2, DSL2 and workflow output publishing |
| Java | 17 |
| HiC-Pro | 3.1.0, installed separately |
| FitHiChIP | 11.0, installed separately |
| MACS2 | 2.2.9.1 |
| cutadapt | 4.9 |
| Bowtie2 / samtools | 2.5.4 / 1.16.1 |
| Python | Python 3; peak preparation uses the standard library |
| R | With the packages required by HiC-Pro and FitHiChIP |
| bedGraphToBigWig | UCSC executable, configured explicitly or detected in the active environment |
| bedToBigBed | Browser-output helper distributed under `bin/` |

Clone the repository and select the pipeline branch/release you intend to run:

```bash
git clone https://github.com/Lantianqianmu/nftide-placseq.git
cd nftide-placseq
mamba env create --file placseq.yml
mamba activate placseq
```

`conda` can replace `mamba`. The environment export records Linux packages. Install [HiC-Pro](https://nservant.github.io/HiC-Pro/INSTALL.html) after preparing the environment so that its installed `config-system.txt` points to the intended dependencies. Install [FitHiChIP](https://github.com/ay-lab/FitHiChIP) separately.

The R package library is recorded in `renv.lock`. When restoring that library, follow the renv setup needed by your installation; a repository-local renv restore does not automatically make that library visible to R scripts launched from the separate FitHiChIP installation. Ensure the required packages are available to those Rscript processes, for example through the active environment's R library or an explicitly configured `R_LIBS_USER`.

Activate the intended environment before launching Nextflow. The workflow uses executable tools on `PATH`, plus the installation paths configured below.

## Reference preparation and configuration

Prepare a Bowtie2 index, chromosome sizes, and restriction-fragment BED from the **same reference genome**. The chromosome-size file contains sequence names and lengths; sequence names must agree across all assets and the supplied peaks.

For example, for MboI:

```bash
mkdir -p reference/bowtie2
bowtie2-build /path/to/genome.fa reference/bowtie2/genome
samtools faidx /path/to/genome.fa
cut -f1,2 /path/to/genome.fa.fai > reference/genome.chrom.sizes
python /path/to/HiC-Pro/bin/utils/digest_genome.py \
  -r '^GATC' -o reference/genome_MboI.bed /path/to/genome.fa
```

Configure paths and analysis settings in `nextflow.config` or an additional `site.config`:

```groovy
params {
    input_csv = '/path/to/samplesheet.csv'
    bowtie2_index = '/path/to/reference/bowtie2/genome'
    chrom_sizes = '/path/to/reference/genome.chrom.sizes'
    genome_fragment = '/path/to/reference/genome_MboI.bed'
    ligation_site = 'GATCGATC'
    hicpro = '/path/to/HiC-Pro/bin/HiC-Pro'
    fithichip_dir = '/path/to/FitHiChIP'
    bedgraph_to_bigwig = '/path/to/bedGraphToBigWig'
    bin_size = '5000 10000 100000 1000000'
    fithichip_bin_size = 10000
}

outputDir = '/path/to/results'
workflow.output.mode = 'copy'
cleanup = false
```

The Bowtie2 path is an **index prefix** with all six `.bt2` or all six `.bt2l` files. Change the restriction fragments and ligation junction together for other enzymes. Adapter sequences must match the library preparation. Reference assets and sequencing data are supplied by the user.

### Main parameters

| Parameter | Default / use |
| --- | --- |
| `input_csv` | Project `samplesheet.csv` |
| `bowtie2_index`, `chrom_sizes`, `genome_fragment` | Set for the reference genome |
| `hicpro`, `fithichip_dir` | Set for the software installations |
| `bedgraph_to_bigwig` | Explicit executable path; set to `''` for active-conda/PATH detection |
| `bin_size` | `5000 10000 100000 1000000` bp; spaces or commas accepted |
| `fithichip_bin_size` | `10000` bp; must also occur in `bin_size` |
| `hicpro_cpus`, `fithichip_cpus` | `16` each |
| `adapter_r1`, `adapter_r2` | `CTGTCTCTTATACACATCT`; change for the library |
| `min_read_length`, `adapter_overlap` | `2` bp per mate, `1` bp overlap |
| `min_mapq` | `10` |
| `rm_dup`, `get_all_interaction_classes` | Must be `1` for this workflow |
| `fithichip_read_length` | `150` bp; used to reconstruct BED read ends for MACS2 |
| `fithichip_peak_genome` | `hs`; set MACS2's effective genome size for the reference |
| `fithichip_macs2_options` | `-q 0.05 --extsize 147 --nomodel` |
| `fithichip_macs2_broad_cutoff` | `0.1` |
| `fithichip_int_type` | `3`, peak-to-all |
| `fithichip_low_dist`, `fithichip_upp_dist` | `20000`, `2000000` bp |
| `fithichip_bias_type`, `fithichip_merge_int` | `1`, `1` |
| `fithichip_qvalue` | `0.05`, loop-significance threshold |

All remaining HiC-Pro and FitHiChIP analysis settings are exposed in `nextflow.config`. The workflow generates each task's native configuration from these parameters. The legacy `config-hicpro.txt` is not read by the workflow.

The default executor is `local`, with `process.cpus=16` and `process.maxForks=4`. `PLAC_SIGNAL` and `SHARED_PEAKS` request two CPUs each; HiC-Pro and FitHiChIP use their dedicated CPU parameters. `maxForks` applies per process, so stages can run concurrently. Adapt executor/resources to the available machine or scheduler.

## Samplesheet: four required columns

```csv
sample,fastq_1,fastq_2,peak_file
control-27ac-1,/reads/control_rep1_lane1_R1.fq.gz,/reads/control_rep1_lane1_R2.fq.gz,
control-27ac-1,/reads/control_rep1_lane2_R1.fq.gz,/reads/control_rep1_lane2_R2.fq.gz,
control-27ac-2,/reads/control_rep2_R1.fq.gz,/reads/control_rep2_R2.fq.gz,
case-27ac-1,/reads/case_rep1_R1.fq.gz,/reads/case_rep1_R2.fq.gz,reference/chip_27ac.bed.gz
case-27ac-2,/reads/case_rep2_R1.fq.gz,/reads/case_rep2_R2.fq.gz,reference/chip_27ac.bed.gz
```

- **sample:** library/sample ID and output-folder name. Start with a letter or number; letters, numbers, underscores, periods, and hyphens are accepted.
- **fastq_1 / fastq_2:** matching gzipped paired-end FASTQs. Absolute paths are recommended. Relative FASTQ paths are resolved against the samplesheet directory.
- **peak_file:** optional peak-file path. **Keep the column and trailing comma when the value is blank.** Relative peak paths are resolved against the samplesheet directory and normalized before grouping.

Rows with the same sample ID are merged in CSV order **before alignment and duplicate removal**. These should be sequencing lanes of the same library. Give separate biological samples different IDs. All rows for one sample must repeat the same `peak_file` value, including blank values. Repeated FASTQ basenames in different lane folders are supported through separate staging directories.

### How the peak_file column works

| `peak_file` value at workflow start | FitHiChIP loop-anchor peaks |
| --- | --- |
| Path to an existing, readable, nonempty file | Use that supplied peak set for the sample |
| Nonempty path that does not exist | Pool peak-calling BED reads from all samples with the same normalized path, call one peak set, and use it for every member |
| Blank field | Use peaks inferred separately from that sample's reads |

**The condition is checked at workflow start.** An existing directory, empty file, or unreadable file is an error rather than a request for inference. A missing path used by only one sample produces a one-sample peak group.

#### A. Use an existing peak file

```csv
sample,fastq_1,fastq_2,peak_file
control-27ac-1,/reads/control_R1.fq.gz,/reads/control_R2.fq.gz,/peaks/chip_27ac.narrowPeak
case-27ac-1,/reads/case_R1.fq.gz,/reads/case_R2.fq.gz,/peaks/chip_27ac.narrowPeak
```

If `/peaks/chip_27ac.narrowPeak` exists, both samples use that peak set. Their FASTQs, HiC-Pro contacts, FitHiChIP calls, and 1D tracks remain separate.

Plain or gzip-compressed BED, narrowPeak, and broadPeak files are accepted. BED intervals are **zero-based, half-open**. Blank lines and `#`, `track`, and `browser` headers are ignored. Peak preparation retains additional columns, clips intervals to the chromosome bounds, and removes intervals on sequences absent from `chrom_sizes`; `validation.json` reports these actions. The original input file is unchanged. A file with no usable intervals fails clearly.

#### B. Infer one shared peak set

```csv
sample,fastq_1,fastq_2,peak_file
control-27ac-1,/reads/control_R1.fq.gz,/reads/control_R2.fq.gz,peaks/common_27ac.bed.gz
case-27ac-1,/reads/case_R1.fq.gz,/reads/case_R2.fq.gz,peaks/common_27ac.bed.gz
```

If `peaks/common_27ac.bed.gz` does not exist, the workflow concatenates the two samples' selected MACS2 input BEDs, sorts them, and calls MACS2 once for the group. Samples are deduplicated individually before pooling; MACS2's configured read-tag duplicate filtering still applies to the pooled input.

The shared result is published both under `results/shared_peaks/group_<path-hash>/peaks/` and at the requested path, creating its parent directories. A requested `.gz` path receives gzip-compressed BED content. Existing destination files are not overwritten by group publication. Use a writable shared filesystem path visible to the tasks.

On a later run, the generated path exists and is handled as a supplied peak file. To request a new pooled call after changing group membership, inputs, or settings, use a new missing output path or deliberately remove the old generated file. The peak-file status is excluded from alignment cache keys, so this transition does not itself force HiC-Pro to rerun.

Pooling here is **only for peak discovery**. The original sample contact matrices, valid pairs, and FitHiChIP analyses remain independent. Share a path only across compatible samples of the same assay target/reference. Groups mixing broad and narrow MACS2 target classes are rejected before FASTQ processing.

#### C. Infer peaks separately

```csv
sample,fastq_1,fastq_2,peak_file
control-27ac-1,/reads/control_R1.fq.gz,/reads/control_R2.fq.gz,
case-27ac-1,/reads/case_R1.fq.gz,/reads/case_R2.fq.gz,
```

Each sample receives its own inferred peak set. Multiple FASTQ rows sharing one sample ID contribute to that sample's inference after their normal merge/alignment/deduplication.

The former `fithichip_peak_file` configuration parameter is replaced by this CSV column. Move any nonempty legacy assignment into the samplesheet; nonempty legacy settings produce an explicit migration error.

## Which reads are used for peaks and 1D signal?

For each sample, the workflow selects:

1. **Dangling-end pairs** (`*.DEPairs`).
2. **Self-circle pairs** (`*.SCPairs`).
3. **Re-ligation pairs** (`*.REPairs`).
4. **Duplicate-removed cis valid pairs with separation strictly below 1,000 bp**, read directly from `<sample>.allValidPairs`.

The DE/SC/RE files can contain duplicates. The combined selected pairs are deduplicated using the two mapped chromosome/position coordinates, with pair orientation canonicalized, within each sample. The workflow does not fall back to pre-deduplication `*.validPairs` files. Both read ends are converted to strand-aware BED intervals using `fithichip_read_length` and clipped to the reference.

Long-range valid contacts and trans valid contacts do not contribute to this inferred ChIP-like 1D signal. They remain available in the original HiC-Pro valid-pair/contact inputs for FitHiChIP. `signal/read_selection.tsv` records the input source, eligible pairs, selected pairs after deduplication, and retained BED read ends.

### Narrow and broad MACS2 modes

The `fithichip_peak_types` map controls the target class:

```groovy
params.fithichip_peak_types = [
    '27ac': 'narrow',
    '9ac': 'narrow',
    '27me3': 'broad',
    '9me3': 'broad',
    'pol2': 'narrow',
    'custom_sample': 'narrow'
]
```

#### How a sample name selects its peak type

The lookup uses the **`sample` column of the samplesheet**. Sample names, map keys, and `broad`/`narrow` values are matched case-insensitively. For each sample:

1. **Look for an exact full sample-name key first.** For example, `'my_library_1': 'broad'` applies to the sample ID `my_library_1`. This exact entry takes priority over a mark-based entry.
2. **If there is no exact match, derive a mark key from the sample name:**
   - Remove one final numeric replicate suffix consisting of a separator and digits, such as `-1`, `_2`, or `.3`.
   - Split the remaining name on hyphens, underscores, or periods, and take the final token.
   - Remove an optional `H3K` prefix from that token. Use the resulting key, such as `27ac` or `27me3`, to look up the map.
3. **Require a matching `broad` or `narrow` value.** If neither lookup succeeds, the workflow stops with `Cannot determine MACS2 peak type ...` and asks for a sample-name or mark entry.

For example, `control_H3K27me3_2` becomes `control_h3k27me3_2`, then `control_h3k27me3`, then `h3k27me3`, and finally the lookup key `27me3`.

With the example map above:

| Samplesheet `sample` | Matching key | Result |
| --- | --- | --- |
| `control-27ac-1` | `27ac` | narrow |
| `control_H3K27ac_2` | `27ac` | narrow |
| `control.H3K9me3.3` | `9me3` | broad |
| `control-pol2` | `pol2` | narrow |
| `custom_sample` | Exact full name `custom_sample` | narrow |
| `custom_sample-1` | Fallback token `sample` has no entry | Error; add the full sample ID |
| `control-27ac-rep1` | Fallback token `rep1` has no entry | Error; add the full sample ID or use a numeric suffix |

#### What `custom_sample` means

**`custom_sample` is a literal example sample ID.** It is not a wildcard, a prefix rule, or a default for all custom libraries. The entry `'custom_sample': 'narrow'` matches `custom_sample` exactly (ignoring case). It does not automatically match `custom_sample-1`, `custom_sample-2`, or `custom_sample_rep1`.

For custom names without a recognized mark token, replace the example key with each complete sample ID. Keep the mark entries for any other samples that use them:

```groovy
params.fithichip_peak_types = [
    '27ac': 'narrow',
    '9ac': 'narrow',
    '27me3': 'broad',
    '9me3': 'broad',
    'pol2': 'narrow',
    'custom_sample-1': 'narrow',
    'custom_sample-2': 'narrow',
    'repressive_library_rep1': 'broad'
]
```

Alternatively, sample IDs such as `custom-27ac-1` and `custom-27ac-2` both use the `27ac` entry automatically. Only numeric replicate suffixes are removed; suffixes such as `rep1`, `r1`, or `lane1` require an explicit full-name entry when they occupy the final token. An exact full-name key always takes priority: adding `'control-27ac-1': 'broad'` would override the `27ac` narrow setting for that one sample.

The mapping is required even when the samplesheet supplies an existing peak file, because per-sample MACS2 signal generation is retained. Samples sharing a missing peak-file path must resolve to the same broad/narrow class.

#### MACS2 options

Shared MACS2 options default to:

```text
-q 0.05 --extsize 147 --nomodel
```

The workflow adds `--broad --broad-cutoff 0.1` for broad targets and adds `-B --SPMR` for pileup output. Configure broad mode through the map and broad cutoff parameter, rather than adding `--broad` to the shared options. `--call-summits` is incompatible with broad targets. MACS2's own read-tag filtering applies in addition to pair-level deduplication; set its options deliberately if changing that behavior.

### 1D bigWig tracks

**Every sample gets its own track in all three peak-file cases:**

```text
results/<sample>/bw/<sample>.plac_1d.cpm.bw
```

MACS2 runs per sample on the selected read-end BED to produce the treatment pileup with `-B --SPMR`. This is signal per million MACS2-retained tags. The pileup is clipped, sorted, and converted with `bedGraphToBigWig`. It is an inferred PLAC-seq/HiChIP 1D signal, rather than an independent ChIP-seq measurement.

For blank `peak_file` values, the sample peak set and track use the same MACS2 run. For a missing shared path, loop anchors use the pooled MACS2 run while each track uses its own sample's selected reads with the same settings. With externally supplied peaks, the supplied file defines loop anchors and the track still represents the sample's inferred MACS2 signal. Per-sample inferred peaks are retained for QC in all cases; they are not substituted for supplied/shared anchors.

Converter selection is:

1. Use the executable path in `bedgraph_to_bigwig`, if nonempty.
2. Otherwise use `$CONDA_PREFIX/bin/bedGraphToBigWig`, if executable.
3. Otherwise search `PATH`.

An explicitly configured invalid path fails; it is not silently replaced by another executable.

## Run

```bash
nextflow run main.nf \
  -c site.config \
  --input_csv /path/to/samplesheet.csv \
  -output-dir /path/to/results \
  -with-report nf_plac_report.html \
  -with-timeline nf_plac_timeline.html \
  -with-trace nf_plac_trace.tsv \
  -resume
```

Add `-bg` for a background run. `-preview` checks compilation and configuration without executing analysis tasks; valid reference settings are still needed. Keep `work/` and `.nextflow/` for resume. Intermediate workflow outputs explicitly use copy mode so publishing cannot move files away from downstream tasks; `workflow.output.mode='copy'` remains the recommended default.

## Expected output

```text
results/
├── <sample>/
│   ├── fastqs/                         # Merged FASTQs and cutadapt log
│   ├── hicpro/
│   │   ├── config-hicpro.txt
│   │   └── hic_results/
│   │       ├── data/<sample>/           # .allValidPairs and interaction classes
│   │       └── matrix/<sample>/         # raw/ and iced/ at configured resolutions
│   ├── bw/<sample>.plac_1d.cpm.bw
│   ├── PLAC_signal/signal/
│   │   ├── MACS2_input.bed
│   │   ├── inferred_peaks.bed
│   │   ├── MACS2_ExtSize/
│   │   ├── peak_calling.log
│   │   ├── read_selection.tsv
│   │   └── plac_1d.cpm.bedGraph
│   ├── FitHiChIP_peaks/peaks/
│   │   ├── peaks.bed                   # Actual anchors passed to FitHiChIP
│   │   ├── source.json                 # sample / provided / shared and group members
│   │   └── validation.json
│   ├── FitHiChIP_L/
│   │   ├── configuration_file
│   │   ├── run.log
│   │   ├── Summary_results_FitHiChIP.html
│   │   └── ...
│   └── FitHiChIP_S/                    # Same structure, strict background
└── shared_peaks/group_<path-hash>/peaks/ # Only for missing nonempty paths
    ├── MACS2_input.bed                 # Pooled selected read ends
    ├── peaks.bed
    ├── group.json
    ├── validation.json
    ├── MACS2_ExtSize/
    └── peak_calling.log
```

FitHiChIP uses each sample's **raw** matrix, matching `_abs.bed`, deduplicated `.allValidPairs`, and assigned `FitHiChIP_peaks/peaks/peaks.bed`. ICE matrices are retained for visualization. Scored/significant loops appear in native FitHiChIP subfolders, for example:

```text
FitHiChIP_L/FitHiChIP_Peak2ALL_b10000_L20000_U2000000/
  P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/
```

`FitHiChIP.interactions_FitHiC.bed` is the full scored contact table. `FitHiChIP.interactions_FitHiC_Q0.05.bed` contains significant contacts. Keep full scored tables and separate biological samples for differential analysis. A successful analysis can have zero significant loops; runtime errors and missing scoring outputs still fail.

## Tests and troubleshooting

After activating the environment:

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
python tests/run_smoke_tests.py --bedgraph-to-bigwig /path/to/bedGraphToBigWig
```

The smoke test creates small synthetic files and runs real cutadapt, MACS2, and bigWig conversion. HiC-Pro and FitHiChIP are explicit test doubles in this routing test; it does not validate their alignment or statistical algorithms. It checks all peak modes, broad/narrow behavior, gzip, duplicate selection, the 999/1,000 bp boundary, separate sample tracks, pooled-read counts, output publication, and L/S input handoffs. It also tests invalid input plans, resumes after shared paths become supplied inputs, and restores a missing requested peak file from cache. Run it in a fresh temporary directory when setting `--work-dir`.

For failures, inspect `.command.err`, `.command.sh`, and `.command.out` in the reported task work directory, plus `.nextflow.log`. For peaks, inspect the MACS2 log, `read_selection.tsv`, `source.json`, and `validation.json`. Empty/unusable supplied peaks, inconsistent peak values across lanes, unknown target types, and mixed broad/narrow shared groups fail explicitly. Use a new missing path when deliberately regenerating a shared peak group.

## References

- [HiC-Pro documentation](https://nservant.github.io/HiC-Pro/)
- [FitHiChIP documentation](https://ay-lab.github.io/FitHiChIP/)
- [FitHiChIP study](https://doi.org/10.1038/s41467-019-11950-y)
- [MACS2](https://github.com/macs3-project/MACS)
- [Nextflow documentation](https://www.nextflow.io/docs/stable/)
