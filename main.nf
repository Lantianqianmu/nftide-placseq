#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def resolvePeakType(String sample, def peakTypes) {
    if (!(peakTypes instanceof Map))
        throw new IllegalArgumentException('fithichip_peak_types must map sample names or mark suffixes to broad or narrow')
    def types = peakTypes.collectEntries { key, value ->
        [(key.toString().toLowerCase()): value.toString().toLowerCase()]
    }
    def sample_key = sample.toLowerCase()
    // A separate library may end in a numeric replicate ID (e.g. dC-27ac-1).
    // Ignore that suffix only for mark lookup; retain the full sample ID and overrides.
    def mark_sample = sample_key.replaceFirst(/[_.-][0-9]+$/, '')
    def mark = mark_sample.tokenize('_.-').last().replaceFirst(/^h3k/, '')
    def peak_type = types[sample_key] ?: types[mark]
    if (!(peak_type in ['broad', 'narrow']))
        throw new IllegalArgumentException("Cannot determine MACS2 peak type for ${sample}; add its sample name or mark suffix to fithichip_peak_types as broad or narrow")
    return peak_type
}

process MERGE_FQ {
    tag "Merging fastq files of ${meta.id}..."

    input:
    tuple val(meta), path(r1s), path(r2s)

    output:
    tuple val(meta), path("*_merged_R1.fq.gz"), path("*_merged_R2.fq.gz"), emit: merged_fq

    script:
    """
    cat ${r1s.join(' ')} > ${meta.id}_merged_R1.fq.gz
    cat ${r2s.join(' ')} > ${meta.id}_merged_R2.fq.gz
    """
}

process CUTADAPT {
    tag "cutadapt on ${meta.id}..."

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("*_cutadapt_R1.fq.gz"), path("*_cutadapt_R2.fq.gz"), emit: trimmed_reads
    tuple val(meta), path("*_cutadapt.log"), emit: cutadapt_log

    script:
    """
    cutadapt -j ${task.cpus} -m ${params.min_read_length}:${params.min_read_length} \
        -a '${params.adapter_r1}' -A '${params.adapter_r2}' \
        --pair-filter=any --overlap ${params.adapter_overlap} \
        -o ${meta.id}_cutadapt_R1.fq.gz -p ${meta.id}_cutadapt_R2.fq.gz \
        ${read1} ${read2} > "${meta.id}_cutadapt.log"
    """
}

process HICPRO {
    tag "Running HiC-Pro on ${meta.id}..."

    cpus params.hicpro_cpus

    input:
    tuple val(meta), path(r1), path(r2)
    path chrom_sizes
    path fragments

    output:
    tuple val(meta), path('hicpro'), emit: results

    script:
    def settings = [TMP_DIR: 'tmp', LOGS_DIR: 'logs', BOWTIE2_OUTPUT_DIR: 'bowtie_results',
        MAPC_OUTPUT: 'hic_results', RAW_DIR: 'rawdata', N_CPU: task.cpus,
        LOGFILE: 'hicpro.log', JOB_NAME: '', JOB_MEM: '', JOB_WALLTIME: '', JOB_QUEUE: '', JOB_MAIL: '',
        PAIR1_EXT: '_R1', PAIR2_EXT: '_R2',
        BOWTIE2_IDX_PATH: file(params.bowtie2_index).parent,
        REFERENCE_GENOME: file(params.bowtie2_index).name,
        GENOME_SIZE: file(params.chrom_sizes).toAbsolutePath(), GENOME_FRAGMENT: file(params.genome_fragment).toAbsolutePath()]
    // Every scientific setting comes from Nextflow params, never the legacy config file.
    ['sort_ram','min_mapq','bowtie2_global_options','bowtie2_local_options','ligation_site',
     'allele_specific_snp','capture_target','report_capture_reporter','min_frag_size','max_frag_size',
     'min_insert_size','max_insert_size','min_cis_dist','get_all_interaction_classes',
     'get_process_sam','rm_singleton','rm_multi','rm_dup','matrix_format','max_iter',
     'filter_low_count_perc','filter_high_count_perc','eps'].each { String key -> settings[key.toUpperCase()] = params[key] }
    settings.BIN_SIZE = params.bin_size.toString().replaceAll(/[,\s]+/, ' ').trim()
    def configText = settings.collect { String key, def value -> "${key} = ${value}" }.join('\n')

    """
    mkdir -p input/${meta.id}
    ln -s "\$(realpath '${r1}')" input/${meta.id}/${meta.id}_R1.fastq.gz
    ln -s "\$(realpath '${r2}')" input/${meta.id}/${meta.id}_R2.fastq.gz
    cat > config-hicpro.txt <<'HICPRO_CONFIG'
${configText}
HICPRO_CONFIG
    '${params.hicpro}' -i input -o hicpro -c config-hicpro.txt
    cp config-hicpro.txt hicpro/config-hicpro.txt
    """
}

process FITHICHIP_PEAKS {
    tag "Calling peaks and generating matched 1D signal for ${meta.id}..."

    cpus 2

    input:
    tuple val(meta), path(hicpro), path(provided_peaks)
    path chrom_sizes

    output:
    tuple val(meta), path('peaks/peaks.bed'), emit: peaks
    tuple val(meta), path('peaks'), emit: results
    tuple val(meta), path('bw/*.bw'), emit: bw

    script:
    if (params.fithichip_read_length.toInteger() <= 0)
        error 'fithichip_read_length must be positive to generate PLAC-seq peaks and 1D signal'
    def peak_type = resolvePeakType(meta.id, params.fithichip_peak_types)
    def macs2_options = params.fithichip_macs2_options.toString()
    if (macs2_options =~ /(^|\s)--broad(?:-cutoff)?(?:=|\s|$)/)
        error 'Set fithichip_peak_types and fithichip_macs2_broad_cutoff instead of adding --broad options to fithichip_macs2_options'
    if (peak_type == 'broad' && (macs2_options =~ /(^|\s)--call-summits(?:\s|$)/))
        error '--call-summits is incompatible with broad peak calling; remove it from the shared fithichip_macs2_options'
    def broad_cutoff = params.fithichip_macs2_broad_cutoff.toDouble()
    if (!(broad_cutoff > 0 && broad_cutoff < 1))
        error 'fithichip_macs2_broad_cutoff must be between 0 and 1'
    def type_options = peak_type == 'broad' ? "--broad --broad-cutoff ${broad_cutoff}" : ''
    def peak_extension = peak_type == 'broad' ? 'broadPeak' : 'narrowPeak'
    def converter = params.bedgraph_to_bigwig ? file(params.bedgraph_to_bigwig).toAbsolutePath().toString() : ''
    def peak_command = provided_peaks ? """
    if [[ '${provided_peaks}' == *.gz ]]; then
        gzip -cd '${provided_peaks}' > peaks/peaks.bed
    else
        cp '${provided_peaks}' peaks/peaks.bed
    fi
    echo 'Loop anchors supplied externally; bigWig represents the inferred MACS2 signal.' >> peaks/peak_calling.log
    """ : """
    cp peaks/MACS2_ExtSize/out_macs2_peaks.${peak_extension} peaks/peaks.bed
    """
    """
    set -euo pipefail
    converter='${converter}'
    if [[ -z "\$converter" ]]; then
        if [[ -n "\${CONDA_PREFIX:-}" && -x "\${CONDA_PREFIX}/bin/bedGraphToBigWig" ]]; then
            converter="\${CONDA_PREFIX}/bin/bedGraphToBigWig"
        else
            converter=\$(command -v bedGraphToBigWig || true)
        fi
    fi
    [[ -n "\$converter" && -x "\$converter" ]] || {
        echo 'Set bedgraph_to_bigwig to an executable path in nextflow.config, or install bedGraphToBigWig in the active conda environment.' >&2
        exit 1
    }

    pair_dir='${hicpro}/hic_results/data/${meta.id}'
    valid_pairs="\$pair_dir/${meta.id}.allValidPairs"
    test -f "\$valid_pairs" || { echo "Missing duplicate-filtered HiC-Pro allValidPairs: \$valid_pairs" >&2; exit 1; }
    shopt -s nullglob
    class_files=( "\$pair_dir"/*.DEPairs "\$pair_dir"/*.SCPairs "\$pair_dir"/*.REPairs )
    mkdir -p peaks/MACS2_ExtSize bw
    printf 'metric\tvalue\nvalid_pairs_source\t%s\n' "\$valid_pairs" > peaks/read_selection.tsv
    printf 'macs2_peak_type\t%s\n' '${peak_type}' >> peaks/read_selection.tsv
    printf 'macs2_broad_cutoff\t%s\n' '${peak_type == 'broad' ? broad_cutoff : 'not_applicable'}' >> peaks/read_selection.tsv

    # FitHiChIP study: DE, SC, RE and duplicate-removed cis valid pairs <1 kb.
    # Select the named allValidPairs directly, avoiding PeakInferHiChIP.sh's
    # fallback to pre-deduplication *.validPairs for non-rawdata sample names.
    # Also deduplicate DE/SC/RE molecules using HiC-Pro's two-position key.
    {
        for pairs in "\${class_files[@]}"; do
            awk 'BEGIN {OFS=sprintf("%c",9)} NF>=7 {print \$2,\$3,\$4,\$5,\$6,\$7}' "\$pairs"
        done
        awk 'BEGIN {OFS=sprintf("%c",9)} NF>=7 && \$2==\$5 {
            distance=\$6-\$3
            if (distance<0) distance=-distance
            if (distance<1000) print \$2,\$3,\$4,\$5,\$6,\$7
        }' "\$valid_pairs"
    } | awk 'BEGIN {OFS=sprintf("%c",9)} NF>=6 {
        if (\$1>\$4 || (\$1==\$4 && \$2+0>\$5+0))
            print \$4,\$5,\$6,\$1,\$2,\$3
        else
            print \$1,\$2,\$3,\$4,\$5,\$6
        n++
    } END {
        print "eligible_pairs_before_dedup",n+0 >> "peaks/read_selection.tsv"
    }' | LC_ALL=C sort -T . -k1,1 -k2,2n -k4,4 -k5,5n -u > selected_pairs.tsv
    test -s selected_pairs.tsv || { echo 'No PLAC-seq pairs available for inferred peaks and 1D signal' >&2; exit 1; }
    printf 'selected_pairs_after_dedup\t%s\n' "\$(wc -l < selected_pairs.tsv)" >> peaks/read_selection.tsv

    # Pair positions are one-based read 5' coordinates; BED uses zero-based starts.
    # Retain strand for MACS2 fragment extension and clip reads to the reference.
    awk -v read_len=${params.fithichip_read_length} -v sizes='${chrom_sizes}' '
        BEGIN {
            OFS=sprintf("%c",9)
            while ((getline < sizes)>0) chrom_len[\$1]=\$2
            close(sizes)
        }
        function emit(chr,pos,strand, start,end) {
            if (!(chr in chrom_len) || pos !~ /^[0-9]+\$/) return
            if (strand=="+") {start=pos-1; end=start+read_len}
            else if (strand=="-") {end=pos; start=end-read_len}
            else return
            if (start<0) start=0
            if (end>chrom_len[chr]) end=chrom_len[chr]
            if (end>start) print chr,start,end,".",0,strand
        }
        {emit(\$1,\$2,\$3); emit(\$4,\$5,\$6)}
    ' selected_pairs.tsv | LC_ALL=C sort -T . -k1,1 -k2,2n > peaks/MACS2_input.bed
    test -s peaks/MACS2_input.bed || { echo 'No selected reads overlap the supplied chromosome sizes' >&2; exit 1; }
    printf 'retained_read_ends\t%s\n' "\$(wc -l < peaks/MACS2_input.bed)" >> peaks/read_selection.tsv

    # Use configured MACS2 settings and add broad mode by target class.
    # -B --SPMR save the same treatment pileup used for peaks, per million tags.
    # MACS2's default single-end tag filtering is shared by peaks and bigWig.
    macs2 callpeak -t peaks/MACS2_input.bed -f BED -n out_macs2 \
        --outdir peaks/MACS2_ExtSize -g '${params.fithichip_peak_genome}' \
        ${macs2_options} ${type_options} -B --SPMR > peaks/peak_calling.log 2>&1
    ${peak_command}
    test -s peaks/peaks.bed || { echo 'No peaks available for FitHiChIP' >&2; exit 1; }

    # MACS2 can extend fragments beyond chromosome ends; clip before conversion.
    awk 'BEGIN {OFS=sprintf("%c",9)} NR==FNR {sizes[\$1]=\$2; next}
        \$1 in sizes {
            start=\$2<0 ? 0 : \$2
            end=\$3>sizes[\$1] ? sizes[\$1] : \$3
            if (end>start) print \$1,start,end,\$4
        }' '${chrom_sizes}' peaks/MACS2_ExtSize/out_macs2_treat_pileup.bdg \
        | LC_ALL=C sort -T . -k1,1 -k2,2n > peaks/plac_1d.cpm.bedGraph
    test -s peaks/plac_1d.cpm.bedGraph
    "\$converter" peaks/plac_1d.cpm.bedGraph '${chrom_sizes}' bw/${meta.id}.plac_1d.cpm.bw
    printf 'bedGraphToBigWig\t%s\n' "\$converter" >> peaks/read_selection.tsv
    """
}

process FITHICHIP {
    tag "Running FitHiChIP ${mode} on ${meta.id}..."

    cpus params.fithichip_cpus

    input:
    tuple val(meta), path(hicpro), path(peaks), val(background), val(mode)
    path chrom_sizes

    output:
    tuple val(meta), path("FitHiChIP_${mode}"), emit: results

    script:
    def settings = [
        CircularGenome: params.fithichip_circular_genome,
        IntType: params.fithichip_int_type,
        BINSIZE: params.fithichip_bin_size,
        LowDistThr: params.fithichip_low_dist,
        UppDistThr: params.fithichip_upp_dist,
        UseP2PBackgrnd: background,
        BiasType: params.fithichip_bias_type,
        MergeInt: params.fithichip_merge_int,
        QVALUE: params.fithichip_qvalue,
        PREFIX: params.fithichip_prefix,
        OverWrite: params.fithichip_overwrite,
        NBins: params.fithichip_nbins,
        Draw: params.fithichip_draw,
        TimeProf: params.fithichip_time_prof,
        RefFasta: params.fithichip_ref_fasta,
        MappabilityFile: params.fithichip_mappability_file,
        REFragFile: params.fithichip_re_fragment,
        GCSize: params.fithichip_gc_size,
        MappSize: params.fithichip_mapp_size
    ]
    def config_text = settings.collect { key, value -> "${key}=${value}" }.join('\n')
    def matrix_prefix = "${hicpro}/hic_results/matrix/${meta.id}/raw/${params.fithichip_bin_size}/${meta.id}_${params.fithichip_bin_size}"

    """
    export PATH='${file(params.hicpro).parent}':\$PATH
    export OMP_NUM_THREADS=${task.cpus}
    mkdir -p FitHiChIP_${mode}
    # FitHiChIP changes into its installation directory. Keep implicit R plots
    # writable and separate for every task/R process without editing FitHiChIP.
    export FITHICHIP_PREVIOUS_R_PROFILE="\${R_PROFILE_USER:-}"
    export FITHICHIP_RPLOT_DIR="\$PWD/FitHiChIP_${mode}/Rplots"
    mkdir -p "\$FITHICHIP_RPLOT_DIR"
    cat > fithichip-graphics.Rprofile <<'R_PROFILE'
previous <- Sys.getenv("FITHICHIP_PREVIOUS_R_PROFILE")
if (!nzchar(previous)) {
    previous <- if (file.exists(".Rprofile")) ".Rprofile" else path.expand("~/.Rprofile")
}
if (file.exists(previous)) sys.source(previous, envir = .GlobalEnv)
options(device = function(...) grDevices::pdf(
    file = file.path(Sys.getenv("FITHICHIP_RPLOT_DIR"), paste0("Rplots.", Sys.getpid(), ".pdf")), ...))
R_PROFILE
    export R_PROFILE_USER="\$PWD/fithichip-graphics.Rprofile"
    cat > FitHiChIP_${mode}/configuration_file <<FITHICHIP_CONFIG
ValidPairs=\$(realpath '${hicpro}/hic_results/data/${meta.id}/${meta.id}.allValidPairs')
Interval=\$(realpath '${matrix_prefix}_abs.bed')
Matrix=\$(realpath '${matrix_prefix}.matrix')
ChrSizeFile=\$(realpath '${chrom_sizes}')
PeakFile=\$(realpath '${peaks}')
OutDir=\$PWD/FitHiChIP_${mode}
${config_text}
FITHICHIP_CONFIG
    test -s '${matrix_prefix}_abs.bed'
    test -s '${matrix_prefix}.matrix'
    bash '${params.fithichip_dir}/FitHiChIP_HiCPro.sh' \
        -C "\$PWD/FitHiChIP_${mode}/configuration_file" > FitHiChIP_${mode}/run.log 2>&1
    # Zero significant loops is a valid outcome; verify scoring rather than a success slogan.
    if grep -Eq 'Execution halted|Traceback \\(most recent call last\\)|^Error in |^ERROR' FitHiChIP_${mode}/run.log; then
        echo 'FitHiChIP reported an internal error; inspect run.log' >&2
        exit 1
    fi
    grep -q 'FINISHED calling significant interactions' FitHiChIP_${mode}/run.log
    test -s FitHiChIP_${mode}/Summary_results_${params.fithichip_prefix}.html
    find FitHiChIP_${mode} -name '${params.fithichip_prefix}.interactions_FitHiC.bed' -size +0c | grep -q .
    """
}

workflow {
    main:
    if (!(params.bin_size.toString() ==~ /[0-9]+([,\s]+[0-9]+)*/))
        error 'bin_size must be positive integer resolutions separated by commas or spaces'
    if (params.bin_size.toString().split(/[,\s]+/).any { String bin -> bin.toLong() <= 0 })
        error 'bin_size resolutions must be positive'
    def index = params.bowtie2_index
    if (!['bt2', 'bt2l'].any { String ext -> ['1','2','3','4','rev.1','rev.2'].every { String part -> file("${index}.${part}.${ext}").exists() } })
        error "Incomplete Bowtie2 index: ${index}"
    ch_read_pairs = channel.fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Invalid sample name: ${row.sample}"
            if (!row.fastq_1 || !row.fastq_2) error 'CSV requires sample,fastq_1,fastq_2 columns'
            // Validate each sample's target assignment before its FASTQ processing.
            resolvePeakType(row.sample, params.fithichip_peak_types)
            [
                row.sample,
                row
            ]
        }
        .groupTuple()
        .map { String _sample, List rows ->
            rows.withIndex().collect { def row, int lane_index ->
                row + [rep: lane_index + 1]
            }
        }
        .flatMap { def item -> item }
        .map { def row ->
            [
                [
                    id: row.sample,
                    rep: row.rep,
                ],
                [
                    file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true)
                ]
            ]
        }
        .map { Map meta, List files -> [meta.subMap(['id']), files] }
        .groupTuple()
        .map { Map meta, List filePairs ->
            [
                meta,
                filePairs.collect { List pair -> pair[0] },
                filePairs.collect { List pair -> pair[1] }
            ]
        }

    MERGE_FQ(ch_read_pairs)
    CUTADAPT(MERGE_FQ.out.merged_fq)
    HICPRO(CUTADAPT.out.trimmed_reads, file(params.chrom_sizes, checkIfExists: true), file(params.genome_fragment, checkIfExists: true))

    if (!params.bin_size.toString().split(/[,\s]+/).contains(params.fithichip_bin_size.toString()))
        error 'fithichip_bin_size must match a HiC-Pro bin_size so Interval and Matrix can be reused'
    if (params.rm_dup.toInteger() != 1)
        error 'rm_dup must be 1: loop calling and short-range 1D peak inference require duplicate-removed allValidPairs'
    if (params.fithichip_read_length.toInteger() <= 0)
        error 'Set a positive fithichip_read_length for PLAC-seq peak inference and 1D signal'

    ch_peak_inputs = HICPRO.out.results.map { meta, hicpro ->
        def peak_file = params.fithichip_peak_file instanceof Map ? params.fithichip_peak_file[meta.id] : params.fithichip_peak_file
        [meta, hicpro, peak_file ? file(peak_file, checkIfExists: true) : []]
    }
    FITHICHIP_PEAKS(ch_peak_inputs, file(params.chrom_sizes, checkIfExists: true))

    ch_fithichip_inputs = HICPRO.out.results.join(FITHICHIP_PEAKS.out.peaks)
        .flatMap { meta, hicpro, peaks ->
            [[meta, hicpro, peaks, 0, 'L'], [meta, hicpro, peaks, 1, 'S']]
        }
    FITHICHIP(ch_fithichip_inputs, file(params.chrom_sizes, checkIfExists: true))

    publish:
    out_merged_fastqs = MERGE_FQ.out.merged_fq
    out_cutadapt_logs = CUTADAPT.out.cutadapt_log
    out_hicpro = HICPRO.out.results
    out_bw = FITHICHIP_PEAKS.out.bw
    out_fithichip_peaks = FITHICHIP_PEAKS.out.results
    out_fithichip = FITHICHIP.out.results
}

output {
    out_fithichip_peaks {
        path { meta, _peaks -> "${meta.id}/FitHiChIP_peaks" }
    }
    out_fithichip {
        path { meta, _results -> "${meta.id}" }
    }
    out_merged_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_cutadapt_logs {
        path { meta, _f1 -> "${meta.id}/fastqs" }
    }
    out_hicpro {
        path { meta, _results -> "${meta.id}" }
    }
    out_bw {
        path { meta, _bw -> "${meta.id}" }
    }
}
