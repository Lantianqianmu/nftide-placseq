#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def shellQuote(def value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def resolvePeakType(String sample, def peakTypes) {
    if (!(peakTypes instanceof Map))
        throw new IllegalArgumentException('fithichip_peak_types must map sample names or mark suffixes to broad or narrow')
    def types = peakTypes.collectEntries { key, value ->
        [(key.toString().trim().toLowerCase()): value.toString().trim().toLowerCase()]
    }
    def sample_key = sample.toLowerCase()
    def mark = sample_key.replaceFirst(/[_.-][0-9]+$/, '').tokenize('_.-').last().replaceFirst(/^h3k/, '')
    def peak_type = types[sample_key] ?: types[mark]
    if (!(peak_type in ['broad', 'narrow']))
        throw new IllegalArgumentException("Cannot determine MACS2 peak type for ${sample}; add its sample name or mark suffix to fithichip_peak_types as broad or narrow")
    return peak_type
}

def macs2Settings(String peak_type, def options, def cutoff) {
    def opts = options.toString()
    if (opts =~ /(^|\s)--broad(?:-cutoff)?(?:=|\s|$)/)
        throw new IllegalArgumentException('Set fithichip_peak_types and fithichip_macs2_broad_cutoff instead of adding --broad options to fithichip_macs2_options')
    if (peak_type == 'broad' && (opts =~ /(^|\s)--call-summits(?:\s|$)/))
        throw new IllegalArgumentException('--call-summits is incompatible with broad peak calling')
    def broad_cutoff = cutoff.toDouble()
    if (!(broad_cutoff > 0 && broad_cutoff < 1))
        throw new IllegalArgumentException('fithichip_macs2_broad_cutoff must be between 0 and 1')
    return [options: opts, type_options: peak_type == 'broad' ? "--broad --broad-cutoff ${broad_cutoff}" : '',
            extension: peak_type == 'broad' ? 'broadPeak' : 'narrowPeak', broad_cutoff: broad_cutoff]
}

def resolveLocalPath(def value, java.nio.file.Path base) {
    def path = java.nio.file.Paths.get(value.toString().trim())
    return (path.isAbsolute() ? path : base.resolve(path)).toAbsolutePath().normalize()
}

def prepareSamples(List rows, java.nio.file.Path csv, def peakTypes) {
    if (!rows) throw new IllegalArgumentException('The samplesheet contains no samples')
    def columns = ['sample', 'fastq_1', 'fastq_2', 'peak_file']
    rows.each { row ->
        if (!row.keySet().containsAll(columns))
            throw new IllegalArgumentException('CSV requires sample,fastq_1,fastq_2,peak_file columns; use an empty peak_file field for per-sample inference')
        row.sample = row.sample?.toString()?.trim()
        if (!row.sample || !(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
            throw new IllegalArgumentException("Invalid sample name: ${row.sample}")
        ['fastq_1','fastq_2'].each { key ->
            if (!row[key]?.toString()?.trim()) throw new IllegalArgumentException("Missing ${key} for ${row.sample}")
            row[key] = resolveLocalPath(row[key], csv.parent)
            if (!java.nio.file.Files.isRegularFile(row[key]) || !java.nio.file.Files.isReadable(row[key]))
                throw new IllegalArgumentException("Missing or unreadable FASTQ: ${row[key]}")
        }
        def requested = row.peak_file?.toString()?.trim()
        row.peak_file = requested ? resolveLocalPath(requested, csv.parent).toString() : ''
    }
    def samples = rows.groupBy { it.sample }.collect { String id, List lanes ->
        def requested = lanes.collect { it.peak_file }.unique()
        if (requested.size() != 1)
            throw new IllegalArgumentException("All rows for sample ${id} must have the same peak_file value, including blank fields")
        def peak_path = requested[0]
        def mode = 'sample'
        if (peak_path) {
            def p = java.nio.file.Paths.get(peak_path)
            if (java.nio.file.Files.exists(p, java.nio.file.LinkOption.NOFOLLOW_LINKS)) {
                if (!java.nio.file.Files.isRegularFile(p) || !java.nio.file.Files.isReadable(p) || java.nio.file.Files.size(p) == 0)
                    throw new IllegalArgumentException("Peak path exists but is not a readable, nonempty file: ${peak_path}")
                mode = 'provided'
            } else mode = 'shared'
        }
        def plan = [path: peak_path, mode: mode, peak_type: resolvePeakType(id, peakTypes)]
        [id, lanes.collect { it.fastq_1 }, lanes.collect { it.fastq_2 }, plan]
    }
    samples.findAll { it[3].mode == 'shared' }.groupBy { it[3].path }.each { String path, List members ->
        if (members.collect { it[3].peak_type }.unique().size() != 1)
            throw new IllegalArgumentException("Shared peak group ${path} mixes broad and narrow targets; assign separate peak_file paths")
    }
    return samples
}

process MERGE_FQ {
    tag "Merging fastq files of ${meta.id}..."

    input:
    tuple val(meta), path(r1s, stageAs: 'r1??/*'), path(r2s, stageAs: 'r2??/*')

    output:
    tuple val(meta), path("*_merged_R1.fq.gz"), path("*_merged_R2.fq.gz"), emit: merged_fq

    script:
    def forward = r1s instanceof Collection ? r1s : [r1s]
    def reverse = r2s instanceof Collection ? r2s : [r2s]
    """
    cat ${forward.collect { shellQuote(it) }.join(' ')} > ${meta.id}_merged_R1.fq.gz
    cat ${reverse.collect { shellQuote(it) }.join(' ')} > ${meta.id}_merged_R2.fq.gz
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
        ${shellQuote(read1)} ${shellQuote(read2)} > "${meta.id}_cutadapt.log"
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

process PLAC_SIGNAL {
    tag "Generating sample peak input and 1D signal for ${meta.id}..."

    cpus 2

    input:
    tuple val(meta), path(hicpro)
    path chrom_sizes

    output:
    tuple val(meta), path('signal/inferred_peaks.bed'), emit: inferred_peaks
    tuple val(meta), path('signal/MACS2_input.bed'), emit: peak_reads
    tuple val(meta), path('signal'), emit: results
    tuple val(meta), path('bw/*.bw'), emit: bw

    script:
    if (params.fithichip_read_length.toInteger() <= 0)
        error 'fithichip_read_length must be positive to generate PLAC-seq peaks and 1D signal'
    def peak_type = meta.peak_type
    def settings = macs2Settings(peak_type, params.fithichip_macs2_options, params.fithichip_macs2_broad_cutoff)
    def macs2_options = settings.options
    def type_options = settings.type_options
    def peak_extension = settings.extension
    def broad_cutoff = settings.broad_cutoff
    def converter = params.bedgraph_to_bigwig ? file(params.bedgraph_to_bigwig).toAbsolutePath().toString() : ''
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
    mkdir -p signal/MACS2_ExtSize bw
    printf 'metric\tvalue\nvalid_pairs_source\t%s\n' "\$valid_pairs" > signal/read_selection.tsv
    printf 'macs2_peak_type\t%s\n' '${peak_type}' >> signal/read_selection.tsv
    printf 'macs2_broad_cutoff\t%s\n' '${peak_type == 'broad' ? broad_cutoff : 'not_applicable'}' >> signal/read_selection.tsv

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
        print "eligible_pairs_before_dedup",n+0 >> "signal/read_selection.tsv"
    }' | LC_ALL=C sort -T . -k1,1 -k2,2n -k4,4 -k5,5n -u > selected_pairs.tsv
    test -s selected_pairs.tsv || { echo 'No PLAC-seq pairs available for inferred peaks and 1D signal' >&2; exit 1; }
    printf 'selected_pairs_after_dedup\t%s\n' "\$(wc -l < selected_pairs.tsv)" >> signal/read_selection.tsv

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
    ' selected_pairs.tsv | LC_ALL=C sort -T . -k1,1 -k2,2n > signal/MACS2_input.bed
    test -s signal/MACS2_input.bed || { echo 'No selected reads overlap the supplied chromosome sizes' >&2; exit 1; }
    printf 'retained_read_ends\t%s\n' "\$(wc -l < signal/MACS2_input.bed)" >> signal/read_selection.tsv

    # Use configured MACS2 settings and add broad mode by target class.
    # -B --SPMR save the same treatment pileup used for peaks, per million tags.
    # MACS2's default single-end tag filtering is shared by peaks and bigWig.
    macs2 callpeak -t signal/MACS2_input.bed -f BED -n out_macs2 \
        --outdir signal/MACS2_ExtSize -g '${params.fithichip_peak_genome}' \
        ${macs2_options} ${type_options} -B --SPMR > signal/peak_calling.log 2>&1
    cp signal/MACS2_ExtSize/out_macs2_peaks.${peak_extension} signal/inferred_peaks.bed

    # MACS2 can extend fragments beyond chromosome ends; clip before conversion.
    awk 'BEGIN {OFS=sprintf("%c",9)} NR==FNR {sizes[\$1]=\$2; next}
        \$1 in sizes {
            start=\$2<0 ? 0 : \$2
            end=\$3>sizes[\$1] ? sizes[\$1] : \$3
            if (end>start) print \$1,start,end,\$4
        }' '${chrom_sizes}' signal/MACS2_ExtSize/out_macs2_treat_pileup.bdg \
        | LC_ALL=C sort -T . -k1,1 -k2,2n > signal/plac_1d.cpm.bedGraph
    test -s signal/plac_1d.cpm.bedGraph
    "\$converter" signal/plac_1d.cpm.bedGraph '${chrom_sizes}' bw/${meta.id}.plac_1d.cpm.bw
    printf 'bedGraphToBigWig\t%s\n' "\$converter" >> signal/read_selection.tsv
    """
}


process SHARED_PEAKS {
    tag "Calling shared peaks for ${group_info.id}..."
    cpus 2

    publishDir { java.nio.file.Paths.get(group_info.path).parent.toString() },
        mode: 'copy', pattern: 'export/requested_peak_file',
        saveAs: { ignored -> java.nio.file.Paths.get(group_info.path).fileName.toString() },
        overwrite: false, failOnError: true

    input:
    tuple val(group_info), path(read_beds, stageAs: 'member??/*')
    path chrom_sizes

    output:
    tuple val(group_info), path('peaks/peaks.bed'), emit: peaks
    tuple val(group_info), path('peaks'), emit: results
    path('export/requested_peak_file'), emit: exported

    script:
    def settings = macs2Settings(group_info.peak_type, params.fithichip_macs2_options, params.fithichip_macs2_broad_cutoff)
    def beds = read_beds instanceof Collection ? read_beds : [read_beds]
    def manifest = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(group_info))
    def export_command = group_info.path.toLowerCase().endsWith('.gz') ?
        'gzip -c peaks/peaks.bed > export/requested_peak_file' : 'cp peaks/peaks.bed export/requested_peak_file'
    """
    set -euo pipefail
    mkdir -p peaks/MACS2_ExtSize export
    cat ${beds.collect { shellQuote(it) }.join(' ')} \
        | LC_ALL=C sort -T . -k1,1 -k2,2n > peaks/MACS2_input.bed
    test -s peaks/MACS2_input.bed
    macs2 callpeak -t peaks/MACS2_input.bed -f BED -n out_macs2 \
        --outdir peaks/MACS2_ExtSize -g '${params.fithichip_peak_genome}' \
        ${settings.options} ${settings.type_options} -B --SPMR > peaks/peak_calling.log 2>&1
    prepare_peaks.py --input peaks/MACS2_ExtSize/out_macs2_peaks.${settings.extension} \
        --chrom-sizes ${shellQuote(chrom_sizes)} --output peaks/peaks.bed \
        --report peaks/validation.json
    cat > peaks/group.json <<'PEAK_GROUP_JSON'
${manifest}
PEAK_GROUP_JSON
    ${export_command}
    """
}

process RESOLVE_PEAKS {
    tag "Assigning ${meta.source} peaks to ${meta.id}..."
    input:
    tuple val(meta), path(selected_peaks, stageAs: 'input_peaks')
    path chrom_sizes

    output:
    tuple val(meta), path('peaks/peaks.bed'), emit: peaks
    tuple val(meta), path('peaks'), emit: results

    script:
    def manifest = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(meta))
    """
    set -euo pipefail
    mkdir -p peaks
    prepare_peaks.py --input ${shellQuote(selected_peaks)} \
        --chrom-sizes ${shellQuote(chrom_sizes)} --output peaks/peaks.bed \
        --report peaks/validation.json
    cat > peaks/source.json <<'PEAK_SOURCE_JSON'
${manifest}
PEAK_SOURCE_JSON
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
    if (!params.bin_size.toString().split(/[,\s]+/).contains(params.fithichip_bin_size.toString()))
        error 'fithichip_bin_size must match a HiC-Pro bin_size so Interval and Matrix can be reused'
    if (params.rm_dup.toInteger() != 1)
        error 'rm_dup must be 1: loop calling and short-range 1D peak inference require duplicate-removed allValidPairs'
    if (params.get_all_interaction_classes.toInteger() != 1)
        error 'get_all_interaction_classes must be 1 to retain DE/SC/RE reads for peaks and 1D signal'
    if (params.fithichip_read_length.toInteger() <= 0)
        error 'Set a positive fithichip_read_length for PLAC-seq peak inference and 1D signal'
    if (params.containsKey('fithichip_peak_file') && params.fithichip_peak_file)
        error 'Move fithichip_peak_file assignments to the samplesheet peak_file column'
    def index = params.bowtie2_index
    if (!['bt2', 'bt2l'].any { String ext -> ['1','2','3','4','rev.1','rev.2'].every { String part -> file("${index}.${part}.${ext}").exists() } })
        error "Incomplete Bowtie2 index: ${index}"
    def csv = file(params.input_csv, checkIfExists: true).toAbsolutePath()
    ch_samples = channel.fromPath(csv, checkIfExists: true)
        .splitCsv(header: true)
        .ifEmpty { error 'The samplesheet contains no samples' }
        .collect(flat: false)
        .flatMap { rows ->
            def samples = prepareSamples(rows, csv, params.fithichip_peak_types)
            samples.each { id, r1s, r2s, plan ->
                macs2Settings(plan.peak_type, params.fithichip_macs2_options, params.fithichip_macs2_broad_cutoff)
            }
            samples
        }

    // Alignment cache keys depend on sample IDs and reads, not peak-file status.
    ch_read_pairs = ch_samples.map { id, r1s, r2s, plan -> [[id: id], r1s, r2s] }
    ch_peak_plans = ch_samples.map { id, r1s, r2s, plan -> [id, plan] }
    MERGE_FQ(ch_read_pairs)
    CUTADAPT(MERGE_FQ.out.merged_fq)
    HICPRO(CUTADAPT.out.trimmed_reads, file(params.chrom_sizes, checkIfExists: true), file(params.genome_fragment, checkIfExists: true))

    ch_signal_inputs = HICPRO.out.results.map { meta, hicpro -> [meta.id, hicpro] }
        .join(ch_peak_plans)
        .map { id, hicpro, plan -> [[id: id, peak_type: plan.peak_type], hicpro] }
    PLAC_SIGNAL(ch_signal_inputs, file(params.chrom_sizes, checkIfExists: true))

    ch_individual_peaks = PLAC_SIGNAL.out.inferred_peaks.map { meta, peaks -> [meta.id, peaks] }
        .join(ch_peak_plans)
        .filter { id, inferred, plan -> plan.mode != 'shared' }
        .map { id, inferred, plan ->
            def meta = [id: id, source: plan.mode, requested_path: plan.path,
                        source_samples: plan.mode == 'sample' ? [id] : []]
            [meta, plan.mode == 'provided' ? java.nio.file.Paths.get(plan.path) : inferred]
        }

    ch_shared_inputs = PLAC_SIGNAL.out.peak_reads.map { meta, bed -> [meta.id, bed] }
        .join(ch_peak_plans)
        .filter { id, bed, plan -> plan.mode == 'shared' }
        .map { id, bed, plan -> [plan.path, plan.peak_type, id, bed] }
        .groupTuple()
        .map { requested, types, ids, beds ->
            def digest = java.security.MessageDigest.getInstance('SHA-256').digest(requested.getBytes('UTF-8')).encodeHex().toString().take(16)
            def entries = [ids, beds].transpose().sort { it[0] }
            [[id: "group_${digest}", path: requested, peak_type: types[0],
              source_samples: entries.collect { it[0] }], entries.collect { it[1] }]
        }
    SHARED_PEAKS(ch_shared_inputs, file(params.chrom_sizes, checkIfExists: true))
    ch_shared_assignments = SHARED_PEAKS.out.peaks.flatMap { group_info, peaks ->
        group_info.source_samples.collect { id ->
            [[id: id, source: 'shared', requested_path: group_info.path,
              source_samples: group_info.source_samples], peaks]
        }
    }
    RESOLVE_PEAKS(ch_individual_peaks.mix(ch_shared_assignments), file(params.chrom_sizes, checkIfExists: true))

    ch_fithichip_inputs = HICPRO.out.results.map { meta, hicpro -> [meta.id, meta, hicpro] }
        .join(RESOLVE_PEAKS.out.peaks.map { meta, peaks -> [meta.id, peaks] })
        .flatMap { id, meta, hicpro, peaks ->
            [[meta, hicpro, peaks, 0, 'L'], [meta, hicpro, peaks, 1, 'S']]
        }
    FITHICHIP(ch_fithichip_inputs, file(params.chrom_sizes, checkIfExists: true))

    publish:
    out_merged_fastqs = MERGE_FQ.out.merged_fq
    out_cutadapt_logs = CUTADAPT.out.cutadapt_log
    out_hicpro = HICPRO.out.results
    out_bw = PLAC_SIGNAL.out.bw
    out_signal = PLAC_SIGNAL.out.results
    out_fithichip_peaks = RESOLVE_PEAKS.out.results
    out_shared_peaks = SHARED_PEAKS.out.results
    out_fithichip = FITHICHIP.out.results
}

output {
    out_fithichip_peaks {
        mode 'copy'
        path { meta, _peaks -> "${meta.id}/FitHiChIP_peaks" }
    }
    out_signal {
        mode 'copy'
        path { meta, _signal -> "${meta.id}/PLAC_signal" }
    }
    out_shared_peaks {
        mode 'copy'
        path { group_info, _peaks -> "shared_peaks/${group_info.id}" }
    }
    out_fithichip {
        path { meta, _results -> "${meta.id}" }
    }
    out_merged_fastqs {
        mode 'copy'
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_cutadapt_logs {
        path { meta, _log -> "${meta.id}/fastqs" }
    }
    out_hicpro {
        mode 'copy'
        path { meta, _results -> "${meta.id}" }
    }
    out_bw {
        path { meta, _bw -> "${meta.id}" }
    }
}
