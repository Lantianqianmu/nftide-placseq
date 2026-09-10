#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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

    publish:
    out_merged_fastqs = MERGE_FQ.out.merged_fq
    out_cutadapt_logs = CUTADAPT.out.cutadapt_log
    out_hicpro = HICPRO.out.results
}

output {
    out_merged_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_cutadapt_logs {
        path { meta, _f1 -> "${meta.id}/fastqs" }
    }
    out_hicpro {
        path { meta, _results -> "${meta.id}" }
    }
}
