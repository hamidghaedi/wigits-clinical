nextflow.enable.dsl = 2

params.manifest = null
params.outdir = null
params.picard_image = null
params.fastp_image = null
params.multiqc_image = null
params.release_tag = 'lab-qc/v0.1.1'
params.release_commit = null

process LAB_FASTQ_QC {
    tag "${sample}:${lane}"
    container params.fastp_image
    cpus 4
    memory '6 GB'
    time { Math.max(1, Math.ceil(input_gib / 12.0)) * 2.h }
    maxForks 12
    publishDir "${params.outdir}/${run}/qc/${sample}/fastp", mode: 'copy'
    input:
    tuple val(run), val(sample), val(lane), val(input_gib), path(r1), path(r2)
    output:
    tuple val(run), path("${sample}_${lane}.fastp.json"), path("${sample}_${lane}.fastp.html"), path("${sample}_${lane}.fastp.version.txt"), emit: files
    script:
    """
    fastp --version > ${sample}_${lane}.fastp.version.txt 2>&1
    fastp --in1 ${r1} --in2 ${r2} --thread ${task.cpus} \
      --disable_adapter_trimming --disable_quality_filtering \
      --disable_length_filtering --disable_trim_poly_g \
      --json ${sample}_${lane}.fastp.json --html ${sample}_${lane}.fastp.html
    """
}

process LAB_WGS_METRICS {
    tag "$sample"
    container params.picard_image
    cpus 4
    memory '16 GB'
    time { Math.max(2, Math.ceil(input_gib / 20.0)) * 2.h }
    maxForks 4
    publishDir "${params.outdir}/${run}/qc/${sample}/picard", mode: 'copy'
    input:
    tuple val(run), val(sample), val(input_gib), path(bam), path(bai), path(reference), path(fai), path(dict), path(intervals)
    output:
    tuple val(run), path("${sample}.wgs_metrics.txt"), path("${sample}.picard.version.txt"), emit: files
    script:
    """
    java -jar /usr/picard/picard.jar CollectWgsMetrics --version > ${sample}.picard.version.txt 2>&1 || test \$? -eq 1
    java -Xmx12g -jar /usr/picard/picard.jar CollectWgsMetrics \
      I=${bam} O=${sample}.wgs_metrics.txt R=${reference} INTERVALS=${intervals} \
      MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=20 COVERAGE_CAP=250 \
      COUNT_UNPAIRED=false USE_FAST_ALGORITHM=true READ_LENGTH=151
    """
}

process LAB_INSERT_METRICS {
    tag "$sample"
    container params.picard_image
    cpus 2
    memory '8 GB'
    time { Math.max(1, Math.ceil(input_gib / 40.0)) * 2.h }
    maxForks 4
    publishDir "${params.outdir}/${run}/qc/${sample}/picard", mode: 'copy'
    input:
    tuple val(run), val(sample), val(input_gib), path(bam), path(bai)
    output:
    tuple val(run), path("${sample}.insert_size_metrics.txt"), path("${sample}.insert_size_histogram.pdf"), emit: files
    script:
    """
    java -Xmx6g -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
      I=${bam} O=${sample}.insert_size_metrics.txt H=${sample}.insert_size_histogram.pdf
    """
}

process LAB_RUN_REPORT {
    tag "$run"
    container params.multiqc_image
    cpus 1
    memory '6 GB'
    time '1h'
    publishDir "${params.outdir}/${run}", mode: 'copy'
    input:
    tuple val(run), path(qc_files, stageAs: 'metrics/*')
    path manifest
    path logo
    path css
    path helper
    output:
    path '*.multiqc_report.html'
    path '*_data.zip'
    path 'provenance'
    path '_SUCCESS'
    script:
    """
    multiqc --version > multiqc.version.txt
    python3 ${helper} prepare --manifest ${manifest} --run '${run}' \
      --metrics metrics --logo ${logo} --css ${css} \
      --release-tag '${params.release_tag}' --release-commit '${params.release_commit}'
    multiqc qc --config multiqc_config.json --force --zip-data-dir \
      --filename '${run}.multiqc_report.html' --outdir .
    python3 ${helper} finalize --manifest ${manifest} --run '${run}' \
      --release-tag '${params.release_tag}' --release-commit '${params.release_commit}'
    """
}

workflow LAB_RUN_QC {
    if (!params.manifest || !params.outdir || !params.release_commit) error 'Supply --manifest, --outdir and --release_commit'
    def manifest = new groovy.json.JsonSlurper().parse(file(params.manifest).toFile())
    def samples = manifest.samples
    if (!samples || samples.collect { it.sample_id }.unique().size() != samples.size()) error 'Expected unique laboratory samples'
    ch_samples = Channel.fromList(samples)
    ch_fastqs = ch_samples.flatMap { s -> s.fastqs.collect { q ->
        tuple(s.run_id, s.sample_id, q.label, q.input_gib, file(q.r1, checkIfExists:true), file(q.r2, checkIfExists:true))
    } }
    ch_wgs = ch_samples.map { s -> tuple(s.run_id, s.sample_id, s.bam_gib,
        file(s.bam, checkIfExists:true), file(s.bai, checkIfExists:true),
        file(manifest.reference, checkIfExists:true), file(manifest.fai, checkIfExists:true),
        file(manifest.dict, checkIfExists:true), file(manifest.intervals, checkIfExists:true)) }
    ch_insert = ch_samples.map { s -> tuple(s.run_id, s.sample_id, s.bam_gib,
        file(s.bam, checkIfExists:true), file(s.bai, checkIfExists:true)) }
    LAB_FASTQ_QC(ch_fastqs)
    LAB_WGS_METRICS(ch_wgs)
    LAB_INSERT_METRICS(ch_insert)
    ch_metrics = LAB_FASTQ_QC.out.files.mix(LAB_WGS_METRICS.out.files).mix(LAB_INSERT_METRICS.out.files)
        .map { entry -> tuple(entry[0], entry[1..-1]) }.groupTuple()
        .map { run, groups -> tuple(run, groups.flatten()) }
    LAB_RUN_REPORT(ch_metrics, file(params.manifest), file(manifest.logo), file(manifest.css), file("${projectDir}/lab_qc_report.py"))
}
