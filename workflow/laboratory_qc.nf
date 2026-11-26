nextflow.enable.dsl = 2

params.manifest = null
params.outdir = null
params.picard_image = null
params.parabricks_image = null
params.slurm_account = null
params.fastp_image = null
params.multiqc_image = null
params.release_tag = 'lab-qc/v0.4.1'
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
    container params.parabricks_image
    cpus 32
    memory { allowance.memory_gib.GB }
    time { allowance.minutes.min }
    // Preserve other running metrics if one fails; still fail the workflow.
    errorStrategy 'finish'
    // No artificial four-library cap: the site queue/scheduler controls capacity.
    publishDir "${params.outdir}/${run}/qc/${sample}/coverage", mode: 'copy'
    input:
    tuple val(run), val(sample), val(input_gib), val(allowance), path(bam), path(bai), path(reference), path(fai)
    output:
    tuple val(run), path("${sample}.wgs_metrics.txt"), path("${sample}.bammetrics.version.txt"), path("${sample}.coverage.bed"), emit: files
    script:
    """
    pbrun bammetrics --version > ${sample}.bammetrics.version.txt 2>&1
    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)\$/ {print \$1 "\\t0\\t" \$2}' ${fai} > ${sample}.coverage.bed
    pbrun bammetrics --ref ${reference} --bam ${bam} \
      --interval-file ${sample}.coverage.bed --out-metrics-file ${sample}.wgs_metrics.txt \
      --minimum-base-quality 20 --minimum-mapping-quality 20 --coverage-cap 250 \
      --num-threads ${task.cpus} --tmp-dir .
    """
}

def reuseScript(sources) {
    def quote = { value -> "'" + value.toString().replace("'", "'\\''") + "'" }
    'mkdir -p reused\n' + sources.collect { target, source ->
        "cp -- ${quote(source)} ${quote('reused/'+target)}"
    }.join('\n')
}

process LAB_REUSE_QC {
    executor 'local'
    cpus 1
    memory '512 MB'
    maxForks 2
    input:
    tuple val(run),val(sample),val(kind),val(sources),val(dependency)
    output:
    tuple val(run),path('reused/*'),emit:files
    script:
    reuseScript(sources)
}

process LAB_REUSE_QC_AFTER_JOB {
    cpus 1
    memory '1 GB'
    time '10min'
    clusterOptions { "--account=${params.slurm_account} --dependency=afterok:${dependency}" }
    input:
    tuple val(run),val(sample),val(kind),val(sources),val(dependency)
    output:
    tuple val(run),path('reused/*'),emit:files
    script:
    reuseScript(sources)
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
    tuple val(run), path("${sample}.insert_size_metrics.txt"), path("${sample}.insert_size_histogram.pdf"), path("${sample}.picard.version.txt"), emit: files
    script:
    """
    java -jar /usr/picard/picard.jar CollectInsertSizeMetrics --version > ${sample}.picard.version.txt 2>&1 || test \$? -eq 1
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
    path style
    output:
    path '*.multiqc_report.html'
    path '*_data.zip'
    path 'provenance'
    path '_SUCCESS', emit:success
    path 'qc', emit:source_metrics
    script:
    """
    mkdir -p assets
    cp ${style} assets/lab_multiqc_style.json
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
    if (!manifest.expected_samples) error 'Laboratory manifest expected_samples is required'
    if (samples.collect { it.run_id }.toSet() != manifest.expected_samples.keySet()) error 'Laboratory run set differs'
    samples.groupBy { it.run_id }.each { run, group ->
        if (group.collect { it.sample_id }.sort() != manifest.expected_samples[run]?.sort())
            error "Laboratory run membership differs: ${run}"
    }
    samples.each { s ->
        (s.reuse ?: [:]).each { kind, item ->
            if (!(kind in ['fastp','insert','coverage'])) error 'Unknown reused metric category'
            def allowed = kind == 'fastp' ? s.fastqs.collectMany { q ->
                ["${s.sample_id}_${q.label}.fastp.json", "${s.sample_id}_${q.label}.fastp.html", "${s.sample_id}_${q.label}.fastp.version.txt"]
            } : kind == 'insert' ? ["${s.sample_id}.insert_size_metrics.txt", "${s.sample_id}.insert_size_histogram.pdf", "${s.sample_id}.picard.version.txt"] :
                ["${s.sample_id}.wgs_metrics.txt", "${s.sample_id}.bammetrics.version.txt", "${s.sample_id}.coverage.interval_list", "${s.sample_id}.coverage.bed"]
            allowed=allowed.collect { it.toString() }
            if (!(item.files.keySet().every { it in allowed })) error "Unexpected reused metric: ${s.sample_id} ${kind}"
            if (kind == 'coverage' && item.producer != 'parabricks.bammetrics') error 'Reused coverage must identify Parabricks'
            if (item.dependency && !(item.dependency.toString() ==~ /[0-9]+/)) error 'Invalid upstream metric job binding'
        }
    }
    def sizes=samples.groupBy { it.run_id }.collectEntries { run, group ->
        [(run): group.sum { s -> (s.reuse?.fastp ? 1 : s.fastqs.size()) + 2 }]
    }
    ch_samples = Channel.fromList(samples)
    ch_fastqs = ch_samples.filter { !it.reuse?.fastp }.flatMap { s -> s.fastqs.collect { q ->
        tuple(s.run_id, s.sample_id, q.label, q.input_gib, file(q.r1, checkIfExists:true), file(q.r2, checkIfExists:true))
    } }
    ch_wgs = ch_samples.filter { !it.reuse?.coverage }.map { s -> tuple(s.run_id, s.sample_id, s.bam_gib,
        LabQcResources.coverage(s.bam_gib, s.bamtools_summary),
        file(s.bam, checkIfExists:true), file(s.bai, checkIfExists:true),
        file(manifest.reference, checkIfExists:true), file(manifest.fai, checkIfExists:true)) }
    ch_insert = ch_samples.filter { !it.reuse?.insert }.map { s -> tuple(s.run_id, s.sample_id, s.bam_gib,
        file(s.bam, checkIfExists:true), file(s.bai, checkIfExists:true)) }
    ch_reuse=ch_samples.flatMap { s -> (s.reuse ?: [:]).collect { kind,item ->
        tuple(s.run_id,s.sample_id,kind,item.files,item.dependency ?: '')
    } }.branch {
        immediate: !it[4]
        pending: true
    }
    LAB_FASTQ_QC(ch_fastqs)
    LAB_WGS_METRICS(ch_wgs)
    LAB_INSERT_METRICS(ch_insert)
    LAB_REUSE_QC(ch_reuse.immediate)
    LAB_REUSE_QC_AFTER_JOB(ch_reuse.pending)
    ch_metrics = LAB_FASTQ_QC.out.files.mix(LAB_WGS_METRICS.out.files).mix(LAB_INSERT_METRICS.out.files)
        .mix(LAB_REUSE_QC.out.files).mix(LAB_REUSE_QC_AFTER_JOB.out.files)
        .map { entry -> tuple(groupKey(entry[0],sizes[entry[0]]), entry[1..-1]) }.groupTuple()
        .map { key, groups -> tuple(key.groupTarget, groups.flatten()) }
    LAB_RUN_REPORT(ch_metrics, file(params.manifest), file(manifest.logo), file(manifest.css), file("${projectDir}/lab_qc_report.py"), file("${projectDir}/assets/lab_multiqc_style.json"))
    LAB_RUN_REPORT.out.success.subscribe { success -> LabQcDelivery.reportReady(success,log) }
}
