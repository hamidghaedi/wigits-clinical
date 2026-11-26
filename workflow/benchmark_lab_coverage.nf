nextflow.enable.dsl=2
params.manifest=null
params.outdir=null
params.parabricks_image=null

process ACCELERATED_COVERAGE {
    tag "$sample"
    container params.parabricks_image
    cpus 32
    memory '64 GB'
    time { Math.max(120,Math.ceil(bam_gib*2)).toInteger().min }
    maxForks 2
    publishDir "${params.outdir}/${sample}", mode:'copy'
    input:
    tuple val(sample), val(bam_gib), path(bam), path(bai), path(intervals), path(baseline)
    path reference
    path fai
    output:
    tuple val(sample),path("${sample}.parabricks.wgs_metrics.txt"),path(baseline),emit: metrics
    path '*.help.txt'
    path '*.version.txt'
    script:
    """
    pbrun bammetrics --help > ${sample}.help.txt 2>&1
    pbrun bammetrics --version > ${sample}.version.txt 2>&1
    pbrun bammetrics --ref ${reference} --bam ${bam} \
      --interval-file ${intervals} --out-metrics-file ${sample}.parabricks.wgs_metrics.txt \
      --minimum-base-quality 20 --minimum-mapping-quality 20 --coverage-cap 250 \
      --num-threads ${task.cpus} --tmp-dir .
    """
}

process COMPARE_COVERAGE {
    tag "$sample"
    cpus 1
    memory '1 GB'
    time '10min'
    publishDir "${params.outdir}/${sample}",mode:'copy'
    input:
    tuple val(sample),path(candidate),path(baseline)
    path helper
    output:
    path '*.comparison.json'
    script:
    """
    python3 ${helper} --baseline ${baseline} --candidate ${candidate} --out ${sample}.comparison.json
    """
}

workflow COVERAGE_BENCHMARK {
    def m=new groovy.json.JsonSlurper().parse(file(params.manifest).toFile())
    ch=Channel.fromList(m.samples).map { s -> tuple(s.sample_id,s.bam_gib,
       file(s.bam,checkIfExists:true),file(s.bai,checkIfExists:true),
       file(s.intervals,checkIfExists:true),file(s.baseline,checkIfExists:true)) }
    ACCELERATED_COVERAGE(ch,file(m.reference),file(m.fai))
    COMPARE_COVERAGE(ACCELERATED_COVERAGE.out.metrics,file("${projectDir}/compare_lab_coverage.py"))
}
