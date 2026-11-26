nextflow.enable.dsl=2

params.plan = null
params.outdir = null

process BLOOD_CH {
    tag "${meta.case_id}"
    cpus { meta.blood_ch.resources.cpus }
    memory { meta.blood_ch.resources.memory_gib.GB }
    time { meta.blood_ch.resources.minutes.min }
    errorStrategy 'finish'
    publishDir "${params.outdir}/${meta.case_id}", mode: 'copy'
    input:
    tuple val(meta), path(plan)
    output:
    tuple val(meta), path('blood_ch'), emit: results
    script:
    """
    python3 '${projectDir}/extension_tasks.py' --plan '${plan}' --case-id '${meta.case_id}' \
      --stage blood_ch --outdir blood_ch --threads ${task.cpus}
    """
    stub:
    """
    mkdir blood_ch
    echo '{"status":"STUB_COMPLETED","synthetic":true}' > blood_ch/status.json
    """
}

process NEO_PREPARE {
    tag "${meta.case_id}"
    cpus { meta.rna_neoantigens.resources.cpus }
    memory { meta.rna_neoantigens.resources.memory_gib.GB }
    time { meta.rna_neoantigens.resources.minutes.min }
    errorStrategy 'finish'
    publishDir "${params.outdir}/${meta.case_id}", mode: 'copy'
    input:
    tuple val(meta), path(plan)
    output:
    tuple val(meta), path('neo_export'), emit: prepared
    script:
    """
    python3 '${projectDir}/extension_tasks.py' --plan '${plan}' --case-id '${meta.case_id}' \
      --stage neo_prepare --outdir neo_export --threads ${task.cpus}
    """
    stub:
    """
    mkdir neo_export
    echo '{"status":"PREPARED","synthetic":true,"resources":{"cpus":1,"memory_gib":1,"minutes":5}}' > neo_export/task_metadata.json
    echo '{"status":"PREPARED","synthetic":true}' > neo_export/status.json
    """
}

process NEO_PREDICT {
    tag "${meta.case_id}"
    cpus { allowance.cpus }
    memory { allowance.memory_gib.GB }
    time { allowance.minutes.min }
    errorStrategy 'finish'
    publishDir "${params.outdir}/${meta.case_id}", mode: 'copy'
    input:
    tuple val(meta), path(prepared), val(allowance), path(plan)
    output:
    tuple val(meta), path('rna_neoantigens'), emit: results
    script:
    """
    python3 '${projectDir}/extension_tasks.py' --plan '${plan}' --case-id '${meta.case_id}' \
      --stage neo_predict --outdir rna_neoantigens --threads ${task.cpus} --prepared '${prepared}'
    """
    stub:
    """
    mkdir rna_neoantigens
    echo '{"status":"STUB_COMPLETED","synthetic":true}' > rna_neoantigens/status.json
    """
}

workflow OPTIONAL_EXTENSIONS {
    if (!params.plan || !params.outdir) error 'Supply --plan and --outdir'
    def planFile = file(params.plan, checkIfExists: true)
    def plan = new groovy.json.JsonSlurper().parse(planFile.toFile())
    def cases = Channel.fromList(plan.cases)
    BLOOD_CH(cases.filter { it.blood_ch.status == 'READY' }.map { tuple(it, planFile) })
    NEO_PREPARE(cases.filter { it.rna_neoantigens.status == 'READY' }.map { tuple(it, planFile) })
    def predictions = NEO_PREPARE.out.prepared.map { meta, prepared ->
        def taskMeta = new groovy.json.JsonSlurper().parse(prepared.resolve('task_metadata.json').toFile())
        tuple(meta, prepared, taskMeta)
    }.filter { meta, prepared, taskMeta -> taskMeta.status == 'PREPARED' }
     .map { meta, prepared, taskMeta -> tuple(meta, prepared, taskMeta.resources, planFile) }
    NEO_PREDICT(predictions)
}
