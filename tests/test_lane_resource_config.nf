// Metadata-only canary: read file sizes, then echo resolved task directives.
// No input content is read and no analytical tool runs.
process BWAMEM2_ALIGN {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    output:
    path "${meta.sample_id}.lane.tsv"
    script:
    """
    printf '%s\t%s\t%s\t%s\n' '${meta.sample_id}' '${task.cpus}' '${task.memory.toGiga()}' '${task.time.toHours()}' > ${meta.sample_id}.lane.tsv
    """
}

process REDUX {
    publishDir params.outdir, mode: 'copy'
    input:
    val meta
    output:
    path "${meta.sample_id}.sample.tsv"
    script:
    """
    printf '%s\t%s\t%s\t%s\n' '${meta.sample_id}' '${task.cpus}' '${task.memory.toGiga()}' '${task.time.toHours()}' > ${meta.sample_id}.sample.tsv
    """
}

workflow {
    inputs = channel.fromList(['N', 'T']).map { sample ->
        tuple([sample_id: sample], file("${params.input_dir}/${sample}1.fastq.gz"), file("${params.input_dir}/${sample}2.fastq.gz"))
    }
    BWAMEM2_ALIGN(inputs)
    REDUX(channel.fromList([[sample_id: 'N'], [sample_id: 'T']]))
}
