// Resource directives only; sparse synthetic inputs are never read as sequences.
process STAR_ALIGN {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), val(rg_lines), path(reads_fwds), path(reads_revs)
    output:
    path "${meta.id}.STAR_ALIGN.tsv"
    script:
    """
    printf '%s\t%s\t%s\n' '${task.cpus}' '${task.memory.toGiga()}' '${task.time.toHours()}' > ${meta.id}.STAR_ALIGN.tsv
    """
}
process SAMTOOLS_SORT {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), path(bam)
    output:
    path "${meta.id}.SAMTOOLS_SORT.tsv"
    script:
    """
    printf '%s\t%s\t%s\n' '${task.cpus}' '${task.memory.toGiga()}' '${task.time.toHours()}' > ${meta.id}.SAMTOOLS_SORT.tsv
    """
}
process GATK4_MARKDUPLICATES {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), path(bam)
    output:
    path "${meta.id}.GATK4_MARKDUPLICATES.tsv"
    script:
    """
    printf '%s\t%s\t%s\n' '${task.cpus}' '${task.memory.toGiga()}' '${task.time.toHours()}' > ${meta.id}.GATK4_MARKDUPLICATES.tsv
    """
}
process ISOFOX {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(meta), path(aln), path(idx)
    output:
    path "${meta.id}.ISOFOX.tsv"
    script:
    """
    printf '%s\t%s\t%s\n' '${task.cpus}' '${task.memory.toGiga()}' '${task.time.toHours()}' > ${meta.id}.ISOFOX.tsv
    """
}
workflow READ_ALIGNMENT_RNA {
    STAR_ALIGN(channel.fromList(['S','L']).map { s -> tuple([id:s], [], [file("${params.inputs}/${s}.R1")], [file("${params.inputs}/${s}.R2")]) })
    SAMTOOLS_SORT(channel.fromList(['S','L']).map { s -> tuple([id:s], file("${params.inputs}/${s}.unsorted.bam")) })
    GATK4_MARKDUPLICATES(channel.fromList(['S','L']).map { s -> tuple([id:s], [file("${params.inputs}/${s}.sorted.bam")]) })
}
workflow ISOFOX_QUANTIFICATION {
    ISOFOX(channel.fromList(['S','L']).map { s -> tuple([id:s], file("${params.inputs}/${s}.md.bam"), file("${params.inputs}/${s}.bai")) })
}
workflow WGTS {
    READ_ALIGNMENT_RNA()
    ISOFOX_QUANTIFICATION()
}
workflow { WGTS() }
