nextflow.enable.dsl=2

process BAM_READY {
    input:
    val meta
    output:
    tuple val(meta), path('bamtools_SYNTHETIC'), topic: bamtools_metrics_dir
    script:
    '''
    mkdir bamtools_SYNTHETIC
    printf 'synthetic completed output\n' > bamtools_SYNTHETIC/event.txt
    '''
}

process CALLER {
    input:
    val start
    output:
    path 'caller_completed.txt'
    script:
    '''
    date -u +%s > caller_started.txt
    sleep 20
    date -u +%s > caller_completed.txt
    '''
}

workflow {
    BAM_READY(channel.value([sample_id: 'SYNTHETIC']))
    CALLER(channel.value(true))
    channel.topic('bamtools_metrics_dir').subscribe { record ->
        LabQcEvents.bamReady(record[1],log)
    }
}
