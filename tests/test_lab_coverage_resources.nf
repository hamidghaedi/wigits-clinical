nextflow.enable.dsl=2

workflow {
    def summary = file("${launchDir}/synthetic-summary.tsv")
    [0:64, 23:64, 50:64, 51:96, 100:96, 101:128, 141:128, 151:160].each { depth, ram ->
        summary.text = "TotalReads\tMeanCoverage\n100\t${depth}\n"
        def result = LabQcResources.coverage(182.4, summary)
        assert result.memory_gib == ram
        assert result.minutes == 183
        assert result.mean_coverage == depth
    }
    assert LabQcResources.coverage(30, null).memory_gib == 128
    assert LabQcResources.coverage(30, '/missing/metadata').memory_gib == 128
    summary.text = 'MeanCoverage\nNaN\n'
    assert LabQcResources.coverage(30, summary).memory_gib == 128
    summary.text = 'MeanCoverage\n23\n'
    assert LabQcResources.coverage(1024, summary).memory_gib == 96
    assert LabQcResources.coverage(30, summary).minutes == 60
    summary.delete()
    log.info 'Coverage depth/size allowances: all boundary, missing and malformed cases passed'
}
