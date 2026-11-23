nextflow.enable.dsl=2

workflow {
    assert OrangeReportEligibility.shouldRun(true, true, true, true, ['PASS'])
    assert OrangeReportEligibility.shouldRun(true, true, true, true, ['WARN_LOW_PURITY'])
    ['FAIL_NO_TUMOR', 'FAIL_CONTAMINATION', 'FAIL_TINC'].each { code ->
        assert OrangeReportEligibility.shouldRun(true, true, true, false, [code])
    }
    assert !OrangeReportEligibility.shouldRun(false, false, false, false, [])
    [['PASS'], ['WARN_LOW_PURITY'], [], ['UNKNOWN']].each { codes ->
        boolean failed = false
        try { OrangeReportEligibility.shouldRun(true, true, true, false, codes) }
        catch (IllegalStateException expected) { failed = true }
        assert failed
    }
    [[false, true], [true, false], [false, false]].each { inputs, vcf ->
        boolean failed = false
        try { OrangeReportEligibility.shouldRun(true, inputs, vcf, false, ['FAIL_NO_TUMOR']) }
        catch (IllegalStateException expected) { failed = true }
        assert failed
    }
    log.info 'ORANGE eligibility: 13 cases passed; QC failure retained; missing required data rejected'
}
