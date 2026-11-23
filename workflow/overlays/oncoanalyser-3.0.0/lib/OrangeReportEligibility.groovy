/** ORANGE eligibility preserves native QC failures without inventing plots. */
class OrangeReportEligibility {
    static final Set<String> QC_FAILURES = ['FAIL_NO_TUMOR', 'FAIL_CONTAMINATION', 'FAIL_TINC'] as Set

    static boolean shouldRun(boolean tumourDna, boolean requiredInputs,
                             boolean somaticVcf, boolean circosPlot,
                             Collection<String> qcStatuses) {
        if (!tumourDna) return false
        if (!requiredInputs || !somaticVcf) {
            throw new IllegalStateException('ORANGE requested but required analytical inputs are missing')
        }
        if (circosPlot) return true
        if (qcStatuses.any { it in QC_FAILURES }) return true
        throw new IllegalStateException('ORANGE requested but the PURPLE circos plot is missing without a supported QC-fail status')
    }
}
