/** Metadata-only initial allowances; never scan BAMs or change coverage filters. */
class LabQcResources {
    static Map coverage(Number bamGiB, Object summaryPath = null) {
        Double depth = null
        if (summaryPath) {
            File source = new File(summaryPath.toString())
            if (source.isFile()) {
                def lines = source.readLines().findAll { it && !it.startsWith('#') }
                if (lines.size() >= 2) {
                    int column = lines[0].split('\t', -1).findIndexOf { it == 'MeanCoverage' }
                    def cells = lines[1].split('\t', -1)
                    if (column >= 0 && column < cells.size()) {
                        try {
                            double measured = Double.parseDouble(cells[column])
                            if (Double.isFinite(measured) && measured >= 0) depth = measured
                        } catch (NumberFormatException ignored) { }
                    }
                }
            }
        }
        // 64 GiB failed on a 141x tumour; a 23x normal used 24 GiB.
        // Depth is already measured by native BAMTOOLS. Missing metadata uses
        // the conservative 128 GiB floor, rather than requiring another scan.
        int depthGiB = depth == null ? 128 : 64 + 32 * Math.max(0, Math.ceil((depth - 50) / 50).toInteger())
        int sizeExtra = 32 * Math.max(0, Math.ceil(bamGiB.doubleValue() / 512).toInteger() - 1)
        [memory_gib: depthGiB + sizeExtra,
         minutes: Math.max(60, Math.ceil(bamGiB.doubleValue())).toInteger(),
         mean_coverage: depth, policy: 'lab-coverage-depth-v1']
    }
}
