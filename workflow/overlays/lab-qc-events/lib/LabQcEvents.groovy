import groovy.json.JsonSlurper

/** Operational notification only: subscribe to completed native task outputs. */
class LabQcEvents {
    static void bamReady(Object metrics, Object log) {
        def configPath = System.getenv('WIGITS_LAB_QC_CONFIG')
        if (!configPath) return
        try {
            def config = new JsonSlurper().parse(new File(configPath))
            def builder = new ProcessBuilder([
                config.python.toString(), config.dispatcher.toString(),
                '--config', configPath, '--bam-ready', metrics.toString()
            ]).redirectErrorStream(true)
            builder.environment().remove('PYTHONPATH')
            builder.environment().put('PYTHONNOUSERSITE', '1')
            def process = builder.start()
            def text = process.inputStream.text
            def status = process.waitFor()
            if (status != 0) log.warn("Laboratory QC BAM event failed (exit ${status}): ${text}")
        } catch (Exception error) {
            log.warn("Laboratory QC BAM event failed: ${error.message}")
        }
    }
}
