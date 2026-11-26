import groovy.json.JsonSlurper

/** Deliver this completed report from its task outputs, without waiting for other runs. */
class LabQcDelivery {
    static void reportReady(Object success, Object log) {
        def configPath=System.getenv('WIGITS_LAB_QC_CONFIG')
        if (!configPath) return
        try {
            def config=new JsonSlurper().parse(new File(configPath))
            def builder=new ProcessBuilder([config.python.toString(),config.dispatcher.toString(),
                '--config',configPath,'--report-ready',success.toString()]).redirectErrorStream(true)
            builder.environment().remove('PYTHONPATH')
            builder.environment().put('PYTHONNOUSERSITE','1')
            def process=builder.start();def output=process.inputStream.text;def code=process.waitFor()
            if (code!=0) log.warn("Laboratory QC delivery event failed (exit ${code}): ${output}")
        } catch (Exception error) {
            log.warn("Laboratory QC delivery event failed: ${error.message}")
        }
    }
}
