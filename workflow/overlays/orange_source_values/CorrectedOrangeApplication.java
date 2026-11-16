package org.qlmp.orange;

import com.google.gson.*;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.report.ReportWriterFactory;
import java.io.*;
import java.nio.file.*;

/** Native ORANGE algorithms and writer, with two source-value mappings. */
public final class CorrectedOrangeApplication {
    public static void main(String[] args) throws Exception {
        ConfigBuilder builder = new ConfigBuilder("Orange");
        OrangeConfig.registerConfig(builder);
        builder.checkAndParseCommandLine(args);
        OrangeConfig config = new OrangeConfig(builder);
        Files.createDirectories(Path.of(config.OutputDir));
        String prefix = config.TumorId + (config.OutputId == null ? "" : "." + config.OutputId);
        Path source = Path.of(config.OutputDir, prefix + ".orange.source.json");
        OrangeRecord record = OrangeAlgo.fromConfig(config).run(config);
        OrangeJson.getInstance().write(record, source);
        JsonObject data;
        try (Reader reader = Files.newBufferedReader(source)) {
            data = JsonParser.parseReader(reader).getAsJsonObject();
        }
        JsonObject receipt = SourceValues.apply(data, Path.of(config.PurpleDataDirectory), config.TumorId);
        Gson gson = new GsonBuilder().serializeNulls().setPrettyPrinting().create();
        Files.writeString(Path.of(config.OutputDir, prefix + ".source_values.json"), gson.toJson(receipt) + "\n");
        OrangeRecord corrected = OrangeJson.getInstance().read(new StringReader(gson.toJson(data)));
        ReportWriterFactory.createToDiskWriter(config).write(corrected);
    }
}
