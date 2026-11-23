package org.qlmp.orange;

import com.google.gson.*;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/** Source-value corrections only. PURPLE 4.4 DeletionDrivers uses CN < 0.5 for DEL.
 * Coordinates in both PURPLE tables are 1-based inclusive. No exon inference,
 * variant selection, purity selection, QC or template changes are performed.
 */
public final class SourceValues {
    public static final String VERSION = "0.1.0";

    static List<Map<String, String>> table(Path path) throws IOException {
        List<Map<String, String>> rows = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(path)) {
            String header = reader.readLine();
            if (header == null) throw new IOException("Empty PURPLE table: " + path);
            String[] names = header.split("\t", -1);
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) continue;
                String[] cells = line.split("\t", -1);
                if (cells.length != names.length) throw new IOException("Invalid PURPLE row: " + path);
                Map<String, String> row = new HashMap<>();
                for (int i = 0; i < names.length; ++i) row.put(names[i], cells[i]);
                rows.add(row);
            }
        }
        return rows;
    }

    public static JsonObject tinc(Path vcf, Path qc) throws IOException {
        JsonObject result = new JsonObject();
        result.add("value", JsonNull.INSTANCE);
        result.addProperty("status", "unavailable");
        result.addProperty("source", vcf.toString());
        result.add("native_qc_value", JsonNull.INSTANCE);
        for (String line : Files.readAllLines(qc)) {
            String[] cells = line.split("\t", -1);
            if (cells.length == 2 && cells[0].equals("TincLevel")) result.addProperty("native_qc_value", cells[1]);
        }
        List<String> headers = new ArrayList<>();
        InputStream stream = Files.newInputStream(vcf);
        if (vcf.toString().endsWith(".gz")) stream = new GZIPInputStream(stream);
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null && line.startsWith("#")) {
                if (line.startsWith("##tincLevel=")) headers.add(line.substring("##tincLevel=".length()));
            }
        }
        result.addProperty("reason", "missing_tincLevel_header");
        if (!headers.isEmpty()) {
            try {
                double value = Double.parseDouble(headers.get(0));
                boolean valid = Double.isFinite(value) && value >= 0 && value <= 1;
                for (String header : headers) valid &= Double.compare(value, Double.parseDouble(header)) == 0;
                if (valid && !result.get("native_qc_value").isJsonNull()) {
                    double nativeValue = Double.parseDouble(result.get("native_qc_value").getAsString());
                    valid = Double.isFinite(nativeValue) && Math.abs(value - nativeValue) <= 0.000050001;
                }
                if (valid) {
                    result.addProperty("value", value);
                    result.addProperty("status", "available");
                    result.addProperty("reason", "vcf_tincLevel_header");
                } else result.addProperty("reason", "invalid_or_conflicting_tincLevel");
            } catch (NumberFormatException e) {
                result.addProperty("reason", "invalid_or_conflicting_tincLevel");
            }
        }
        return result;
    }

    static String chromosome(String name) {
        return name.startsWith("chr") ? name.substring(3) : name;
    }

    static JsonObject extent(Map<String, String> gene, List<Map<String, String>> segments) {
        int start = Integer.parseInt(gene.get("start"));
        int end = Integer.parseInt(gene.get("end"));
        JsonObject result = new JsonObject();
        result.addProperty("chromosome", gene.get("chromosome"));
        result.addProperty("start", start);
        result.addProperty("end", end);
        result.addProperty("range", "N/A");
        result.addProperty("reason", "incomplete_or_ambiguous_segment_coverage");
        List<Map<String, String>> overlapping = new ArrayList<>();
        for (Map<String, String> segment : segments) {
            if (chromosome(gene.get("chromosome")).equals(chromosome(segment.get("chromosome")))
                    && Integer.parseInt(segment.get("start")) <= end && Integer.parseInt(segment.get("end")) >= start) {
                overlapping.add(segment);
            }
        }
        overlapping.sort(Comparator.comparingInt(s -> Integer.parseInt(s.get("start"))));
        long next = start;
        long deleted = 0;
        JsonArray evidence = new JsonArray();
        for (Map<String, String> segment : overlapping) {
            int left = Math.max(start, Integer.parseInt(segment.get("start")));
            int right = Math.min(end, Integer.parseInt(segment.get("end")));
            double copies = Double.parseDouble(segment.get("copyNumber"));
            if (left != next || !Double.isFinite(copies)) return result;
            JsonObject span = new JsonObject();
            span.addProperty("start", Integer.parseInt(segment.get("start")));
            span.addProperty("end", Integer.parseInt(segment.get("end")));
            span.addProperty("copy_number", copies);
            evidence.add(span);
            if (copies < 0.5) deleted += (long) right - left + 1;
            next = (long) right + 1;
        }
        result.add("segments", evidence);
        if (start > end || next != (long) end + 1 || deleted == 0) return result;
        result.addProperty("range", deleted == (long) end - start + 1 ? "FULL" : "PARTIAL");
        result.addProperty("reason", "gene_interval_covered_by_observed_segments");
        result.addProperty("deleted_bases", deleted);
        return result;
    }

    public static JsonObject apply(JsonObject record, Path purple, String sample) throws IOException {
        List<Map<String, String>> genes = table(purple.resolve(sample + ".purple.cnv.gene.tsv"));
        List<Map<String, String>> segments = table(purple.resolve(sample + ".purple.cnv.somatic.tsv"));
        JsonObject receipt = new JsonObject();
        receipt.addProperty("mapping_version", VERSION);
        receipt.add("tinc", tinc(purple.resolve(sample + ".purple.somatic.vcf.gz"), purple.resolve(sample + ".purple.qc")));
        JsonArray ranges = new JsonArray();
        JsonArray calls = record.getAsJsonObject("purple").getAsJsonArray("somaticGainsDels");
        for (JsonElement element : calls) {
            JsonObject call = element.getAsJsonObject();
            JsonObject driver = call.getAsJsonObject("driver");
            if (!driver.get("type").getAsString().equals("DEL")) continue;
            String name = driver.get("gene").getAsString();
            String transcript = driver.get("transcript").getAsString();
            List<Map<String, String>> matches = new ArrayList<>();
            for (Map<String, String> gene : genes) {
                if (name.equals(gene.get("gene")) && transcript.equals(gene.get("transcriptId"))
                        && chromosome(call.get("chromosome").getAsString()).equals(chromosome(gene.get("chromosome")))) matches.add(gene);
            }
            JsonObject range;
            if (matches.size() == 1) range = extent(matches.get(0), segments);
            else {
                range = new JsonObject();
                range.addProperty("range", "N/A");
                range.addProperty("reason", "missing_or_ambiguous_gene_transcript_interval");
            }
            range.addProperty("gene", name);
            range.addProperty("transcript", transcript);
            range.add("native_range", call.get("geneRange").deepCopy());
            call.addProperty("geneRange", range.get("range").getAsString());
            ranges.add(range);
        }
        receipt.add("deletion_ranges", ranges);
        return receipt;
    }
}
