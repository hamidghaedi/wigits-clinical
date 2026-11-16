package org.qlmp.orange;

import com.google.gson.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

public final class SourceValuesTest {
    static int checks = 0;
    static void equal(Object actual, Object expected) {
        ++checks;
        if (!Objects.equals(actual, expected)) throw new AssertionError(actual + " != " + expected);
    }
    static Map<String, String> segment(String chromosome, int start, int end, double cn) {
        return Map.of("chromosome", chromosome, "start", "" + start, "end", "" + end, "copyNumber", "" + cn);
    }
    static String range(Map<String, String>... segments) {
        return SourceValues.extent(Map.of("chromosome", "chr1", "start", "10", "end", "20"), Arrays.asList(segments)).get("range").getAsString();
    }
    public static void main(String[] args) throws Exception {
        equal(range(segment("1", 1, 30, 0.1)), "FULL");
        equal(range(segment("chr1", 1, 14, 0.1), segment("1", 15, 30, 2)), "PARTIAL");
        equal(range(segment("1", 1, 14, 0.1), segment("1", 15, 30, 0.4)), "FULL");
        equal(range(segment("1", 10, 20, 0.4999)), "FULL");
        equal(range(segment("1", 10, 20, 0.5)), "N/A");
        equal(range(segment("1", 1, 18, 0.1)), "N/A");
        equal(range(segment("1", 1, 13, 0.1), segment("1", 15, 30, 0.1)), "N/A");
        equal(range(segment("1", 1, 15, 0.1), segment("1", 15, 30, 0.1)), "N/A");
        equal(range(segment("2", 1, 30, 0.1)), "N/A");
        Path directory = Files.createTempDirectory("source-values-test");
        Path vcf = directory.resolve("T.purple.somatic.vcf.gz");
        Path qc = directory.resolve("T.purple.qc");
        Files.writeString(qc, "TincLevel\t0.0000\n");
        for (String header : List.of("", "0", "0.02504", "NaN", "0.2\n##tincLevel=0.3")) {
            Files.writeString(qc, "TincLevel\t" + (header.equals("0.02504") ? "0.0250" : "0.0000") + "\n");
            try (GZIPOutputStream out = new GZIPOutputStream(Files.newOutputStream(vcf))) {
                out.write(("##fileformat=VCFv4.2\n" + (header.isEmpty() ? "" : "##tincLevel=" + header + "\n") + "#CHROM\tPOS\n").getBytes());
            }
            JsonObject result = SourceValues.tinc(vcf, qc);
            equal(result.get("status").getAsString(), header.equals("0") || header.equals("0.02504") ? "available" : "unavailable");
            if (header.isEmpty()) equal(result.get("value").isJsonNull(), true);
            if (header.equals("0")) equal(result.get("value").getAsDouble(), 0.0);
        }
        Files.writeString(directory.resolve("T.purple.cnv.gene.tsv"), "chromosome\tstart\tend\tgene\ttranscriptId\n1\t10\t20\tA\tTX1\n");
        Files.writeString(directory.resolve("T.purple.cnv.somatic.tsv"), "chromosome\tstart\tend\tcopyNumber\n1\t1\t30\t0.1\n");
        JsonObject record = JsonParser.parseString("{\"purple\":{\"somaticGainsDels\":["
                + "{\"driver\":{\"gene\":\"A\",\"transcript\":\"TX1\",\"type\":\"DEL\"},\"chromosome\":\"chr1\",\"geneRange\":\"PARTIAL\"},"
                + "{\"driver\":{\"gene\":\"A\",\"transcript\":\"TX2\",\"type\":\"DEL\"},\"chromosome\":\"1\",\"geneRange\":\"PARTIAL\"},"
                + "{\"driver\":{\"gene\":\"A\",\"transcript\":\"TX1\",\"type\":\"AMP\"},\"chromosome\":\"1\",\"geneRange\":\"FULL\"}]}}").getAsJsonObject();
        SourceValues.apply(record, directory, "T");
        JsonArray calls = record.getAsJsonObject("purple").getAsJsonArray("somaticGainsDels");
        equal(calls.get(0).getAsJsonObject().get("geneRange").getAsString(), "FULL");
        equal(calls.get(1).getAsJsonObject().get("geneRange").getAsString(), "N/A");
        equal(calls.get(2).getAsJsonObject().get("geneRange").getAsString(), "FULL");
        try (var paths = Files.walk(directory)) {
            for (Path path : paths.sorted(Comparator.reverseOrder()).toList()) Files.delete(path);
        }
        System.out.println("SourceValues tests passed: " + checks);
    }
}
