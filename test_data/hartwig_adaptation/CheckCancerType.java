import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;

/** Probe the installed NEO cache, using a single public reference transcript. */
public class CheckCancerType {
    public static void main(String[] args) throws Exception {
        String transcript = Files.readAllLines(Path.of(args[0])).get(1).split(",", -1)[2];
        Class<?> type = Class.forName("com.hartwig.hmftools.neo.score.TpmMediansCache");
        Object cache = type.getConstructor(String.class).newInstance(args[0]);
        Method has = type.getMethod("hasCancerType", String.class);
        Method get = type.getMethod("getTranscriptTpm", String.class, String.class);
        boolean carcinoid = (boolean) has.invoke(cache, "Carcinoid");
        boolean smallIntestine = (boolean) has.invoke(cache, "Small intestine");
        boolean underscored = (boolean) has.invoke(cache, "Small_intestine");
        double[] unknown = (double[]) get.invoke(cache, transcript, "Carcinoid");
        double[] omitted = (double[]) get.invoke(cache, transcript, "");
        double[] exact = (double[]) get.invoke(cache, transcript, "Small intestine");
        if (carcinoid || !smallIntestine || underscored
                || unknown[0] >= 0 || unknown[1] < 0 || exact[0] < 0
                || unknown[0] != omitted[0] || unknown[1] != omitted[1]) {
            throw new AssertionError("Native NEO terminology/fallback behavior differs from the recorded contract");
        }
        System.out.printf("{\"Carcinoid_match\":%s,\"Small intestine_match\":%s,"
                + "\"Small_intestine_match\":%s,\"unmatched_cancer_tpm\":%s,"
                + "\"pan_cancer_tpm_available\":true,\"unmatched_equals_omitted\":true}%n",
                carcinoid, smallIntestine, underscored, unknown[0]);
    }
}
