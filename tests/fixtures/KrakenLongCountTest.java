import au.edu.wehi.idsv.kraken.KrakenReportLine;

public class KrakenLongCountTest {
    public static void main(String[] args) {
        for (long count : new long[]{0, 50, 2147483647L, 2147483648L, 3576733308L, 5000000000L}) {
            String line = "93.27\t" + count + "\t" + count + "\tU\t0\tunclassified";
            KrakenReportLine r = new KrakenReportLine(line);
            if (r.countAssignedToTree != count || r.countAssignedDirectly != count
                    || r.taxonomyId != 0 || !r.line.equals(line)) {
                throw new AssertionError("Count or original report text changed");
            }
        }
        System.out.println("Kraken long-count regression passed");
    }
}
