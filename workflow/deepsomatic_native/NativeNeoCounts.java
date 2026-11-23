import java.io.PrintWriter;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.neo.score.SomaticVariant;

/** Call the installed NEO 1.3 reader on genuine tumour and RNA measurements.
 * This is a consumer compatibility test, not neoantigen prediction. */
public class NativeNeoCounts {
    public static void main(String[] args) throws Exception {
        int count=0;
        try (var reader=new VcfFileReader(args[0]); var out=new PrintWriter(args[3])) {
            out.println("chrom\tpos\tref\talt\ttumor_fragments\ttumor_depth\trna_fragments\trna_depth");
            for(var context : reader.iterator()) {
                var call=SomaticVariant.fromContext(context,args[1],args[2]);
                if(call==null)throw new IllegalStateException("Missing tumour/RNA genotype");
                if(call.TumorFragments<0 || call.TumorDepth<0 || call.RnaFragments<0 || call.RnaDepth<0)
                    throw new IllegalStateException("Negative measurement");
                out.printf("%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d%n",call.Chromosome,call.Position,
                    call.Ref,call.Alt,call.TumorFragments,call.TumorDepth,call.RnaFragments,call.RnaDepth);
                count++;
            }
        }
        System.out.println("native_neo_records="+count);
    }
}
