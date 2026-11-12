import java.io.File;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/** Exercise the actual SAGE decoder, not a locally invented RC_INFO parser. */
public class NativeContextCheck {
    public static void main(String[] args) throws Exception {
        int count=0;
        try (var fasta=new IndexedFastaSequenceFile(new File(args[1]));
             var reader=new VcfFileReader(args[0])) {
            var ref=new RefGenomeSource(fasta);
            for (var variant : reader.iterator()) {
                var region=new ChrBaseRegion(variant.getContig(),
                    Math.max(1,variant.getStart()-500),variant.getEnd()+500);
                var candidate=CandidateSerialisation.toCandidate(variant,new RefSequence(region,ref));
                if(candidate==null || !candidate.readContext().isValid())
                    throw new IllegalStateException("Invalid context "+variant.getContig()+":"+variant.getStart());
                count++;
            }
        }
        System.out.println("native_context_roundtrips="+count);
    }
}
