import java.io.File;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

/** Obtain actual native MH/repeat/TNC annotations by native context roundtrip. */
public class EnrichSeedContext {
    public static void main(String[] args) throws Exception {
        try(var fasta=new IndexedFastaSequenceFile(new File(args[1]));
            var input=new VCFFileReader(new File(args[0]),false);
            var output=new VariantContextWriterBuilder().setOutputFile(args[2]).clearOptions().build()) {
            output.writeHeader(input.getFileHeader());
            var ref=new RefGenomeSource(fasta);
            for(var v:input) {
                var seq=new RefSequence(new ChrBaseRegion(v.getContig(),Math.max(1,v.getStart()-1000),v.getEnd()+1000),ref);
                var c=CandidateSerialisation.toCandidate(v,seq);
                if(c==null || !c.readContext().isValid())throw new IllegalStateException("Invalid observed seed "+v);
                var b=new VariantContextBuilder(v);
                CandidateSerialisation.toContext(c).make().getAttributes().forEach((k,value)->b.attribute(k,value));
                output.add(b.make());
            }
        }
    }
}
