import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.candidate.RefContext;
import com.hartwig.hmftools.sage.vcf.ReadContextVcfInfo;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.NumberEvents.calcAdjustedNumMutations;

/** Diagnose missing indels in observed primary reads. Haplotype substrings are
 * diagnostic only; seed export requires the exact observed CIGAR event.
 * No altered qualities, constructed reads or invented allele-support vectors. */
public class InspectMissingContexts {
    public static void main(String[] args) throws Exception {
        try(var fasta=new IndexedFastaSequenceFile(new File(args[1]));
            // Google diagnostic BAM serialization can omit paired-read flags.
            // Permit inspection explicitly; never use these BAMs for counting.
            var bam=SamReaderFactory.makeDefault().validationStringency(
                args.length>4 && args[4].equals("--diagnostic-bam") ? ValidationStringency.LENIENT : ValidationStringency.STRICT)
                .open(new File(args[0]));
            var out=new PrintWriter(args[3]);
            var primary=new PrintWriter(args[3]+".primary.tsv")) {
            var builder=new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);
            var ref=new RefGenomeSource(fasta);
            out.println("requested\tobserved\tmapq\tcontext\tbase_qualities\tread_name\tnative_events\tmethod");
            primary.println("requested\tcontext\tfull_matches\tcore_matches\tnative_events\trole");
            for(var line:Files.readAllLines(Path.of(args[2])).subList(1,Files.readAllLines(Path.of(args[2])).size())) {
                var x=line.split(",");String chr=x[0];int pos=Integer.parseInt(x[1]);
                var variant=new SimpleVariant(chr,pos,x[2],x[3]);
                var observed=new AltContext(new RefContext(chr,pos),x[2],x[3]);
                var seq=new RefSequence(new ChrBaseRegion(chr,Math.max(1,pos-2000),pos+x[2].length()+2000),ref);
                try(var it=bam.queryOverlapping(chr,Math.max(1,pos-1000),pos+x[2].length()+1000)) {
                    while(it.hasNext()) {
                        var read=it.next();
                        if(read.isSecondaryOrSupplementary() || read.getDuplicateReadFlag() || read.getReadFailsVendorQualityCheckFlag())continue;
                        for(int flank:new int[]{10,12,15,20,30,40}) {
                            String left=new String(fasta.getSubsequenceAt(chr,pos-flank,pos-1).getBases());
                            String right=new String(fasta.getSubsequenceAt(chr,pos+x[2].length(),pos+x[2].length()+flank-1).getBases());
                            String hap=left+x[3]+right;
                            String referenceWindow=new String(fasta.getSubsequenceAt(chr,
                                Math.max(1,pos-1000),pos+x[2].length()+1000).getBases());
                            if(referenceWindow.contains(hap))continue;
                            int index=read.getReadString().indexOf(hap);
                            if(index<0 || read.getReadString().indexOf(hap,index+1)>=0)continue;
                            var rc=builder.createContext(variant,read,index+flank,seq);
                            String context=rc==null?"NATIVE_BUILDER_REJECTED":
                                (rc.isValid()?new ReadContextVcfInfo(rc).toVcfTag():"INVALID");
                            out.printf("%s:%s:%s:%s\t%s:%s:%s:%s\t%d\t%s\t%s\t%s\t%d\tALT_HAPLOTYPE_%d%n",
                                chr,x[1],x[2],x[3],chr,x[1],x[2],x[3],read.getMappingQuality(),context,
                                Arrays.toString(read.getBaseQualities()),read.getReadName(),calcAdjustedNumMutations(read,seq),flank);
                        }
                        int rp=read.getAlignmentStart(),ri=0;
                        for(var ce:read.getCigar().getCigarElements()) {
                            var op=ce.getOperator();int len=ce.getLength();String r=null,a=null;
                            if(op==CigarOperator.I && ri>0) {
                                r=new String(fasta.getSubsequenceAt(chr,rp-1,rp-1).getBases());
                                a=read.getReadString().substring(ri-1,ri+len);
                            } else if(op==CigarOperator.D && ri>0) {
                                r=new String(fasta.getSubsequenceAt(chr,rp-1,rp+len-1).getBases());
                                a=r.substring(0,1);
                            }
                            if(r!=null && Math.abs(rp-1-pos)<100) {
                                String context="NOT_EXACT";
                                if(rp-1==pos && r.equals(x[2]) && a.equals(x[3])) {
                                    observed.addReadContext(calcAdjustedNumMutations(read,seq),read,ri-1,builder,seq);
                                    var rc=builder.createContext(variant,read,ri-1,seq);
                                    context=rc==null?"NATIVE_BUILDER_REJECTED":
                                        (rc.isValid()?new ReadContextVcfInfo(rc).toVcfTag():"INVALID");
                                }
                                out.printf("%s:%s:%s:%s\t%s:%d:%s:%s\t%d\t%s\t%s\t%s\t%d\tCIGAR%n",chr,x[1],x[2],x[3],
                                    chr,rp-1,r,a,read.getMappingQuality(),context,Arrays.toString(read.getBaseQualities()),
                                    read.getReadName(),calcAdjustedNumMutations(read,seq));
                            }
                            if(op.consumesReadBases())ri+=len;
                            if(op.consumesReferenceBases())rp+=len;
                        }
                    }
                }
                observed.selectCandidates();
                if(observed.hasValidCandidate()) {
                    var c=observed.candidate();
                    primary.printf("%s:%s:%s:%s\t%s\t%d\t%d\t%d\tPRIMARY%n",chr,x[1],x[2],x[3],
                        new ReadContextVcfInfo(c.readContext()).toVcfTag(),c.FullMatch,c.CoreMatch,c.MinNumberOfEvents);
                    if(observed.hasSecondCandidate()) {
                        c=observed.secondCandidate();
                        primary.printf("%s:%s:%s:%s\t%s\t%d\t%d\t%d\tSECONDARY%n",chr,x[1],x[2],x[3],
                            new ReadContextVcfInfo(c.readContext()).toVcfTag(),c.FullMatch,c.CoreMatch,c.MinNumberOfEvents);
                    }
                }
            }
        }
    }
}
