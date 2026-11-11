import java.io.*;
import java.util.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;

/** Reuses VariantReadContextBuilder.setMaxRefRepeat's search, without reads.
 * Ported from DNA adapter 7093033; TNC uses the same reference-centred three
 * bases as SAGE 5.0.2 VariantReadContext.trinucleotide(). No read context/MH
 * or phase information is inferred from reference sequence.
 */
public class ReferenceRepeatAnnotator {
    static RepeatInfo repeat(byte[] bases,int index) {
        return RepeatInfo.findMaxRepeat(bases,index,index,MAX_REPEAT_LENGTH,MIN_REPEAT_COUNT+1,false,index);
    }
    public static void main(String[] args) throws Exception {
        if(args.length==1 && args[0].equals("--self-test")) {
            RepeatInfo a=repeat("AAAAAAAAACG".getBytes(),0);
            RepeatInfo b=repeat("ATATATATATATATATCG".getBytes(),0);
            if(a==null || a.Count!=9 || !a.Bases.equals("A")) throw new AssertionError("Homopolymer");
            if(b==null || b.Count!=8 || !b.Bases.equals("AT")) throw new AssertionError("Dinucleotide repeat");
            if(repeat("ACGTGCTACGATC".getBytes(),0)!=null) throw new AssertionError("Nonrepeat");
            System.out.println("Native repeat tests passed; max unit="+MAX_REPEAT_LENGTH+" min count="+(MIN_REPEAT_COUNT+1));return;
        }
        if(args.length!=3) throw new IllegalArgumentException("input.vcf output.vcf reference.fa");
        Map<String,Integer> lengths=new HashMap<>();
        try(BufferedReader fai=new BufferedReader(new FileReader(args[2]+".fai"))) {
            for(String line;(line=fai.readLine())!=null;) {String[] f=line.split("\t");lengths.put(f[0],Integer.parseInt(f[1]));}
        }
        int repeats=0,none=0;
        try(IndexedFastaSequenceFile ref=new IndexedFastaSequenceFile(new File(args[2]));VCFFileReader in=new VCFFileReader(new File(args[0]),false)) {
            VCFHeader h=new VCFHeader(in.getFileHeader());
            h.addMetaDataLine(new VCFInfoHeaderLine("REP_C",1,VCFHeaderLineType.Integer,"Native Hartwig qualifying reference repeat count;0 means measured search found none"));
            h.addMetaDataLine(new VCFInfoHeaderLine("REP_S",1,VCFHeaderLineType.String,"Native Hartwig qualifying reference repeat unit"));
            h.addMetaDataLine(new VCFInfoHeaderLine("PA_REPEAT_STATUS",1,VCFHeaderLineType.String,"Reference-only native search; no read-context repeat or phasing inferred"));
            h.addMetaDataLine(new VCFInfoHeaderLine("TNC",1,VCFHeaderLineType.String,"Three reference bases centred on variant POS, as native SAGE VariantReadContext.trinucleotide"));
            h.addMetaDataLine(new VCFInfoHeaderLine("PA_TNC_STATUS",1,VCFHeaderLineType.String,"Explicit reference context availability"));
            h.addMetaDataLine(new VCFHeaderLine("purpleAdapterRepeatMethod","SAGE5.0.2-RepeatInfo.findMaxRepeat-native-parameters-reference-only"));
            try(VariantContextWriter out=new VariantContextWriterBuilder().setOutputFile(args[1]).clearOptions().build()) {
                out.writeHeader(h);
                for(VariantContext v:in) {
                    String bases=v.getReference().getDisplayString();
                    String observed=ref.getSubsequenceAt(v.getContig(),v.getStart(),v.getStart()+bases.length()-1).getBaseString();
                    if(!observed.equalsIgnoreCase(bases)) throw new IllegalArgumentException("REF mismatch at "+v.getContig()+":"+v.getStart());
                    int target=v.getStart()+(v.isIndel()?1:0), length=lengths.get(v.getContig());
                    if(target>length) throw new IllegalArgumentException("No reference context after contig end");
                    int window=256;RepeatInfo found;
                    while(true) {
                        int end=Math.min(length,target+window-1);
                        byte[] sequence=ref.getSubsequenceAt(v.getContig(),target,end).getBases();
                        for(int i=0;i<sequence.length;i++)sequence[i]=(byte)Character.toUpperCase((char)sequence[i]);
                        found=repeat(sequence,0);
                        if(found==null || found.endIndex()<sequence.length-MAX_REPEAT_LENGTH || end==length)break;
                        if(window>=65536)throw new IllegalArgumentException("Repeat exceeds bounded reference search");
                        window*=2;
                    }
                    VariantContextBuilder b=new VariantContextBuilder(v).attribute("REP_C",found==null?0:found.Count)
                        .attribute("PA_REPEAT_STATUS",found==null?"MEASURED_NO_QUALIFYING_REPEAT":"MEASURED_REFERENCE_REPEAT");
                    if(v.getStart()>1 && v.getStart()<length) {
                        String tnc=ref.getSubsequenceAt(v.getContig(),v.getStart()-1,v.getStart()+1).getBaseString().toUpperCase();
                        b.attribute("TNC",tnc).attribute("PA_TNC_STATUS",tnc.matches("[ACGT]{3}")?"MEASURED":"AMBIGUOUS_REFERENCE_BASE");
                    } else {
                        b.rmAttribute("TNC").attribute("PA_TNC_STATUS","CONTIG_BOUNDARY");
                    }
                    if(found!=null){b.attribute("REP_S",found.Bases);repeats++;}else{b.rmAttribute("REP_S");none++;}
                    out.add(b.make());
                }
            }
        }
        System.out.println("Reference repeats="+repeats+" measured no qualifying repeat="+none);
    }
}
