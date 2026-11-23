import java.io.File;
import java.util.*;
import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.select.TierSelector;

/** Use SAGE's own ReferenceData loader, without constructing/running a caller.
 * Ported from DNA adapter 7093033; native libraries remain SAGE 5.0.2.
 */
public class NativeSageRegions {
    public static void main(String[] args) throws Exception {
        if(args.length<3)throw new IllegalArgumentException("input.vcf output.vcf native-SAGE-reference-options...");
        ConfigBuilder config=new ConfigBuilder("Sage");
        SageCallConfig.registerConfig(config);
        config.checkAndParseCommandLine(Arrays.copyOfRange(args,2,args.length));
        SageCallConfig sage=new SageCallConfig("5.0.2",config);
        ReferenceData reference=new ReferenceData(sage,config);
        if(!reference.load())throw new IllegalArgumentException("Native SAGE reference loading failed");
        Map<String,TierSelector> selectors=new HashMap<>();
        Map<String,Integer> counts=new TreeMap<>();
        try(VCFFileReader in=new VCFFileReader(new File(args[0]),false)) {
            VCFHeader header=new VCFHeader(in.getFileHeader());
            header.addMetaDataLine(new VCFInfoHeaderLine("TIER",1,VCFHeaderLineType.String,"Native SAGE ReferenceData plus TierSelector; DeepSomatic remains the caller"));
            header.addMetaDataLine(new VCFInfoHeaderLine("HOTSPOT",0,VCFHeaderLineType.Flag,"Exact match to native SAGE-loaded hotspot"));
            header.addMetaDataLine(new VCFInfoHeaderLine("PA_ORIG_FILTER",1,VCFHeaderLineType.String,"Original caller FILTER"));
            header.addMetaDataLine(new VCFHeaderLine("purpleAdapterRegionMethod","native-SAGE5.0.2-ReferenceData-load-and-TierSelector"));
            try(VariantContextWriter out=new VariantContextWriterBuilder().setOutputFile(args[1]).clearOptions().build()) {
                out.writeHeader(header);
                for(VariantContext v:in) {
                    if(v.getAlternateAlleles().size()!=1)throw new IllegalArgumentException("Explicit biallelic projection required");
                    TierSelector selector=selectors.computeIfAbsent(v.getContig(),chr->{
                        if(!HumanChromosome.contains(chr))return new TierSelector(List.of(),List.of(),List.of());
                        var c=HumanChromosome.fromString(chr);
                        return new TierSelector(reference.Hotspots.get(c),reference.PanelWithHotspots.getOrDefault(c,List.of()),reference.HighConfidence.getOrDefault(c,List.of()));
                    });
                    SimpleVariant allele=new SimpleVariant(v.getContig(),v.getStart(),v.getReference().getDisplayString(),v.getAlternateAllele(0).getDisplayString());
                    String tier=selector.tier(allele).toString();counts.merge(tier,1,Integer::sum);
                    VariantContextBuilder b=new VariantContextBuilder(v).attribute("TIER",tier)
                        .attribute("PA_ORIG_FILTER",v.isFiltered()?String.join("&",v.getFilters()):(v.filtersWereApplied()?"PASS":"."));
                    if(selector.isHotspot(allele))b.attribute("HOTSPOT",true);
                    out.add(b.make());
                }
            }
        } finally {reference.RefGenome.close();}
        System.out.println("Native SAGE ReferenceData tier counts: "+counts);
    }
}
