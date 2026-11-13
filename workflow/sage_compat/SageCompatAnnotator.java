import java.io.*;
import java.nio.file.*;
import java.util.*;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import htsjdk.variant.vcf.*;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.common.Microhomology;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.select.TierSelector;

import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;

/**
 * sage-compat adapter annotator (wigits-clinical workflow/sage_compat).
 *
 * Runs with SAGE 5.0.2's own jar on the class path (java --class-path sage.jar SageCompatAnnotator.java ...)
 * and reuses SAGE's native classes: ReferenceData/TierSelector for INFO/TIER, RepeatInfo.findMaxRepeat for
 * REP_C/REP_S (as VariantReadContextBuilder.setMaxRefRepeat, reference-only), Microhomology.findHomology for MH.
 * TNC is the reference trinucleotide centred on POS (VariantReadContext.trinucleotide). Normal AD/DP/AF come
 * from an htsjdk pileup of the normal REDUX BAM. Nothing read-context based (RC_*, LPS, quality fields) is emitted.
 *
 * Modes:
 *   panel-bed <out.bed> <sage reference options...>
 *   annotate <candidates.tsv> <out.vcf.gz> <annotations.tsv> <normal.bam> <normal_id> <tumor_id> <header_extra.txt>
 *            <min_mapq> <min_baseq> <sage reference options...>
 *   --self-test
 */
public class SageCompatAnnotator
{
    static final String VERSION = "0.1.0";

    // ---------------------------------------------------------------- reference-derived annotations
    static RepeatInfo repeatAt(byte[] bases, int index)
    {
        // identical call to VariantReadContextBuilder.setMaxRefRepeat (sage-v5.0.2)
        return RepeatInfo.findMaxRepeat(bases, index, index, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1, false, index);
    }

    static byte[] upper(byte[] b) { for(int i = 0; i < b.length; i++) b[i] = (byte)Character.toUpperCase((char)b[i]); return b; }

    /** Reference-only repeat search starting at the variant base (indel: following base), unbounded to the right. */
    static RepeatInfo refRepeat(IndexedFastaSequenceFile ref, int contigLength, String chr, int pos, boolean indel)
    {
        int target = pos + (indel ? 1 : 0);
        if(target > contigLength) return null;
        int window = 256;
        while(true)
        {
            int end = Math.min(contigLength, target + window - 1);
            byte[] seq = upper(ref.getSubsequenceAt(chr, target, end).getBases());
            RepeatInfo found = repeatAt(seq, 0);
            if(found == null || found.endIndex() < seq.length - MAX_REPEAT_LENGTH || end == contigLength) return found;
            if(window >= 65536) throw new IllegalStateException("repeat exceeds bounded reference search at " + chr + ":" + pos);
            window *= 2;
        }
    }

    static String tnc(IndexedFastaSequenceFile ref, int contigLength, String chr, int pos)
    {
        if(pos - 1 < 1 || pos + 1 > contigLength) return null;
        return new String(upper(ref.getSubsequenceAt(chr, pos - 1, pos + 1).getBases()));
    }

    /** MH as SAGE: longer of findHomology on the ALT haplotype (read) and on the reference, both starting at POS. */
    static String microhomology(IndexedFastaSequenceFile ref, int contigLength, String chr, int pos, String refAllele, String altAllele)
    {
        SimpleVariant v = new SimpleVariant(chr, pos, refAllele, altAllele);
        if(!v.isIndel()) return null;
        int flank = Math.max(200, 4 * Math.max(refAllele.length(), altAllele.length()));
        int end = Math.min(contigLength, pos + refAllele.length() - 1 + flank);
        byte[] refBases = upper(ref.getSubsequenceAt(chr, pos, end).getBases()); // starts at POS
        String refStr = new String(refBases);
        if(!refStr.startsWith(refAllele.toUpperCase())) throw new IllegalStateException("REF mismatch at " + chr + ":" + pos);
        byte[] altHaplotype = (altAllele.toUpperCase() + refStr.substring(refAllele.length())).getBytes();
        Microhomology readHom = Microhomology.findHomology(v, altHaplotype, 0, v.isInsert());
        Microhomology refHom = Microhomology.findHomology(v, refBases, 0, v.isDelete());
        Microhomology hom = readHom;
        if(refHom != null && (hom == null || refHom.Length > hom.Length)) hom = refHom;
        return hom == null ? null : hom.Bases;
    }

    // ---------------------------------------------------------------- normal pileup
    /** 0 = not evaluable, 1 = ref, 2 = alt, 3 = other allele. */
    static int classify(SAMRecord read, int pos, String refAllele, String altAllele, int minBaseQ)
    {
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        int refLen = refAllele.length(), altLen = altAllele.length();
        if(refLen == altLen)
        {
            int first = read.getReadPositionAtReferencePosition(pos);
            if(first == 0) return 0;
            StringBuilder sb = new StringBuilder();
            for(int i = 0; i < refLen; i++)
            {
                int idx = read.getReadPositionAtReferencePosition(pos + i);
                if(idx == 0 || idx != first + i) return 0;
                if(quals[idx - 1] < minBaseQ) return 0;
                sb.append((char)bases[idx - 1]);
            }
            String s = sb.toString().toUpperCase();
            if(s.equals(refAllele.toUpperCase())) return 1;
            if(s.equals(altAllele.toUpperCase())) return 2;
            return 3;
        }
        // indel / complex: bases strictly between the aligned anchors POS-1 and POS+len(REF)
        int left = read.getReadPositionAtReferencePosition(pos - 1);
        int right = read.getReadPositionAtReferencePosition(pos + refLen);
        if(left == 0 || right == 0 || right <= left) return 0;
        String s = new String(bases, left, right - left - 1).toUpperCase();
        if(s.equals(refAllele.toUpperCase())) return 1;
        if(s.equals(altAllele.toUpperCase())) return 2;
        return 3;
    }

    static boolean usable(SAMRecord r, int minMapQ)
    {
        return !r.getReadUnmappedFlag() && !r.isSecondaryOrSupplementary() && !r.getDuplicateReadFlag()
            && !r.getReadFailsVendorQualityCheckFlag() && r.getMappingQuality() >= minMapQ;
    }

    /** returns {DP, refCount, altCount, otherCount}; DP = usable reads whose alignment spans POS. */
    static int[] pileup(Iterable<SAMRecord> reads, int pos, String refAllele, String altAllele, int minMapQ, int minBaseQ)
    {
        int[] out = new int[4];
        for(SAMRecord r : reads)
        {
            if(!usable(r, minMapQ)) continue;
            if(r.getAlignmentStart() > pos || r.getAlignmentEnd() < pos) continue;
            out[0]++;
            int cls = classify(r, pos, refAllele, altAllele, minBaseQ);
            if(cls > 0) out[cls]++;
        }
        return out;
    }

    static int[] pileup(SamReader reader, String chr, int pos, String refAllele, String altAllele, int minMapQ, int minBaseQ)
    {
        List<SAMRecord> reads = new ArrayList<>();
        try(SAMRecordIterator it = reader.query(chr, Math.max(1, pos - 1), pos + refAllele.length(), false))
        {
            while(it.hasNext()) reads.add(it.next());
        }
        return pileup(reads, pos, refAllele, altAllele, minMapQ, minBaseQ);
    }

    // ---------------------------------------------------------------- SAGE reference data
    static ReferenceData loadReference(String[] sageArgs)
    {
        ConfigBuilder config = new ConfigBuilder("Sage");
        SageCallConfig.registerConfig(config);
        config.checkAndParseCommandLine(sageArgs);
        SageCallConfig sage = new SageCallConfig("5.0.2", config);
        ReferenceData reference = new ReferenceData(sage, config);
        if(!reference.load()) throw new IllegalStateException("native SAGE reference loading failed");
        return reference;
    }

    static TierSelector selector(ReferenceData reference, String chr)
    {
        if(!HumanChromosome.contains(chr)) return new TierSelector(List.of(), List.of(), List.of());
        var c = HumanChromosome.fromString(chr);
        return new TierSelector(reference.Hotspots.get(c), reference.PanelWithHotspots.getOrDefault(c, List.of()),
                reference.HighConfidence.getOrDefault(c, List.of()));
    }

    static String contigName(String chromosome, SAMSequenceDictionary d)
    {
        String name = chromosome.startsWith("_") ? chromosome.substring(1) : chromosome;
        for(String cand : new String[]{"chr" + name, name, name.equals("MT") ? "chrM" : "chr" + name, name.equals("MT") ? "M" : name})
            if(d.getSequence(cand) != null) return cand;
        return null;
    }

    static void panelBed(String out, String[] sageArgs) throws Exception
    {
        ReferenceData reference = loadReference(sageArgs);
        try(PrintWriter w = new PrintWriter(new FileWriter(out)))
        {
            int n = 0;
            SAMSequenceDictionary d = reference.RefGenome.getSequenceDictionary();
            for(var e : reference.PanelWithHotspots.entrySet())
            {
                String chr = contigName(e.getKey().toString(), d);
                if(chr == null) { System.err.println("no contig for chromosome " + e.getKey()); continue; }
                for(BaseRegion r : e.getValue()) { w.println(chr + "\t" + (r.start() - 1) + "\t" + r.end()); n++; }
            }
            System.out.println("panel+hotspot regions written: " + n);
        }
        finally { reference.RefGenome.close(); }
    }

    // ---------------------------------------------------------------- annotate
    static void annotate(String[] a) throws Exception
    {
        String candidates = a[0], outVcf = a[1], outTsv = a[2], normalBam = a[3], normalId = a[4], tumorId = a[5], headerExtra = a[6];
        int minMapQ = Integer.parseInt(a[7]), minBaseQ = Integer.parseInt(a[8]);
        String[] sageArgs = Arrays.copyOfRange(a, 9, a.length);
        String fasta = null;
        for(int i = 0; i < sageArgs.length - 1; i++) if(sageArgs[i].equals("-ref_genome")) fasta = sageArgs[i + 1];
        if(fasta == null) throw new IllegalArgumentException("-ref_genome required");

        ReferenceData reference = loadReference(sageArgs);
        IndexedFastaSequenceFile ref = reference.RefGenome;
        SAMSequenceDictionary dict = ref.getSequenceDictionary();
        if(dict == null) throw new IllegalStateException("reference .dict not found");

        Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype (tumour: DeepSomatic call; normal: not called, ./.)"));
        lines.add(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Allelic depths (tumour: DeepSomatic AD; normal: pileup of the normal REDUX BAM)"));
        lines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read depth (tumour: DeepSomatic DP; normal: usable reads spanning POS in the normal REDUX BAM)"));
        lines.add(new VCFFormatHeaderLine("AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allelic frequency (tumour: DeepSomatic VAF; normal: pileup ALT/DP)"));
        lines.add(new VCFInfoHeaderLine("TIER", 1, VCFHeaderLineType.String, "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE] (SAGE 5.0.2 TierSelector with native resources)"));
        lines.add(new VCFInfoHeaderLine("REP_C", 1, VCFHeaderLineType.Integer, "Repeat sequence count (reference-only RepeatInfo.findMaxRepeat as SAGE setMaxRefRepeat)"));
        lines.add(new VCFInfoHeaderLine("REP_S", 1, VCFHeaderLineType.String, "Repeat sequence"));
        lines.add(new VCFInfoHeaderLine("TNC", 1, VCFHeaderLineType.String, "Tri-nucleotide context (reference, centred on POS)"));
        lines.add(new VCFInfoHeaderLine("MH", 1, VCFHeaderLineType.String, "Microhomology (SAGE Microhomology.findHomology on reference/ALT haplotype)"));
        lines.add(new VCFInfoHeaderLine("DS_FILTER", 1, VCFHeaderLineType.String, "Original DeepSomatic FILTER for retained non-PASS records (HOTSPOT/PANEL tier only)"));
        lines.add(new VCFInfoHeaderLine("DS_ORIGINAL_ALLELE", 1, VCFHeaderLineType.String, "Original DeepSomatic record when the representation changed (bcftools norm old-rec-tag): CHR|POS|REF|ALT|USED_ALT_IDX"));
        lines.add(new VCFInfoHeaderLine("DS_NAD", VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "DeepSomatic FORMAT/NAD of the original record projected to this allele (normal counts as computed by DeepSomatic; informational)"));
        lines.add(new VCFInfoHeaderLine("DS_NDP", 1, VCFHeaderLineType.Integer, "DeepSomatic FORMAT/NDP of the original record (informational)"));
        lines.add(new VCFFilterHeaderLine("PASS", "All filters passed"));

        // candidates: tab-separated, header row; sorted in dictionary order
        List<String[]> rows = new ArrayList<>();
        String[] cols = null;
        try(BufferedReader br = Files.newBufferedReader(Paths.get(candidates)))
        {
            for(String line; (line = br.readLine()) != null;)
            {
                if(line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if(cols == null) { cols = f; continue; }
                rows.add(f);
            }
        }
        Map<String, Integer> ci = new HashMap<>();
        for(int i = 0; i < cols.length; i++) ci.put(cols[i], i);
        for(String need : new String[]{"record_id","chrom","pos","ref","alt","qual","ds_filter","t_gt","t_ad_ref","t_ad_alt","t_dp","t_af","orig_allele","ds_nad","ds_ndp"})
            if(!ci.containsKey(need)) throw new IllegalArgumentException("candidates.tsv lacks column " + need);

        Set<String> filterIds = new TreeSet<>();
        for(String[] f : rows) for(String x : f[ci.get("ds_filter")].split(";")) if(!x.equals("PASS") && !x.isEmpty() && !x.equals(".")) filterIds.add(x);
        Map<String, String> dsFilterDesc = Map.of(
                "RefCall", "DeepSomatic: Genotyping model thinks this site is reference.",
                "LowQual", "DeepSomatic: Confidence in this variant being real is below calling threshold.",
                "NoCall", "DeepSomatic: Site has depth=0 resulting in no call.",
                "GERMLINE", "DeepSomatic: Non somatic variants");
        for(String id : filterIds) lines.add(new VCFFilterHeaderLine(id, dsFilterDesc.getOrDefault(id, "DeepSomatic filter " + id)));

        try(BufferedReader br = Files.newBufferedReader(Paths.get(headerExtra)))
        {
            for(String line; (line = br.readLine()) != null;)
            {
                if(!line.startsWith("##") || !line.contains("=")) continue;
                int eq = line.indexOf('=');
                lines.add(new VCFHeaderLine(line.substring(2, eq), line.substring(eq + 1)));
            }
        }
        lines.add(new VCFHeaderLine("sage_compat_adapter_annotator_version", VERSION));
        lines.add(new VCFHeaderLine("sage_compat_adapter_normal_pileup", "bam=" + normalBam + ";min_mapq=" + minMapQ + ";min_baseq_substitutions=" + minBaseQ
                + ";excluded=unmapped,secondary,supplementary,duplicate,qcfail;DP=usable reads spanning POS;indel/complex alleles matched between aligned anchors POS-1 and POS+len(REF) without base-quality filter"));
        lines.add(new VCFHeaderLine("sage_compat_adapter_not_emitted", "LPS,LPS_RC,RC_INFO,RC_CNT,RC_QUAL,RC_JIT,RC_IPC,RC_NM,RC_MH,RC_REPC,RC_REPS,UMI_CNT,AED,AFBQ,AFMQ,ARCBQ,ARMQ,ASBQ,FSB,MUC,RSB,SAC,PON_COUNT,PON_MAX,PON_MULTI,PON_READS,GND_FREQ,MQF,TQP,MED,MSG,NEARBY_INDEL,TINC_RECOVERED (not fabricated; see README)"));

        VCFHeader header = new VCFHeader(lines, Arrays.asList(normalId, tumorId));
        header.setSequenceDictionary(dict);

        Map<String, Integer> contigLen = new HashMap<>();
        for(SAMSequenceRecord s : dict.getSequences()) contigLen.put(s.getSequenceName(), s.getSequenceLength());

        SamReader bam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(fasta)).open(new File(normalBam));

        Map<String, Integer> tierCounts = new TreeMap<>();
        int written = 0, droppedNonPassTier = 0;
        String lastChr = null; int lastPos = -1; TierSelector sel = null; int lastContigIdx = -1;

        try(VariantContextWriter out = new VariantContextWriterBuilder().setOutputFile(outVcf).setReferenceDictionary(dict)
                .setOption(Options.INDEX_ON_THE_FLY).build();
            PrintWriter tsv = new PrintWriter(new FileWriter(outTsv)))
        {
            out.writeHeader(header);
            tsv.println("record_id\tchrom\tpos\tref\talt\tTIER\tREP_C\tREP_S\tTNC\tMH\tn_dp\tn_ad_ref\tn_ad_alt\tn_other\tn_af\tadapter_status");
            for(String[] f : rows)
            {
                String chr = f[ci.get("chrom")]; int pos = Integer.parseInt(f[ci.get("pos")]);
                String refA = f[ci.get("ref")].toUpperCase(), altA = f[ci.get("alt")].toUpperCase();
                if(!contigLen.containsKey(chr)) throw new IllegalArgumentException("contig not in reference dictionary: " + chr);
                int cidx = dict.getSequenceIndex(chr);
                if(!chr.equals(lastChr))
                {
                    if(cidx < lastContigIdx) throw new IllegalStateException("candidates not in dictionary order at " + chr);
                    sel = selector(reference, chr); lastChr = chr; lastPos = -1; lastContigIdx = cidx;
                }
                else if(pos < lastPos) throw new IllegalStateException("candidates not position-sorted at " + chr + ":" + pos);
                lastPos = pos;
                int clen = contigLen.get(chr);

                String observed = new String(upper(ref.getSubsequenceAt(chr, pos, pos + refA.length() - 1).getBases()));
                if(!observed.equals(refA)) throw new IllegalStateException("REF mismatch at " + chr + ":" + pos + " vcf=" + refA + " fasta=" + observed);

                SimpleVariant sv = new SimpleVariant(chr, pos, refA, altA);
                String tier = sel.tier(sv).toString();
                String dsFilter = f[ci.get("ds_filter")];
                boolean isPass = dsFilter.equals("PASS");
                if(!isPass && !(tier.equals("HOTSPOT") || tier.equals("PANEL")))
                {
                    droppedNonPassTier++;
                    tsv.println(String.join("\t", f[ci.get("record_id")], chr, String.valueOf(pos), refA, altA, tier, ".", ".", ".", ".", ".", ".", ".", ".", ".", "dropped_nonpass_tier_" + tier));
                    continue;
                }
                RepeatInfo rep = refRepeat(ref, clen, chr, pos, sv.isIndel());
                String tncS = tnc(ref, clen, chr, pos);
                String mh = microhomology(ref, clen, chr, pos, refA, altA);
                int[] np = pileup(bam, chr, pos, refA, altA, minMapQ, minBaseQ);
                double nAf = np[0] > 0 ? (double)np[2] / np[0] : 0.0;

                Allele refAllele = Allele.create(refA, true), altAllele = Allele.create(altA, false);
                Map<String, Object> info = new LinkedHashMap<>();
                info.put("TIER", tier);
                if(rep != null) { info.put("REP_C", rep.Count); info.put("REP_S", rep.Bases); }
                if(tncS != null) info.put("TNC", tncS);
                if(mh != null) info.put("MH", mh);
                if(!isPass) info.put("DS_FILTER", dsFilter.replace(';', '&'));
                String orig = f[ci.get("orig_allele")];
                if(!orig.isEmpty() && !orig.equals(".")) info.put("DS_ORIGINAL_ALLELE", orig);
                String nad = f[ci.get("ds_nad")], ndp = f[ci.get("ds_ndp")];
                if(!nad.isEmpty() && !nad.equals(".")) info.put("DS_NAD", Arrays.asList(nad.split(",")));
                if(!ndp.isEmpty() && !ndp.equals(".")) info.put("DS_NDP", Integer.parseInt(ndp));

                List<Allele> tGtAlleles = new ArrayList<>();
                for(String g : f[ci.get("t_gt")].split("[/|]")) tGtAlleles.add(g.equals("0") ? refAllele : g.equals("1") ? altAllele : Allele.NO_CALL);
                Genotype tumour = new GenotypeBuilder(tumorId, tGtAlleles)
                        .AD(new int[]{Integer.parseInt(f[ci.get("t_ad_ref")]), Integer.parseInt(f[ci.get("t_ad_alt")])})
                        .DP(Integer.parseInt(f[ci.get("t_dp")])).attribute("AF", f[ci.get("t_af")]).make();
                Genotype normal = new GenotypeBuilder(normalId, Arrays.asList(Allele.NO_CALL, Allele.NO_CALL))
                        .AD(new int[]{np[1], np[2]}).DP(np[0]).attribute("AF", String.format("%.4f", nAf)).make();

                VariantContextBuilder vb = new VariantContextBuilder("deepsomatic", chr, pos, pos + refA.length() - 1, Arrays.asList(refAllele, altAllele))
                        .attributes(info).genotypes(normal, tumour);
                String qual = f[ci.get("qual")];
                if(!qual.equals(".")) vb.log10PError(-Double.parseDouble(qual) / 10.0);
                if(isPass) vb.passFilters(); else for(String x : dsFilter.split(";")) vb.filter(x);
                out.add(vb.make());
                written++;
                tierCounts.merge(tier, 1, Integer::sum);
                tsv.println(String.join("\t", f[ci.get("record_id")], chr, String.valueOf(pos), refA, altA, tier,
                        rep == null ? "0" : String.valueOf(rep.Count), rep == null ? "." : rep.Bases, tncS == null ? "." : tncS, mh == null ? "." : mh,
                        String.valueOf(np[0]), String.valueOf(np[1]), String.valueOf(np[2]), String.valueOf(np[3]), String.format("%.4f", nAf),
                        isPass ? "retained_pass" : "retained_nonpass_" + tier));
            }
        }
        finally { bam.close(); reference.RefGenome.close(); }
        System.out.println("written=" + written + " dropped_nonpass_lowtier=" + droppedNonPassTier + " tiers=" + tierCounts);
    }

    // ---------------------------------------------------------------- self test
    static SAMRecord read(SAMFileHeader h, int start, String cigar, String bases, int mapq)
    {
        SAMRecord r = new SAMRecord(h);
        r.setReadName("r" + start + cigar); r.setReferenceName("chr1"); r.setAlignmentStart(start); r.setCigarString(cigar);
        r.setReadString(bases); StringBuilder q = new StringBuilder(); for(int i = 0; i < bases.length(); i++) q.append('I'); r.setBaseQualityString(q.toString());
        r.setMappingQuality(mapq); r.setReadPairedFlag(false);
        return r;
    }

    static void check(boolean ok, String msg) { if(!ok) throw new AssertionError(msg); }

    static void selfTest()
    {
        // repeats (as Codex ReferenceRepeatAnnotator self-test)
        RepeatInfo a = repeatAt("AAAAAAAAACG".getBytes(), 0); check(a != null && a.Count == 9 && a.Bases.equals("A"), "homopolymer");
        RepeatInfo b = repeatAt("ATATATATATATATATCG".getBytes(), 0); check(b != null && b.Count == 8 && b.Bases.equals("AT"), "dinucleotide");
        check(repeatAt("ACGTGCTACGATC".getBytes(), 0) == null, "non-repeat");
        check(repeatAt("AAACGT".getBytes(), 0) == null, "3 copies is below MIN_REPEAT_COUNT+1");
        // microhomology: deletion of CAG in ...A CAG CAG T...  (REF=ACAG ALT=A at pos0) -> MH CAG, read/ref agree
        SimpleVariant del = new SimpleVariant("chr1", 1, "ACAG", "A");
        Microhomology hDel = Microhomology.findHomology(del, "ACAGCAGT".getBytes(), 0, true); // ref haplotype, applyIndelLength for delete
        check(hDel != null && hDel.Bases.equals("CAG"), "deletion homology on reference: " + hDel);
        Microhomology hDelRead = Microhomology.findHomology(del, "ACAGT".getBytes(), 0, false);
        check(hDelRead != null && hDelRead.Bases.equals("CAG"), "deletion homology on alt haplotype: " + hDelRead);
        SimpleVariant ins = new SimpleVariant("chr1", 1, "A", "AT");
        Microhomology hIns = Microhomology.findHomology(ins, "ATTGG".getBytes(), 0, true); // alt haplotype of ref ATGG
        check(hIns != null && hIns.Bases.equals("T"), "insertion homology on alt haplotype: " + hIns);
        Microhomology hInsRef = Microhomology.findHomology(ins, "ATGG".getBytes(), 0, false);
        check(hInsRef != null && hInsRef.Bases.equals("T"), "insertion homology on reference: " + hInsRef);
        check(Microhomology.findHomology(ins, "ATGGG".getBytes(), 0, true) == null, "no insertion homology");
        // pileup classification
        SAMFileHeader h = new SAMFileHeader(); h.addSequence(new SAMSequenceRecord("chr1", 1000));
        // reference around pos 100: ... we define REF=A ALT=G at 100 for SNV reads
        SAMRecord snvRef = read(h, 91, "20M", "CCCCCCCCCACCCCCCCCCC", 60); // base at 100 = index 9 -> 'A'
        SAMRecord snvAlt = read(h, 91, "20M", "CCCCCCCCCGCCCCCCCCCC", 60);
        SAMRecord lowMq = read(h, 91, "20M", "CCCCCCCCCGCCCCCCCCCC", 5);
        check(classify(snvRef, 100, "A", "G", 20) == 1, "snv ref");
        check(classify(snvAlt, 100, "A", "G", 20) == 2, "snv alt");
        int[] p = pileup(List.of(snvRef, snvAlt, lowMq), 100, "A", "G", 10, 20);
        check(p[0] == 2 && p[1] == 1 && p[2] == 1, "snv pileup counts " + Arrays.toString(p));
        // deletion REF=AC ALT=A at pos 100 (deletes base 101)
        SAMRecord delAlt = read(h, 91, "10M1D10M", "CCCCCCCCCACCCCCCCCCC", 60);
        SAMRecord delRef = read(h, 91, "20M", "CCCCCCCCCACCCCCCCCCC", 60);
        check(classify(delAlt, 100, "AC", "A", 20) == 2, "deletion alt");
        check(classify(delRef, 100, "AC", "A", 20) == 1, "deletion ref");
        // insertion REF=A ALT=AT at pos 100 (T inserted after base 100)
        SAMRecord insAlt = read(h, 91, "10M1I10M", "CCCCCCCCCATCCCCCCCCCC", 60);
        check(classify(insAlt, 100, "A", "AT", 20) == 2, "insertion alt");
        check(classify(delRef, 100, "A", "AT", 20) == 1, "insertion ref");
        check(classify(delAlt, 100, "A", "AT", 20) == 0, "deletion across the insertion anchor is not evaluable (depth only)");
        SAMRecord otherIns = read(h, 91, "10M1I10M", "CCCCCCCCCAGCCCCCCCCCC", 60);
        check(classify(otherIns, 100, "A", "AT", 20) == 3, "other inserted allele at insertion site");
        // read ending before anchor is not evaluable
        SAMRecord shortRead = read(h, 91, "10M", "CCCCCCCCCA", 60);
        check(classify(shortRead, 100, "AC", "A", 20) == 0, "truncated read not evaluable");
        System.out.println("SageCompatAnnotator self-test passed (MAX_REPEAT_LENGTH=" + MAX_REPEAT_LENGTH + ", MIN_REPEAT_COUNT+1=" + (MIN_REPEAT_COUNT + 1)
                + ", DEFAULT_MIN_MAP_QUALITY=" + DEFAULT_MIN_MAP_QUALITY + ", MATCHING_BASE_QUALITY=" + MATCHING_BASE_QUALITY + ")");
    }

    public static void main(String[] args) throws Exception
    {
        if(args.length == 1 && args[0].equals("--self-test")) { selfTest(); return; }
        if(args.length >= 2 && args[0].equals("panel-bed")) { panelBed(args[1], Arrays.copyOfRange(args, 2, args.length)); return; }
        if(args.length >= 10 && args[0].equals("annotate")) { annotate(Arrays.copyOfRange(args, 1, args.length)); return; }
        throw new IllegalArgumentException("usage: --self-test | panel-bed <out.bed> <sage opts> | annotate <candidates.tsv> <out.vcf.gz> <annotations.tsv> <normal.bam> <normal_id> <tumor_id> <header_extra.txt> <min_mapq> <min_baseq> <sage opts>");
    }
}
