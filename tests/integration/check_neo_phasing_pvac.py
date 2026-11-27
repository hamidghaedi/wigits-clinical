#!/usr/bin/env python3
"""Run under Slurm inside the pinned pVACtools 5.4 image, with --outdir.

Synthetic positive/negative controls exercise pVACseq's real parser, converter
and protein generator. No models, network, reference data or real reads needed.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow'))
import neo_phasing


def check(outdir):
    import pysam
    import vcfpy
    from pvactools.lib.proximal_variant import ProximalVariant
    from pvactools.lib.input_file_converter import VcfConverter
    from pvactools.tools.pvacseq import generate_protein_fasta
    from Bio import SeqIO

    out = Path(outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    csq = 'Allele|Consequence|Feature|SYMBOL|Gene|Protein_position|Amino_acids|Codons|TSL|BIOTYPE|WildtypeProtein|FrameshiftSequence'
    header = ('##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n'
              '##INFO=<ID=LPS,Number=.,Type=Integer,Description="Local phase set">\n'
              '##INFO=<ID=LPS_RC,Number=.,Type=Integer,Description="Phase support">\n'
              f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: {csq}">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOUR\n')

    def row(pos, protein, alt='T', aa='S', gt='0/1', lps='7', support='4', transcript='ENST_SYNTHETIC'):
        pos += 1000
        codon = 'Gcc/' + ('Acc' if alt == 'A' else 'Tcc')
        info = f'CSQ={alt}|missense_variant|{transcript}|SYNTHETIC|ENSG_SYNTHETIC|{protein}|A/{aa}|{codon}|1|protein_coding|' + 'A'*30 + '|'
        if lps:
            info += f';LPS={lps};LPS_RC={support}'
        return f'chr1\t{pos}\t.\tG\t{alt}\t50\tPASS\t{info}\tGT:AD\t{gt}:20,10\n'

    anchor = row(100,10,'A','T')
    cis = row(103,11)
    unknown = row(106,12,'A','T',lps='')
    homozygous = row(109,13,gt='1/1',lps='')
    other_group = row(112,14,lps='8')
    other_transcript = row(115,15,transcript='ENST_OTHER')
    source = out/'source.vcf'
    source.write_text(header+(anchor+cis+unknown+homozygous+other_group+other_transcript).replace('0/1:20,10','./.:20,10'))
    selected = out/'selected.vcf';selected.write_text(header+anchor+cis)
    adapted = out/'context.vcf'
    somatic = out/'somatic.vcf';somatic.write_text(header+anchor+cis)
    germline = out/'germline.vcf';germline.write_text(header+unknown+homozygous+other_group+other_transcript)
    neo_phasing.export_context(str(somatic),str(germline),str(source),'TUMOUR',
        {('chr1','1100','G','A'),('chr1','1103','G','T')},adapted,out/'evidence.json')

    def index(path):
        target = str(path)+'.gz'
        pysam.tabix_compress(str(path),target,force=True)
        pysam.tabix_index(target,preset='vcf',force=True)
        return target

    selected_gz, context_gz = index(selected), index(adapted)
    main = next(vcfpy.Reader.from_path(selected_gz))
    parser = ProximalVariant(context_gz,True,100)
    positions = sorted(r.POS for r, csq, alt in parser.extract(main,'A','ENST_SYNTHETIC'))
    assert positions == [1103,1109], positions
    converted, proximal = out/'converted.tsv', out/'proximal.tsv'
    VcfConverter(input_file=selected_gz,output_file=str(converted),sample_name='TUMOUR',
        proximal_variants_vcf=context_gz,proximal_variants_tsv=str(proximal),
        flanking_bases=100,pass_only=True).execute()
    rows = list(csv.DictReader(proximal.open(),delimiter='\t'))
    anchor_rows = [r for r in rows if r['main_somatic_variant'].startswith('1.')]
    assert sorted(int(r['stop']) for r in anchor_rows) == [1103,1109], rows
    assert {r['type'] for r in anchor_rows if r['stop']=='1103'} == {'somatic'}
    assert {r['type'] for r in anchor_rows if r['stop']=='1109'} == {'germline'}
    fasta = out/'peptides.fasta'
    generate_protein_fasta.main([selected_gz,'10',str(fasta),'-s','TUMOUR','-p',context_gz,'--pass-only'])
    peptides = list(SeqIO.parse(str(fasta),'fasta'))
    mutant = [str(p.seq) for p in peptides if p.id.startswith('MT.1.')]
    wildtype = [str(p.seq) for p in peptides if p.id.startswith('WT.1.')]
    # Mutation at 10 plus cis at 11 and actual homozygote at 13. Unknown at
    # 12, different group at 14 and different transcript at 15 are not applied.
    assert mutant == ['A'*9+'TSAS'+'A'*8], mutant
    assert wildtype == ['A'*12+'S'+'A'*8], wildtype

    # Missing HP on two unknown records must not create a cis relationship.
    unresolved = out/'unresolved_source.vcf'
    unresolved.write_text(header+row(100,10,'A','T',lps='')+unknown)
    unresolved_out = out/'unresolved.vcf'
    neo_phasing.export_context(str(unresolved),None,str(unresolved),'TUMOUR',
        {('chr1','1100','G','A')},unresolved_out,out/'unresolved_evidence.json')
    assert ProximalVariant(index(unresolved_out),True,100).extract(main,'A','ENST_SYNTHETIC') == []
    result = dict(status='PASS',tool='pVACtools 5.4.0',accepted_positions=positions,
                  converter_rows=len(rows),mutant=mutant,wildtype=wildtype,
                  unknown_phase_not_linked=True,other_group_not_linked=True,other_transcript_not_applied=True)
    (out/'completion.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',required=True)
    check(ap.parse_args().outdir)
