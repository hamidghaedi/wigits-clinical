#!/usr/bin/env python3
"""Inventory neoantigen-relevant oncoanalyser outputs for a set of runs. Read-only. Prints a markdown table + JSON."""
import csv, glob, gzip, json, os, sys
RUNS = {
 "KQ0004": "/global/scratch/hpc6049/wigits_runs/KQ0004_full_native_20261114_onco3.0.0",
 "KQ0005": "/global/scratch/hpc6049/wigits_runs/KQ0005_full_native_20261116_onco3.0.0",
 "KQ0012": "/global/scratch/hpc6049/wigits_runs/KQ0012_full_native_20261116_onco3.0.0",
 "KQ0013": "/global/scratch/hpc6049/wigits_runs/KQ0013_full_native_20261116_onco3.0.0",
 "KQ0017": "/global/scratch/hpc6049/wigits_runs/KQ0017_full_native_20261116_onco3.0.0",
 "KQ0205": "/global/scratch/hpc6049/wigits_runs/KQ0205_full_native_20261109_onco3.0.0",
 "subject_a": "/global/scratch/hpc6049/oncoanalyser/test_3.0.0",
}
def first(pat):
    g = sorted(glob.glob(pat)); return g[0] if g else None
def nrows(path, skip_header=True):
    if not path: return None
    op = gzip.open if path.endswith(".gz") else open
    n = 0
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            n += 1
    return n - (1 if skip_header and not path.endswith(".vcf.gz") else 0)
def tsv(path):
    with open(path) as fh: return list(csv.DictReader(fh, delimiter="\t"))
out = {}
for case, root in RUNS.items():
    res = first(f"{root}/results/*/neo") and os.path.dirname(first(f"{root}/results/*/neo"))
    d = {"root": root, "results_case_dir": res}
    if not res:
        d["status"] = "no results/<case>/neo"; out[case] = d; continue
    ps = first(f"{res}/neo/scorer/*.neo.peptide_scores.tsv"); nd = first(f"{res}/neo/finder/*.neo_data.tsv")
    d["neo_peptide_rows"] = nrows(ps); d["neo_candidates"] = nrows(nd)
    if nd:
        vt = {}
        for r in tsv(nd): vt[r["VariantType"]] = vt.get(r["VariantType"], 0) + 1
        d["neo_variant_types"] = vt
    if ps:
        rows = tsv(ps); d["neo_alleles"] = sorted({r["Allele"] for r in rows})
        d["neo_allele_disrupted_rows"] = sum(1 for r in rows if r.get("AlleleDisrupted","false").lower()=="true")
    li = first(f"{res}/lilac/*.lilac.tsv"); lq = first(f"{res}/lilac/*.lilac.qc.tsv")
    if li:
        rows = tsv(li); d["lilac_alleles"] = [r["Allele"] for r in rows]
        d["lilac_tumor_cn"] = [r.get("TumorCopyNumber") for r in rows]
        d["lilac_rna_total"] = [r.get("RnaTotal") for r in rows]
    if lq:
        q = tsv(lq); d["lilac_qc"] = {k: q[0][k] for k in q[0] if k in ("Status","ScoreMargin","TotalFragments","FittedFragments","DiscardedAlignmentFragments","HlaYAllele","Status")}
    sa = first(f"{res}/sage_append/*tumor*/*.sage.append.vcf.gz") or first(f"{res}/sage_append/*/*.sage.append.vcf.gz")
    d["sage_append_tumor_vcf"] = os.path.basename(sa) if sa else None
    if sa:
        n=0; rna_cov=0; rna_alt=0
        with gzip.open(sa,"rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    if line.startswith("#CHROM"): samples=line.rstrip().split("\t")[9:]
                    continue
                n+=1; f=line.rstrip().split("\t"); fmt=f[8].split(":")
                # RNA sample is the appended (last) column
                vals=dict(zip(fmt,f[-1].split(":")))
                try:
                    dp=int(vals.get("DP","0")); ad=vals.get("AD","0,0").split(","); alt=int(ad[1]) if len(ad)>1 else 0
                except ValueError: dp=0; alt=0
                rna_cov += dp>=10; rna_alt += alt>=2
        d["sage_append_records"]=n; d["sage_append_samples"]=samples; d["rna_dp_ge10"]=rna_cov; d["rna_alt_ge2"]=rna_alt
    gd = first(f"{res}/isofox/*.isf.gene_data.tsv"); d["isofox_genes"] = nrows(gd)
    pf = first(f"{res}/isofox/*.isf.pass_fusions.tsv"); d["isofox_pass_fusions"] = nrows(pf)
    ln = first(f"{res}/linx/somatic_annotations/*.linx.neoepitope.tsv"); d["linx_neoepitopes"] = nrows(ln)
    pp = first(f"{res}/purple/*.purple.purity.tsv")
    if pp:
        r = tsv(pp)[0]; d["purple"] = {k: r.get(k) for k in ("purity","ploidy","status","qcStatus","msStatus","tmbPerMb","tml","wholeGenomeDuplication")}
    sv = first(f"{res}/purple/*.purple.somatic.vcf.gz"); d["purple_somatic_records"] = nrows(sv)
    oj = first(f"{res}/orange/*.orange.json"); d["orange_json"] = bool(oj)
    out[case] = d
json.dump(out, open(sys.argv[1] if len(sys.argv)>1 else "case_inventory.json","w"), indent=1)
hdr = ["case","purity","ploidy","PURPLE qc","somatic recs","NEO cands","NEO peptide rows","NEO types","LILAC alleles","LILAC qc","RNA recs (DP>=10 / alt>=2)","ISOFOX genes","pass fusions","LINX neoep","ORANGE"]
print("| " + " | ".join(hdr) + " |"); print("|" + "---|"*len(hdr))
for case,d in out.items():
    if "status" in d: print(f"| {case} | {d['status']} |" + " |"*(len(hdr)-2)); continue
    p=d.get("purple",{}) or {}
    la=",".join(d.get("lilac_alleles",[])[:6]); 
    print(f"| {case} | {p.get('purity')} | {p.get('ploidy')} | {p.get('qcStatus')} | {d.get('purple_somatic_records')} | {d.get('neo_candidates')} | {d.get('neo_peptide_rows')} | {d.get('neo_variant_types')} | {la} | {(d.get('lilac_qc') or {}).get('Status')} | {d.get('sage_append_records')} ({d.get('rna_dp_ge10')} / {d.get('rna_alt_ge2')}) | {d.get('isofox_genes')} | {d.get('isofox_pass_fusions')} | {d.get('linx_neoepitopes')} | {d.get('orange_json')} |")
