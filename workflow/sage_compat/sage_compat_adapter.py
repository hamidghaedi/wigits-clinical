#!/usr/bin/env python3
"""sage-compat adapter: DeepSomatic paired VCF -> SAGE-shaped somatic VCF for oncoanalyser `sage_dir` substitution.

Stages (see README.md and run_sage_compat_adapter.sbatch):
  build     normalised PASS records (+ non-PASS records inside the SAGE panel/hotspot regions) -> candidates.tsv
            with the per-allele ledger; allele selection uses the called tumour GT only.
  finalize  merge the Java annotations into ALLELE_LEDGER.tsv / ledger_summary.json and validate the final VCF
            (header declarations, contig order, sample columns, required INFO/FORMAT, tier parity vs native SAGE).
  purple-summary  append PURPLE outcome counts (PASS, REPORTED, PON-filtered) to ledger_summary.json after pass 2.

Nothing in this module fabricates SAGE read-context fields. Representation changes (multiallelic split,
left-alignment/trimming by bcftools norm) are recorded, never hidden.
"""
import argparse
import collections
import gzip
import json
import sys
from pathlib import Path

VERSION = "0.1.0"
LOW_ALT_SUPPORT = 3  # informational threshold only (SomaReport P0 uses ALT>=3); PASS records are retained regardless


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_dict_order(dict_path):
    order = []
    with open(dict_path) as f:
        for line in f:
            if line.startswith("@SQ"):
                fields = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:] if ":" in x)
                order.append(fields["SN"])
    if not order:
        raise ValueError("no @SQ lines in " + str(dict_path))
    return {name: i for i, name in enumerate(order)}


def vcf_records(path):
    samples = None
    with opener(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                if samples is None:
                    raise ValueError("records before #CHROM in " + str(path))
                yield samples, line.rstrip("\n").split("\t")


def info_dict(text):
    out = {}
    for item in text.split(";"):
        if not item or item == ".":
            continue
        k, _, v = item.partition("=")
        out[k] = v if _ else True
    return out


def gt_calls_alt(gt):
    return any(a == "1" for a in gt.replace("|", "/").split("/"))


def classify_representation(chrom, pos, ref, alt, orig_tag):
    """DS_ORIGINAL_ALLELE = CHR|POS|REF|ALT|USED_ALT_IDX from bcftools norm --old-rec-tag."""
    if not orig_tag or orig_tag == ".":
        return "unchanged", 1
    parts = orig_tag.split("|")
    o_chrom, o_pos, o_ref, o_alt = parts[0], parts[1], parts[2], parts[3]
    used_idx = int(parts[4]) if len(parts) > 4 and parts[4] else 1
    changes = []
    if "," in o_alt:
        changes.append("multiallelic_split")
    if o_pos != pos or o_ref != ref or o_alt.split(",")[used_idx - 1] != alt:
        changes.append("realigned_or_trimmed")
    return "+".join(changes) if changes else "unchanged", used_idx


def project_row(samples, row, tumor_id, source_class):
    """Return (candidate dict or None, ledger dict) for a biallelic normalised DeepSomatic record."""
    chrom, pos, _, ref, alt, qual, filt, info, fmt = row[:9]
    values = dict(zip(fmt.split(":"), row[9 + samples.index(tumor_id)].split(":")))
    inf = info_dict(info)
    orig = inf.get("DS_ORIGINAL_ALLELE", ".")
    representation, used_idx = classify_representation(chrom, pos, ref, alt, orig if orig is not True else ".")
    record_id = (orig if orig not in (".", True) else f"{chrom}|{pos}|{ref}|{alt}|1")
    ledger = dict(record_id=record_id, chrom=chrom, pos=pos, ref=ref, alt=alt, ds_filter=filt, ds_qual=qual,
                  t_gt=values.get("GT", "./."), t_ad=values.get("AD", "."), t_dp=values.get("DP", "."),
                  t_vaf=values.get("VAF", "."), ds_nad=values.get("NAD", "."), ds_ndp=values.get("NDP", "."),
                  representation=representation, source_class=source_class)
    if alt.startswith("<") or alt in ("*", ".", "N"):
        ledger["disposition"] = "dropped_symbolic_or_missing_alt"
        return None, ledger
    ad = ledger["t_ad"].split(",")
    if len(ad) != 2 or "." in ad:
        ledger["disposition"] = "dropped_missing_tumour_ad"
        return None, ledger
    ledger["low_alt_support"] = int(ad[1]) < LOW_ALT_SUPPORT
    if filt == "PASS":
        # PASS: the allele must be in the called tumour GT (uncalled ALTs of multiallelic sites are dropped)
        if not gt_calls_alt(ledger["t_gt"]):
            ledger["disposition"] = "dropped_alt_not_in_called_gt"
            return None, ledger
        ledger["disposition"] = "candidate_pass"
    else:
        # non-PASS (panel regions only): DeepSomatic writes GERMLINE/LowQual records with the somatic GT 0/0, so the
        # called-GT rule cannot apply; RefCall/NoCall are not variant calls and are dropped, the rest need ALT support.
        filters = set(filt.split(";"))
        if filters & {"RefCall", "NoCall"}:
            ledger["disposition"] = "dropped_nonpass_no_variant_call"
            return None, ledger
        if int(ad[1]) == 0:
            ledger["disposition"] = "dropped_nonpass_no_alt_support"
            return None, ledger
        ledger["disposition"] = "candidate_nonpass"
    cand = dict(record_id=record_id, chrom=chrom, pos=int(pos), ref=ref, alt=alt, qual=qual, ds_filter=filt,
                t_gt=ledger["t_gt"], t_ad_ref=ad[0], t_ad_alt=ad[1], t_dp=ledger["t_dp"],
                t_af=ledger["t_vaf"] if ledger["t_vaf"] != "." else (f"{int(ad[1])/int(ledger['t_dp']):.4f}" if ledger["t_dp"] not in (".", "0") else "0"),
                orig_allele=orig if orig not in (".", True) else ".", ds_nad=ledger["ds_nad"], ds_ndp=ledger["ds_ndp"])
    return cand, ledger


CAND_COLS = ["record_id", "chrom", "pos", "ref", "alt", "qual", "ds_filter", "t_gt", "t_ad_ref", "t_ad_alt", "t_dp", "t_af",
             "orig_allele", "ds_nad", "ds_ndp"]
LEDGER_COLS = ["record_id", "chrom", "pos", "ref", "alt", "source_class", "ds_filter", "ds_qual", "t_gt", "t_ad", "t_dp", "t_vaf",
               "ds_nad", "ds_ndp", "representation", "low_alt_support", "disposition"]


def build(args):
    order = read_dict_order(args.dict)
    candidates, ledger, counts = [], [], collections.Counter()
    seen = set()
    for path, source_class in ((args.pass_vcf, "deepsomatic_pass"), (args.panel_vcf, "deepsomatic_nonpass_in_panel")):
        if not path:
            continue
        for samples, row in vcf_records(path):
            if source_class == "deepsomatic_nonpass_in_panel" and row[6] == "PASS":
                counts["panel_extract_pass_skipped"] += 1  # already handled through the PASS stream
                continue
            if "," in row[4]:
                raise ValueError(f"multiallelic record survived normalisation at {row[0]}:{row[1]}")
            cand, led = project_row(samples, row, args.tumor_id, source_class)
            key = (row[0], row[1], row[3], row[4])
            if key in seen:
                led["disposition"] = "dropped_duplicate_allele"
                cand = None
            seen.add(key)
            counts[led["disposition"]] += 1
            counts["representation:" + led["representation"]] += 1
            if cand:
                candidates.append(cand)
                if led.get("low_alt_support"):
                    counts["candidate_low_alt_support_lt%d" % LOW_ALT_SUPPORT] += 1
            ledger.append(led)
    missing = {c["chrom"] for c in candidates} - set(order)
    if missing:
        raise ValueError("contigs absent from reference dictionary: " + ",".join(sorted(missing)))
    candidates.sort(key=lambda c: (order[c["chrom"]], c["pos"], c["ref"], c["alt"]))
    with open(args.out_candidates, "w") as f:
        f.write("\t".join(CAND_COLS) + "\n")
        for c in candidates:
            f.write("\t".join(str(c[k]) for k in CAND_COLS) + "\n")
    with open(args.out_ledger, "w") as f:
        f.write("\t".join(LEDGER_COLS) + "\n")
        for l in ledger:
            f.write("\t".join(str(l.get(k, ".")) for k in LEDGER_COLS) + "\n")
    summary = dict(version=VERSION, candidates=len(candidates), ledger_rows=len(ledger), counts=dict(counts))
    Path(args.out_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def parse_filter_counts(path):
    counts = {}
    if not path:
        return counts
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                counts[parts[1]] = int(parts[0])
    return counts


def validate_vcf(path, dict_path, normal_id, tumor_id):
    order = read_dict_order(dict_path)
    declared = dict(INFO=set(), FORMAT=set(), FILTER=set())
    contigs_in_header = []
    samples = None
    used = dict(INFO=set(), FORMAT=set(), FILTER=set())
    problems, n = [], 0
    last = None
    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                for kind in declared:
                    if line.startswith(f"##{kind}=<ID="):
                        declared[kind].add(line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">\n"))
                if line.startswith("##contig=<ID="):
                    contigs_in_header.append(line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">\n"))
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            row = line.rstrip("\n").split("\t")
            n += 1
            key = (order.get(row[0], -1), int(row[1]))
            if key[0] < 0:
                problems.append(f"contig not in dict: {row[0]}")
            if last and key < last:
                problems.append(f"out of order at {row[0]}:{row[1]}")
            last = key
            inf = info_dict(row[7])
            used["INFO"].update(inf)
            used["FORMAT"].update(row[8].split(":"))
            used["FILTER"].update(x for x in row[6].split(";"))
            if "TIER" not in inf:
                problems.append(f"missing TIER at {row[0]}:{row[1]}")
            if "TNC" not in inf:
                problems.append(f"missing TNC at {row[0]}:{row[1]}")
            keys = row[8].split(":")
            for col in row[9:]:
                g = dict(zip(keys, col.split(":")))
                if "." in g.get("AD", ".").split(",") or g.get("DP", ".") == ".":
                    problems.append(f"missing AD/DP at {row[0]}:{row[1]}")
            if len(problems) > 50:
                break
    if samples != [normal_id, tumor_id]:
        problems.append(f"sample columns {samples} != [{normal_id}, {tumor_id}]")
    for kind in declared:
        undeclared = used[kind] - declared[kind] - {"PASS"}
        if undeclared:
            problems.append(f"undeclared {kind}: {sorted(undeclared)}")
    if [c for c in contigs_in_header] != sorted(contigs_in_header, key=lambda c: order.get(c, 1e9)):
        problems.append("header contig order differs from dictionary")
    return dict(records=n, problems=problems[:50], declared={k: sorted(v) for k, v in declared.items()},
                used={k: sorted(v) for k, v in used.items()}, samples=samples)


def tier_parity(final_vcf, sage_vcf):
    """TIER agreement on the exact-allele intersection with the native SAGE somatic VCF (all SAGE records)."""
    ours = {}
    for _, row in vcf_records(final_vcf):
        ours[(row[0], row[1], row[3], row[4])] = info_dict(row[7]).get("TIER")
    agree = disagree = 0
    sage_pass_total = sage_pass_shared = 0
    changes = collections.Counter()
    for _, row in vcf_records(sage_vcf):
        key = (row[0], row[1], row[3], row[4])
        if row[6] == "PASS":
            sage_pass_total += 1
            if key in ours:
                sage_pass_shared += 1
        if key in ours:
            st = info_dict(row[7]).get("TIER")
            if st == ours[key]:
                agree += 1
            else:
                disagree += 1
                changes[f"{st}->{ours[key]}"] += 1
    return dict(intersection=agree + disagree, tier_agree=agree, tier_disagree=disagree, tier_changes=dict(changes),
                sage_pass_total=sage_pass_total, sage_pass_in_adapter=sage_pass_shared, adapter_records=len(ours))


def finalize(args):
    cand_ledger = []
    with open(args.candidates_ledger) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            cand_ledger.append(dict(zip(cols, line.rstrip("\n").split("\t"))))
    ann = {}
    with open(args.annotations) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(cols, line.rstrip("\n").split("\t")))
            ann[(d["chrom"], d["pos"], d["ref"], d["alt"])] = d
    counts, flags = collections.Counter(), collections.Counter()
    out_cols = LEDGER_COLS + ["TIER", "REP_C", "REP_S", "TNC", "MH", "n_dp", "n_ad_ref", "n_ad_alt", "n_af", "adapter_status", "final_disposition"]
    with open(args.out_ledger, "w") as f:
        f.write("\t".join(out_cols) + "\n")
        for l in cand_ledger:
            a = ann.get((l["chrom"], l["pos"], l["ref"], l["alt"]), {})
            if l["disposition"].startswith("candidate"):
                status = a.get("adapter_status", "missing_from_annotator_output")
            else:
                status = l["disposition"]
            l.update({k: a.get(k, ".") for k in ["TIER", "REP_C", "REP_S", "TNC", "MH", "n_dp", "n_ad_ref", "n_ad_alt", "n_af", "adapter_status"]})
            l["final_disposition"] = status
            counts[status] += 1
            if status == "retained_pass" and l.get("low_alt_support") == "True":
                flags["retained_pass_low_alt_support_lt%d" % LOW_ALT_SUPPORT] += 1
            if status.startswith("retained") and a.get("n_ad_alt", "0") not in ("0", "."):
                flags["retained_with_normal_pileup_alt_support"] += 1
            if status.startswith("retained") and l.get("representation", "unchanged") != "unchanged":
                flags["retained_representation_" + l["representation"]] += 1
            f.write("\t".join(str(l.get(k, ".")) for k in out_cols) + "\n")
    filter_counts = parse_filter_counts(args.filter_counts)
    validation = validate_vcf(args.final_vcf, args.dict, args.normal_id, args.tumor_id)
    parity = tier_parity(args.final_vcf, args.sage_vcf) if args.sage_vcf else None
    tiers = collections.Counter(v["TIER"] for v in ann.values() if v["adapter_status"].startswith("retained"))
    summary = dict(version=VERSION,
                   deepsomatic_records_by_filter=filter_counts,
                   deepsomatic_records_total=sum(filter_counts.values()) if filter_counts else None,
                   ledger_rows=len(cand_ledger), dispositions=dict(counts), flags=dict(flags),
                   retained_total=sum(v for k, v in counts.items() if k.startswith("retained")),
                   retained_tiers=dict(tiers),
                   bulk_classes_note="RefCall/NoCall records and GERMLINE records outside the SAGE panel/hotspot regions are "
                                     "represented by the FILTER counts only (dropped_nonpass_filter); every PASS record and "
                                     "every non-PASS record inside the panel regions has its own ledger row.",
                   final_vcf=str(args.final_vcf), validation=validation, tier_parity_vs_sage=parity)
    Path(args.out_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("dispositions", "flags", "retained_total", "retained_tiers", "tier_parity_vs_sage")}, indent=2))
    print("validation problems:", validation["problems"])
    if validation["problems"]:
        sys.exit(2)


def purple_summary(args):
    summary = json.loads(Path(args.summary_json).read_text())
    counts = collections.Counter()
    adapter_alleles = set()
    with open(args.ledger) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(cols, line.rstrip("\n").split("\t")))
            if d["final_disposition"].startswith("retained"):
                adapter_alleles.add((d["chrom"], d["pos"], d["ref"], d["alt"]))
    in_purple = set()
    for _, row in vcf_records(args.purple_vcf):
        key = (row[0], row[1], row[3], row[4])
        in_purple.add(key)
        inf = info_dict(row[7])
        filters = set(row[6].split(";"))
        counts["purple_records"] += 1
        if row[6] == "PASS":
            counts["purple_pass"] += 1
            if "REPORTED" in inf:
                counts["purple_pass_reported"] += 1
            if "HOTSPOT" in inf:
                counts["purple_pass_hotspot"] += 1
        for fl in filters:
            if fl != "PASS":
                counts["purple_filter:" + fl] += 1
        if "REPORTED" in inf:
            counts["purple_reported_any_filter"] += 1
        if key not in adapter_alleles:
            counts["purple_records_not_in_adapter_retained"] += 1
    counts["adapter_retained_missing_from_purple"] = len(adapter_alleles - in_purple)
    drivers = []
    if args.driver_catalog:
        with open(args.driver_catalog) as f:
            cols = f.readline().rstrip("\n").split("\t")
            for line in f:
                d = dict(zip(cols, line.rstrip("\n").split("\t")))
                drivers.append({k: d.get(k) for k in ("gene", "driver", "category", "likelihoodMethod", "driverLikelihood", "reportedStatus") if k in d})
    summary["purple"] = dict(purple_vcf=str(args.purple_vcf), counts=dict(counts), driver_catalog_rows=len(drivers), drivers=drivers)
    Path(args.summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["purple"]["counts"], indent=2))


def merge_bed(args):
    """Sort and merge BED intervals (dictionary contig order) so bcftools -R never emits a record twice."""
    order = read_dict_order(args.dict)
    by_chrom = collections.defaultdict(list)
    with open(args.bed) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            c, s, e = line.rstrip("\n").split("\t")[:3]
            by_chrom[c].append((int(s), int(e)))
    n = 0
    with open(args.out, "w") as out:
        for c in sorted(by_chrom, key=lambda x: order.get(x, 1e9)):
            merged = []
            for s, e in sorted(by_chrom[c]):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            for s, e in merged:
                out.write(f"{c}\t{s}\t{e}\n")
                n += 1
    print(f"merged intervals: {n}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--pass-vcf", required=True)
    b.add_argument("--panel-vcf")
    b.add_argument("--dict", required=True)
    b.add_argument("--tumor-id", required=True)
    b.add_argument("--out-candidates", required=True)
    b.add_argument("--out-ledger", required=True)
    b.add_argument("--out-json", required=True)
    b.set_defaults(func=build)
    fz = sub.add_parser("finalize")
    fz.add_argument("--candidates-ledger", required=True)
    fz.add_argument("--annotations", required=True)
    fz.add_argument("--filter-counts")
    fz.add_argument("--final-vcf", required=True)
    fz.add_argument("--dict", required=True)
    fz.add_argument("--normal-id", required=True)
    fz.add_argument("--tumor-id", required=True)
    fz.add_argument("--sage-vcf")
    fz.add_argument("--out-ledger", required=True)
    fz.add_argument("--out-json", required=True)
    fz.set_defaults(func=finalize)
    ps = sub.add_parser("purple-summary")
    ps.add_argument("--ledger", required=True)
    ps.add_argument("--summary-json", required=True)
    ps.add_argument("--purple-vcf", required=True)
    ps.add_argument("--driver-catalog")
    ps.set_defaults(func=purple_summary)
    mb = sub.add_parser("merge-bed")
    mb.add_argument("--bed", required=True)
    mb.add_argument("--dict", required=True)
    mb.add_argument("--out", required=True)
    mb.set_defaults(func=merge_bed)
    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
