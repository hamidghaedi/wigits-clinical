#!/usr/bin/env python3
"""Plan and record independently selected post-native blood CH and RNA neo analysis.

Biological eligibility comes from the bound laboratory metadata. Native outputs
are read-only prerequisites. No specimen type is inferred from a sample suffix.
The caller, curation and neo ranking remain the feature branch implementations.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import sys
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parent))
import neo_interface as ni

SCHEMA = "wigits.optional-extensions.v1"


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n")
    temp.replace(path)


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def resource_allowance(stage, workload):
    """Conservative initial bounds; no full BAM scans or CPU inflation.

    CH uses paired BAM bytes (GiB); neo annotation/prediction uses selected
    variants x class-I alleles, measured by the completed export task. These
    initial margins are explicit allowances, not benchmarked throughput claims.
    """
    if stage == "blood_ch":
        return dict(cpus=8, memory_gib=48, minutes=120 + 30 * math.ceil(workload / 100))
    if stage == "neo_prepare":
        return dict(cpus=1, memory_gib=8, minutes=30 + 30 * math.ceil(workload / 0.25))
    if stage == "neo_predict":
        return dict(cpus=8, memory_gib=32, minutes=120 + 60 * math.ceil(workload / 500))
    raise ValueError("Unknown stage")


def require_identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError("Invalid sample/case identifier")


def disposition(enabled, eligible, reason, missing):
    if not enabled:
        return dict(status="DISABLED", reason="Not selected")
    if not eligible:
        return dict(status="INELIGIBLE", reason=reason)
    if missing:
        return dict(status="MISSING_INPUTS", reason="Missing prerequisites", missing=missing)
    return dict(status="READY", reason="Selected with laboratory eligibility and native prerequisites")


def build_plan(config):
    if config.get("schema_version") != SCHEMA or not config.get("cases"):
        raise ValueError("Expected optional extension v1 configuration and cases")
    resources = config.get("resources", {})
    records, seen = [], set()
    common = ("reference_fasta", "samtools_image")
    for original in config["cases"]:
        case = dict(original)
        for key in ("case_id", "tumour_sample", "normal_sample"):
            require_identifier(case[key])
        if case.get("rna_sample"):
            require_identifier(case["rna_sample"])
        if case["case_id"] in seen:
            raise ValueError("Duplicate case")
        seen.add(case["case_id"])
        if case.get("native_genome") != "GRCh38_hmf":
            raise ValueError("These extensions currently require the native GRCh38_hmf route")
        selected = case.get("options", {})
        if any(type(selected.get(key, False)) is not bool for key in ("blood_ch", "rna_neoantigens")):
            raise ValueError("Options must be booleans")
        lab = case.get("laboratory", {})
        root = Path(case["results_case_dir"])
        if not root.is_absolute():
            raise ValueError("Native case results path must be absolute")
        normal_bam = root / "alignments" / case["normal_sample"] / (case["normal_sample"] + ".redux.bam")
        tumour_bam = root / "alignments" / case["tumour_sample"] / (case["tumour_sample"] + ".redux.bam")
        purity = root / "purple" / (case["tumour_sample"] + ".purple.purity.tsv")
        is_blood = str(lab.get("normal_specimen_type", "")).strip().casefold() == "blood" and bool(lab.get("source"))
        rna_available = bool(case.get("rna_sample")) and lab.get("tumour_rna_available") is True and bool(lab.get("source"))

        def missing_resources(keys):
            return [key for key in keys if not resources.get(key) or not Path(resources[key]).exists()]

        ch_inputs = [normal_bam, Path(str(normal_bam) + ".bai"), tumour_bam, Path(str(tumour_bam) + ".bai"), purity]
        ch_missing = [str(p) for p in ch_inputs if not p.is_file()]
        ch_missing += missing_resources(common + ("sage_image", "pave_image", "hmf_resources"))
        if not resources.get("bind_roots"):
            ch_missing.append("bind_roots")
        case["blood_ch"] = disposition(selected.get("blood_ch", False), is_blood,
                                        "Normal is not confirmed blood in the bound laboratory record", ch_missing)
        case["blood_ch"]["resources"] = resource_allowance("blood_ch", sum(p.stat().st_size for p in (normal_bam, tumour_bam) if p.is_file()) / 2**30)

        neo_missing = missing_resources(common + ("vep_image", "vep_cache", "vep_plugins", "pvactools_image", "genes_of_interest"))
        if not resources.get("bind_roots"):
            neo_missing.append("bind_roots")
        try:
            layout = ni.CaseLayout.from_results_dir(str(root))
            if (layout.tumour_dna, layout.normal_dna) != (case["tumour_sample"], case["normal_sample"]):
                raise ValueError("Native VCF sample identities disagree with the configured case")
            if rna_available and layout.tumour_rna != case["rna_sample"]:
                neo_missing.append("Recorded RNA sample is absent from the native SAGE_APPEND header")
            for key in ("sage_append_vcf", "isofox_gene_tsv", "isofox_transcript_tsv", "lilac_tsv", "lilac_qc_tsv"):
                if not layout.files.get(key):
                    neo_missing.append(key)
            case["neo_somatic_bytes"] = Path(layout.files["purple_somatic_vcf"]).stat().st_size
        except FileNotFoundError as exc:
            neo_missing.append(str(exc))
            case["neo_somatic_bytes"] = 0
        case["rna_neoantigens"] = disposition(selected.get("rna_neoantigens", False), rna_available,
                                               "No tumour RNA confirmed in the bound laboratory record", neo_missing)
        case["rna_neoantigens"]["resources"] = resource_allowance("neo_prepare", case["neo_somatic_bytes"] / 2**30)
        records.append(case)
    return dict(schema_version=SCHEMA, generated_at=timestamp(), resources=resources, cases=records,
                native_results_modified=False, report_template_modified=False)


def collect_status(plan, result_root, trace=None):
    result_root = Path(result_root)
    terminal = {}
    if trace and Path(trace).is_file():
        with Path(trace).open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                process = row.get("process", "").split(":")[-1]
                if process in ("BLOOD_CH", "NEO_PREPARE", "NEO_PREDICT"):
                    terminal[(row.get("tag"), process)] = row
    cases = []
    for case in plan["cases"]:
        record = dict(case_id=case["case_id"])
        for feature in ("blood_ch", "rna_neoantigens"):
            status = dict(case[feature])
            path = result_root / case["case_id"] / feature / "status.json"
            if feature == "rna_neoantigens" and not path.is_file():
                path = result_root / case["case_id"] / "neo_export/status.json"
            if path.is_file():
                status = json.loads(path.read_text())
                status["status_file"] = str(path)
                if status["status"] == "PREPARED":
                    status.update(status="INCOMPLETE", reason="Export finished; prediction has no successful completion receipt")
            elif status["status"] == "READY":
                status.update(status="INCOMPLETE", reason="Selected task has no successful completion receipt")
            for process in (("BLOOD_CH",) if feature == "blood_ch" else ("NEO_PREPARE", "NEO_PREDICT")):
                row = terminal.get((case["case_id"], process), {})
                if row.get("status") in ("FAILED", "ABORTED"):
                    folder = {"BLOOD_CH": "blood_ch", "NEO_PREPARE": "neo_export", "NEO_PREDICT": "rna_neoantigens"}[process]
                    task_status = Path(row["workdir"]) / folder / "status.json"
                    status = json.loads(task_status.read_text()) if task_status.is_file() else dict(status="FAILED", reason="Nextflow task failed before a completion receipt")
                    status.update(task_workdir=row["workdir"], task_exit=row.get("exit"), process=process)
            record[feature] = status
        cases.append(record)
    return dict(schema_version=SCHEMA, updated_at=timestamp(), cases=cases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("plan")
    prepare.add_argument("--config", type=Path, required=True)
    prepare.add_argument("--out", type=Path, required=True)
    status = sub.add_parser("status")
    status.add_argument("--plan", type=Path, required=True)
    status.add_argument("--results", type=Path, required=True)
    status.add_argument("--trace", type=Path)
    status.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "plan":
        write_json(args.out, build_plan(json.loads(args.config.read_text())))
    else:
        write_json(args.out, collect_status(json.loads(args.plan.read_text()), args.results, args.trace))


if __name__ == "__main__":
    main()
