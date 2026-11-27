#!/usr/bin/env python3
"""Execute the completed feature implementations from optional Nextflow tasks."""
import argparse
import json
import os
import re
from pathlib import Path
import subprocess
import sys

from optional_extensions import resource_allowance, timestamp, write_json
from neo_phasing import summarise_usage

HERE = Path(__file__).resolve().parent


class Commands:
    def __init__(self, resources, outdir):
        self.resources, self.outdir = resources, Path(outdir).resolve()
        self.calls = []

    def run(self, argv, env=None):
        self.calls.append([str(arg) for arg in argv])
        write_json(self.outdir / "commands.json", self.calls)
        subprocess.run(self.calls[-1], check=True, env=env)

    def container(self, image, *args):
        binds = list(self.resources["bind_roots"]) + [str(self.outdir), os.environ.get("SLURM_TMPDIR", str(self.outdir))]
        self.run(["apptainer", "exec", "--cleanenv", "--env", "PYTHONNOUSERSITE=1",
                  "--env", "PYTHONPATH=", "-B", ",".join(dict.fromkeys(binds)), self.resources[image], *args])


def sort_annotation_vcf(path):
    """Order task-owned VEP output for tabix without changing record contents.

    Somatic and germline proximal records arrive as separate ordered groups.
    Sorting after VEP also covers completed exports reused by recovery runs.
    """
    path = Path(path)
    lines = path.read_text().splitlines(keepends=True)
    headers = [line for line in lines if line.startswith('#')]
    records = [line for line in lines if not line.startswith('#')]
    contigs = [re.search(r'ID=([^,>]+)', line).group(1)
               for line in headers if line.startswith('##contig=<ID=')]
    ranks = {chrom: index for index, chrom in enumerate(contigs)}

    def key(line):
        chrom, pos = line.split('\t', 2)[:2]
        return ranks.get(chrom, len(ranks)), chrom, int(pos)

    records.sort(key=key)
    temporary = path.with_suffix(path.suffix + '.sorted')
    temporary.write_text(''.join(headers + records))
    temporary.replace(path)


def chip(case, resources, out, threads, commands):
    env = dict(os.environ, THREADS=str(threads), REF=resources["reference_fasta"], RES=resources["hmf_resources"],
               SAGE_SIF=resources["sage_image"], PAVE_SIF=resources["pave_image"], SAMTOOLS_SIF=resources["samtools_image"],
               BIND=",".join(resources["bind_roots"] + [str(out)]))
    args = ["bash", HERE / "chip_call.sh", case["case_id"], case["normal_sample"], case["tumour_sample"], case["results_case_dir"], out]
    if case.get("mean_blood_depth") is not None:
        args.append(str(case["mean_blood_depth"]))
    commands.run(args, env)
    # chip_call stages native REDUX inputs as links for SAGE. These are task
    # inputs, not analytical output. publishDir(copy) follows links inside a
    # directory, so retire just these links before emitting the output folder.
    stage = out / "stage"
    if stage.is_dir():
        links = list(stage.iterdir())
        if any(not path.is_symlink() for path in links):
            raise ValueError("Unexpected non-link in CH input staging; preserve it and stop publication")
        for path in links:
            path.unlink()
        stage.rmdir()
    return dict(status="COMPLETED", reason="Dedicated blood SAGE/PAVE call, curation and review table completed")


def neo_prepare(case, resources, out, threads, commands):
    commands.run([sys.executable, HERE / "neo_export.py", "--results-case-dir", case["results_case_dir"], "--outdir", out, "--case", case["case_id"]])
    manifest = json.loads((out / (case["case_id"] + ".interface_manifest.json")).read_text())
    if manifest["hla"]["lilac_qc_status"] is None or "FAIL" in manifest["hla"]["lilac_qc_status"] or not manifest["hla"]["alleles"]:
        status = dict(status="NOT_ASSESSED", reason="Class-I HLA typing is absent or failed; no predictions made")
    elif manifest["counts"]["selected"] == 0:
        status = dict(status="NO_CANDIDATES", reason="Export assessed native somatic variants; no coding candidates selected")
    else:
        workload = manifest["counts"]["selected"] * len(manifest["hla"]["alleles"])
        status = dict(status="PREPARED", reason="Native RNA/DNA/HLA export ready", resources=resource_allowance("neo_predict", workload))
    write_json(out / "task_metadata.json", status)
    return status


def neo_predict(case, resources, out, threads, commands, prepared):
    sample = case["case_id"]
    prepared = Path(prepared).resolve()
    manifest = json.loads((prepared / (sample + ".interface_manifest.json")).read_text())
    # Keep the source export read-only. Native sequence alleles and measured
    # RNA support are consumed unchanged, including the single-sample proximal VCF.
    for name, source in (("selected", prepared / (sample + ".somatic.selected.vcf")),
                         ("proximal", prepared / (sample + ".proximal.vcf"))):
        target = out / (sample + "." + name + ".vep.vcf")
        commands.container("vep_image", "vep", "--input_file", source, "--output_file", target,
                           "--format", "vcf", "--vcf", "--force_overwrite", "--offline", "--cache",
                           "--dir_cache", resources["vep_cache"], "--cache_version", "115", "--merged", "--assembly", "GRCh38",
                           "--fasta", resources["reference_fasta"], "--use_given_ref", "--symbol", "--terms", "SO",
                           "--tsl", "--biotype", "--transcript_version", "--hgvs", "--flag_pick",
                           "--pick_order", "canonical,tsl,biotype,rank,ccds,length", "--plugin", "Wildtype", "--plugin", "Frameshift",
                           "--dir_plugins", resources["vep_plugins"], "--fork", str(threads))
        sort_annotation_vcf(target)
        commands.container("samtools_image", "bgzip", "-kf", target)
        commands.container("samtools_image", "tabix", "-fp", "vcf", str(target) + ".gz")
    prediction_dir = out / "pvacseq"
    coverage = ["--normal-cov", "5", "--tdna-cov", "10", "--trna-cov", "5", "--normal-vaf", "0.02", "--tdna-vaf", "0.05", "--trna-vaf", "0.0", "--expn-val", "0.0"]
    commands.container("pvactools_image", "pvacseq", "run", out / (sample + ".selected.vep.vcf.gz"), case["tumour_sample"],
                       ",".join(manifest["hla"]["alleles"]), "MHCflurry", "MHCflurryEL", "MHCnuggetsI", prediction_dir,
                       "--normal-sample-name", case["normal_sample"], "-e1", "8,9,10,11", "-p", out / (sample + ".proximal.vep.vcf.gz"),
                       "--n-threads", str(threads), "--pass-only", "--genes-of-interest-file", resources["genes_of_interest"], *coverage)
    results = prediction_dir / "MHC_Class_I"
    all_epitopes = results / (case["tumour_sample"] + ".all_epitopes.tsv")
    # Preserve the complete prediction product and the feature's explicit review
    # cascade. These filters never change the original somatic/report variant lists.
    binding, covered, supported, top = [results / (sample + "." + suffix + ".tsv") for suffix in ("binding", "coverage", "tsl", "top")]
    commands.container("pvactools_image", "pvacseq", "binding_filter", all_epitopes, binding, "-b", "500", "-p", "2")
    commands.container("pvactools_image", "pvacseq", "coverage_filter", binding, covered, *coverage)
    commands.container("pvactools_image", "pvacseq", "transcript_support_level_filter", covered, supported)
    commands.container("pvactools_image", "pvacseq", "top_score_filter", supported, top)
    phase_evidence = prepared / (sample + ".phasing_evidence.json")
    if phase_evidence.is_file():
        summarise_usage(prediction_dir,phase_evidence,out / "phasing_usage.json")
    write_json(out / "prediction_complete.json", dict(case_id=sample, completed_at=timestamp()))
    return neo_review(case, resources, out, threads, commands, prepared, out)


def neo_review(case, resources, out, threads, commands, prepared, prediction):
    sample = case["case_id"]
    results = Path(prediction).resolve() / "pvacseq/MHC_Class_I"
    aggregate = out / (sample + ".review.aggregated.tsv")
    # Preserve the feature's standalone aggregation, whose recorded metrics
    # use RNA VAF0.25, depth10 and allele-expression2.5. The repaired pVACseq
    # run now also produces an automatic aggregate using its discovery flags;
    # keep that as separate evidence instead of silently changing review labels.
    commands.container("pvactools_image", "pvacseq", "generate_aggregated_report",
                       results / (case["tumour_sample"] + ".all_epitopes.tsv"), aggregate,
                       # pVACtools multiplies this gene-expression argument by
                       # 2.5 for the allele-expression cutoff in its metrics.
                       "--trna-vaf", "0.25", "--trna-cov", "10", "--expn-val", "1.0")
    commands.run([sys.executable, HERE / "neo_candidates.py", "--results-case-dir", case["results_case_dir"],
                  "--pvacseq-dir", results, "--outdir", out, "--case", sample,
                  "--selected-vcf", Path(prepared) / (sample + ".somatic.selected.vcf"),
                  "--aggregate-report", aggregate, "--genes-of-interest-file", resources["genes_of_interest"]])
    return dict(status="COMPLETED", reason="RNA-informed independent class-I predictions and native NEO comparison completed")


def execute(plan, case_id, stage, outdir, threads, prepared=None, prediction=None):
    case = next(c for c in plan["cases"] if c["case_id"] == case_id)
    out = Path(outdir).resolve()
    if stage == "neo_predict" and (out == Path(prepared).resolve() or Path(prepared).resolve() in out.parents):
        raise ValueError("Prediction output must be separate from its preserved export")
    if stage == "neo_review" and (out == Path(prediction).resolve() or Path(prediction).resolve() in out.parents):
        raise ValueError("Review output must be separate from the preserved prediction directory")
    out.mkdir(parents=True, exist_ok=True)
    commands = Commands(plan["resources"], out)
    feature = "blood_ch" if stage == "blood_ch" else "rna_neoantigens"
    if case[feature]["status"] != "READY":
        raise ValueError("Refusing to execute a disabled/ineligible/unprepared feature")
    status = dict(case_id=case_id, feature=feature, started_at=timestamp(), status="RUNNING")
    write_json(out / "status.json", status)
    write_json(out / "source_plan.json", plan)
    try:
        if stage == "neo_predict":
            prepared_status = json.loads((Path(prepared) / "task_metadata.json").read_text())
            if prepared_status["status"] != "PREPARED":
                raise ValueError("Prediction requires a successful eligible export")
            export = json.loads((Path(prepared) / (case_id + '.interface_manifest.json')).read_text())
            samples = export['samples']
            if export['case'] != case_id or any(samples[key] != case[value] for key, value in
                    (('tumour_dna', 'tumour_sample'), ('normal_dna', 'normal_sample'), ('tumour_rna', 'rna_sample'))):
                raise ValueError("Export source disagrees with the configured case")
            result = neo_predict(case, plan["resources"], out, threads, commands, prepared)
        elif stage == "neo_review":
            source = Path(prediction)
            source_plan = json.loads((source / "source_plan.json").read_text())
            original = next(c for c in source_plan["cases"] if c["case_id"] == case_id)
            for key in ("tumour_sample", "normal_sample", "rna_sample", "results_case_dir", "native_genome"):
                if original.get(key) != case.get(key):
                    raise ValueError("Prediction source disagrees with the configured case")
            source_status = json.loads((source / "status.json").read_text())
            if source_status.get("status") != "COMPLETED" and not (source / "prediction_complete.json").is_file():
                raise ValueError("Prediction stage is not complete")
            result = neo_review(case, plan["resources"], out, threads, commands, prepared, prediction)
            result["reused_prediction_directory"] = str(source.resolve())
        else:
            result = {"blood_ch": chip, "neo_prepare": neo_prepare}[stage](case, plan["resources"], out, threads, commands)
        status.update(result)
    except Exception as exc:
        status.update(status="FAILED", reason=str(exc), completed_at=timestamp())
        write_json(out / "status.json", status)
        raise
    status["completed_at"] = timestamp()
    write_json(out / "status.json", status)
    write_json(out / "source_plan.json", plan)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--stage", choices=("blood_ch", "neo_prepare", "neo_predict", "neo_review"), required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--prediction", type=Path)
    args = parser.parse_args()
    execute(json.loads(args.plan.read_text()), args.case_id, args.stage, args.outdir, args.threads, args.prepared, args.prediction)


if __name__ == "__main__":
    main()
