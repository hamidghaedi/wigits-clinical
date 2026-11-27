# Releases and upstream Oncoanalyser updates

Owner direction November23,2026: maintain reviewed work on `main`, publish
numbered releases, and use short-lived branches for changes. The original
`hartwig_adoptation` branch is retained as historical source; its commits are
part of the first consolidated release. `main` is not a deployment instruction
for existing runs.

## One maintained source, frozen executions

1. Branch from current `main` for the change. Keep experimental routes explicitly
   selected and the native routine entry unchanged unless the owner authorizes
   a route change.
2. Run the applicable offline tests and source/config checks. Include the exact
   commands, outcomes and relevant prior HPC evidence in the pull request.
3. Update `VERSION`, `VERSIONS.md`, current instructions and release notes.
   Use a new numbered tag; never move an existing published tag.
4. Merge the reviewed pull request to `main`, preserving its published history.
   Tag the merged commit and push the tag. The Git tag plus committed release
   notes are the release; a GitHub release page may additionally mirror them.
5. New execution packages come from that tag/exact commit and record the actual
   workflow, containers, resource binding and selected route. Existing heads,
   work/cache and frozen packages retain their original versions. Deploying a
   new release is distinct from merely merging its source.

Current routine original-FASTQ entry:
[`run_full_native.sbatch`](../workflow/run_full_native.sbatch), using a frozen
package prepared according to [FULL_NATIVE_EXECUTION.md](FULL_NATIVE_EXECUTION.md).
The coordinator's maintained staged queue prepares the existing packages;
merging this source does not replace its deployed runtime. Keep all known
lanes/mates and existing tumour/normal bindings. No new bulk counting or
hashing gates are introduced.

## Current release

[v0.3.0](releases/v0.3.0.md) adds independently optional blood CH and
RNA-informed neoantigen execution, plus explicit Fir site selection and an
owner-operated launcher. The owner explicitly authorized direct Git merge/tag/
push for this promotion on November26 after the passing real-data tests. This
is a separate one-time exception; future changes retain the normal PR policy.

## Previous consolidation

[v0.2.0](releases/v0.2.0.md) combines the lab-qc/v0.4.1 operational work and
input-v5 RNA sizing without changing upstream source/resource or native caller
bindings. The source review, tests and publication state are recorded there.
Fir site configuration/deployment is a separate step.

## When upstream publishes a new version

Use a branch such as `codex/update-oncoanalyser-X.Y.Z` from our current `main`.
Our repository is a wrapper around separately pinned oncoanalyser source;
do not merge the upstream repository's entire history into this repository.

1. Read the official release notes, input/schema changes, native module/container
   bindings and compatible reference-bundle instructions. Check the coordinator's
   `docs/HARTWIG_RELEASE_POLICY.md` and maintained release checker for the new
   package. Record URLs/date and exact versions; never silently choose a moving
   upstream development branch.
2. Stage the new source/resources under new versioned directories, reusing
   compatible existing assets where documented. Keep earlier versions needed
   by running or reproducible historical attempts. Site paths remain in the
   existing site configuration/directory contract.
3. Review every local overlay/config against the new upstream source. Remove
   a local fix when upstream provides it. Reapply only changes still needed,
   as an explicit versioned patch/overlay. The current inventory is below.
4. Run offline tests plus a small workflow test, then compare a representative
   full case with the preceding release. Assess inputs, caller FILTER/counts,
   variant retention, SV/CN/QC, source-value semantics and report generation.
   Expected analytical differences must be explained; agreement is not itself
   validation. Record limitations instead of claiming universal equivalence.
5. Update the bindings and publish a new numbered wrapper release through a
   pull request. Only new authorized packages select it. Existing runs are
   never repointed to a different source/reference mid-execution.

| Local change | Current scope / upgrade review |
|---|---|
| `prepare_oncoanalyser_source_values.py`, `overlays/orange_source_values/` | Default3.0.0-qlmp.2.2: distinguish unavailable TINC and compute deletion extents from source coordinates. Check native implementation before carrying it forward. |
| `conf/lane_alignment_resources.config`, `build_input_resources.py` | Site resource sizing; verify upstream process names and task input shapes. Preserve lane scope and BWA64CPUs. |
| `conf/rna_duplicate_marking.config` | Explicitly disables optical comparisons for RNA while retaining duplicate flags. Recheck upstream Picard invocation. |
| `conf/esvee_runtime_recovery.config`, `esvee_jvm_recovery.env` | Diagnosed runtime recovery only; recheck JVM/HTSJDK versions before applying to a new source. |
| `patches/gridss-2.13.2-kraken-long-counts.patch` | Explicit overflow recovery only; remove when the selected upstream implementation handles large counts. |
| `sage_compat/`, `patches/sage-5.0.2-evidence-export.patch` | Experimental DeepSomatic substitution/native context. Verify every downstream consumed field; never synthesize missing measurements. |
| `patches/oncoanalyser-3.0.0-qc-report.patch`, older renderer overlay | Explicit historical QC-failure rendering experiment; not the default source-value release or permission to bypass scientific QC. |

## Existing test entry points

Use an existing Python with installed PyYAML. Do not use `python -S` for the
complete suite: it hides PyYAML and prevents the YAML-intake module loading.

```bash
python -m unittest discover -s tests -v
python -m unittest discover -s workflow/sage_compat/tests -v
git diff --check
```

Run `bash -n` over tracked `.sh`/`.sbatch` sources. The Java source-value test
uses the already pinned ORANGE image under Slurm as documented in
[ORANGE_SOURCE_VALUES.md](ORANGE_SOURCE_VALUES.md). Skipped integration tests
remain explicitly skipped unless their existing tool/runtime prerequisites
are available; do not weaken them to manufacture a green result.

The [Nextflow revision documentation](https://docs.seqera.io/nextflow/cli)
describes `-r` for remote pipeline releases/commits. Our installed-source launch
uses the exact versioned local checkout instead; record its SHA explicitly.

## Reusable native-intake component

`native-intake/v1.1.5` supplies the sample-independent laboratory manifest
preparer, operator guide and owner-run QC enrollment on the v0.3.0 native base.
The maintained source is `codex/reusable-native-intake`; main promotion still
follows the normal PR policy. Issued tool packages pin this component tag and
its exact commit, not a moving branch. No new scientific workflow or per-case
implementation is introduced. The wrapper's native tool/reference/resource
bindings remain those in v0.3.0. See [manifest preparation](NATIVE_MANIFEST_PREPARATION.md).
The [v1.1.5 fix](releases/native-intake-v1.1.5.md) adds optional-prediction recovery
without invoking native analysis or blood CH; previous executions remain pinned.
