# Contributing

How work moves through this repository.

## Branches

`main` holds reviewed work only. Everything else happens on a branch named for what it does, for example `docs/clinical-layer`, `evidence/serve-directory`, `reporter/qc-fail-templates`. Branches merge to `main` through a pull request.

**Owner direction November23: `main` is the maintained integration branch.**
The work consolidated in `hartwig_adoptation` is promoted through a pull request;
that branch remains historical rather than the default destination for new work.
Use short-lived branches from current `main` (for example
`codex/update-oncoanalyser-X.Y.Z`), run the relevant tests, and merge through a
pull request. Preserve existing commit history; do not rebase or force-push
published releases. A sibling DNA repository owns its own workflow and releases.

Production packages use numbered tags or exact commits, not the moving `main`
branch. Experimental DeepSomatic/BAM/QC-renderer paths remain explicitly selected;
merging their source does not enable them. See [release and upstream update
procedure](docs/RELEASING.md), [VERSION](VERSION) and [VERSIONS.md](VERSIONS.md).

## Where things are tested

Development and documentation happen here. Anything that needs real data, a reference bundle, containers or a scheduler runs on the HPC side. The hand-off is the pull request.

A pull request that touches `workflow/` or any runnable code should say in its description what needs to be run to test it, which inputs it expects, and what a pass looks like. The HPC side runs it and reports back on the pull request with the command used, the oncoanalyser and tool versions, and the outcome. Test outputs stay on the HPC side; only summaries and small diagnostic files come back into the repository, under `tests/results/<branch>/`.

## Commits

Plain commit messages that say what changed and why. One concern per commit where practical. No tooling trailers or generated signatures.

## Documentation

Markdown with Mermaid for diagrams. Every claim about Hartwig or OICR behaviour should point to a file and, where it helps, a line in the source it was read from, with the commit or tag named once at the top of the document.

## Gap register

When you find something missing, discontinued or private, add a row to `docs/gap_register.md` before writing code for it. When a gap is filled, change its status and link the commit or directory.

## Versions

`VERSION` identifies the repository-wide release. `VERSIONS.md` records its
upstream, local-source, reference and maintained component bindings. Existing
component tags retain their independent versions. Every execution package must
record the actual selected source/commit, containers, resource version and panel
hash; optional experiments retain their own explicit inputs and bindings.

## Reusable operations (owner direction November26)

Operational code and the human running guide must work across cases by changing
arguments and shared configuration. Do not write case-specific scripts, launch
helpers or operational guides. Generate case-specific samplesheets, option/site
bindings and provenance as runtime data from the maintained preparer; they are
inputs/receipts, never a second implementation. Laboratory QC retains each
complete original run manifest independently of the case selected for analysis.
