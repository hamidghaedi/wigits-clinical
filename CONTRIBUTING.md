# Contributing

How work moves through this repository.

## Branches

`main` holds reviewed work only. Everything else happens on a branch named for what it does, for example `docs/clinical-layer`, `evidence/serve-directory`, `reporter/qc-fail-templates`. Branches merge to `main` through a pull request.

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

Every tool version this project depends on is pinned in `VERSIONS.md` once that file exists. Until then, state versions in the document or pull request that relies on them.
