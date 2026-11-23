#!/usr/bin/env bash
# Build an isolated, versioned source checkout from a pinned local upstream.
# No shared Nextflow cache or installed container is changed.
set -euo pipefail
native_repo=${1:?upstream local repository required}
destination=${2:?new destination required}
package_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
upstream=7c74c87a43749952b38c9a18915947570f0595a0
test ! -e "$destination"
git clone --shared --no-checkout "$native_repo" "$destination"
git -C "$destination" -c checkout.workers=8 checkout --detach "$upstream"
git -C "$destination" apply "$package_dir/patches/oncoanalyser-3.0.0-qc-report.patch"
cp -R "$package_dir/overlays/oncoanalyser-3.0.0/." "$destination/"
git -C "$destination" add subworkflows/local/orange_reporting/main.nf modules/local/orange/main.nf nextflow.config lib/OrangeReportEligibility.groovy assets/orange_qc_report
git -C "$destination" diff --cached --check
git -C "$destination" commit -m 'fix: render QC-failed ORANGE reports with unavailable circos plots'
git -C "$destination" tag v3.0.0-qlmp.1
git -C "$destination" rev-parse HEAD
