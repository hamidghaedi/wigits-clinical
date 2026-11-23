#!/usr/bin/env bash
# Build small interface fixtures from a real oncoanalyser case (default: the subject_a test run).
# Keeps header lines plus a bounded number of records so tests run without the full outputs.
set -euo pipefail
SRC="${1:-/global/scratch/hpc6049/oncoanalyser/test_3.0.0/results/subject_a}"
DST="${2:-$(dirname "$0")/../tests/fixtures/neo_interface/subject_a}"
N="${3:-200}"
mkdir -p "$DST"/{purple,sage_append,lilac,isofox,neo/finder,neo/scorer,linx/somatic_annotations,orange}
copy_vcf() { # keep header + first N records
  local src="$1" dst="$2"
  { zcat "$src" | grep "^#"; zcat "$src" | grep -v "^#" | head -n "$N"; } | gzip > "$dst"
}
for f in "$SRC"/purple/*.purple.somatic.vcf.gz;  do copy_vcf "$f" "$DST/purple/$(basename "$f")"; done
for f in "$SRC"/purple/*.purple.germline.vcf.gz; do copy_vcf "$f" "$DST/purple/$(basename "$f")"; done
cp "$SRC"/purple/*.purple.purity.tsv "$DST/purple/"
cp "$SRC"/purple/*.purple.qc "$DST/purple/" 2>/dev/null || true
for d in "$SRC"/sage_append/*/; do
  s=$(basename "$d"); mkdir -p "$DST/sage_append/$s"
  for f in "$d"*.sage.append.vcf.gz; do copy_vcf "$f" "$DST/sage_append/$s/$(basename "$f")"; done
done
cp "$SRC"/lilac/*.lilac.tsv "$SRC"/lilac/*.lilac.qc.tsv "$DST/lilac/" 2>/dev/null || true
for f in "$SRC"/isofox/*.isf.gene_data.tsv "$SRC"/isofox/*.isf.transcript_data.tsv; do
  [ -e "$f" ] && head -n "$N" "$f" > "$DST/isofox/$(basename "$f")"
done
for f in "$SRC"/neo/finder/*.neo_data.tsv;            do [ -e "$f" ] && head -n "$N" "$f" > "$DST/neo/finder/$(basename "$f")"; done
for f in "$SRC"/neo/scorer/*.neo.peptide_scores.tsv;  do [ -e "$f" ] && head -n "$N" "$f" > "$DST/neo/scorer/$(basename "$f")"; done
for f in "$SRC"/neo/scorer/*.neo.neoepitope.tsv;      do [ -e "$f" ] && head -n "$N" "$f" > "$DST/neo/scorer/$(basename "$f")"; done
for f in "$SRC"/linx/somatic_annotations/*.linx.neoepitope.tsv; do [ -e "$f" ] && head -n "$N" "$f" > "$DST/linx/somatic_annotations/$(basename "$f")"; done
# ORANGE: keep only the immuneEscape block and identifiers (the full json is large and case-identifying)
python3 - "$SRC" "$DST" <<'PY'
import glob, json, os, sys
src, dst = sys.argv[1], sys.argv[2]
hits = glob.glob(os.path.join(src, "orange", "*.orange.json"))
if hits:
    d = json.load(open(hits[0]))
    keep = {k: d[k] for k in ("sampleId", "refGenomeVersion", "pipelineVersion", "experimentType", "immuneEscape") if k in d}
    json.dump(keep, open(os.path.join(dst, "orange", os.path.basename(hits[0])), "w"), indent=1)
PY
echo "fixtures written to $DST"
du -sh "$DST"
