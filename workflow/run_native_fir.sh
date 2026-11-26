#!/usr/bin/env bash
# Owner-operated submission of a prepared original-FASTQ native case.
set -euo pipefail
manifest= config= blood_ch=false resume=false
while (($#)); do
    case "$1" in
        --manifest) manifest="${2:?CSV required}"; shift 2 ;;
        --config) config="${2:?site.env required}"; shift 2 ;;
        --blood-ch) blood_ch=true; shift ;;
        --resume) resume=true; shift ;;
        -h|--help)
            echo 'Usage: run_native_fir.sh --manifest samplesheet.csv --config site.env [--blood-ch] [--resume]'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
manifest="$(readlink -f "${manifest:?Use --manifest}")"
config="$(readlink -f "${config:?Use --config}")"
source "$config"
run_dir="$RUN_ROOT/${WIGITS_RUN_ID:?Set WIGITS_RUN_ID in the supplied configuration}"
package="$run_dir/package"
mkdir -p "$package" "$run_dir/logs" "$run_dir/state"

# Resume reuses the same manifest, release and selected options.
copy_once() {
    if [[ -e "$2" ]]; then cmp -s "$1" "$2" || { echo "Existing package differs: $2" >&2; exit 2; }
    else cp "$1" "$2"; fi
}
copy_once "$manifest" "$package/samplesheet.csv"
copy_once "$config" "$package/site.env"
if [[ -e "$package/workflow" || -L "$package/workflow" ]]; then
    [[ "$(readlink -f "$package/workflow")" == "$(readlink -f "$WIGITS_WORKFLOW_DIR")" ]] || {
        echo 'Existing package uses a different release' >&2; exit 2;
    }
else ln -s "$WIGITS_WORKFLOW_DIR" "$package/workflow"; fi
if $blood_ch; then
    copy_once "${WIGITS_BLOOD_CH_CONFIG:?Bind WIGITS_BLOOD_CH_CONFIG in site.env}" "$package/optional_extensions.json"
    python3 - "$package/optional_extensions.json" <<'PY'
import json, sys
d=json.load(open(sys.argv[1]))
if not d['cases'] or any(c['options'].get('blood_ch') is not True or
    c['laboratory'].get('normal_specimen_type','').strip().casefold() != 'blood'
    for c in d['cases']):
    raise SystemExit('Blood CH requires a selected branch and a recorded blood normal')
PY
elif [[ -f "$package/optional_extensions.json" ]]; then
    echo 'This package already selects extensions; repeat its original options' >&2; exit 2
fi
if [[ -f "$run_dir/state/head_job.txt" ]]; then
    $resume || { echo 'Already submitted; use --resume only after the previous head ends' >&2; exit 2; }
    previous="$(cat "$run_dir/state/head_job.txt")"
    if [[ -n "$(squeue -h -j "$previous" -o '%i')" ]]; then
        echo "Previous head $previous is still active" >&2; exit 2
    fi
fi
job="$(sbatch --parsable --account="$SLURM_ACCOUNT_CPU" --time=6-23:00:00 \
    --job-name="$WIGITS_RUN_ID" --chdir="$run_dir" \
    --output="$run_dir/logs/head-%j.out" --error="$run_dir/logs/head-%j.err" \
    --export="ALL,WIGITS_SITE_ENV=$package/site.env,WIGITS_WORKFLOW_DIR=$WIGITS_WORKFLOW_DIR" \
    "$WIGITS_WORKFLOW_DIR/run_full_native.sbatch" "$package" "$WIGITS_RUN_ID")"
printf '%s\n' "${job%%;*}" > "$run_dir/state/head_job.txt"
echo "Submitted ${job%%;*}"
echo "Run: $run_dir"
