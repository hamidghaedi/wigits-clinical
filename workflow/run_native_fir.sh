#!/usr/bin/env bash
# Owner-operated submission of a prepared original-FASTQ native case.
set -euo pipefail
manifest= config= controller_time= dependency= blood_ch=false resume=false
entry=native original_run= recovery_workflow= existing_predictions= existing_export=
while (($#)); do
    case "$1" in
        --manifest) manifest="${2:?CSV required}"; shift 2 ;;
        --config) config="${2:?site.env required}"; shift 2 ;;
        --entry) entry="${2:?native, neo-prepare, neo-predict or neo-review required}"; shift 2 ;;
        --run-dir) original_run="${2:?original native run root required}"; shift 2 ;;
        --workflow-dir) recovery_workflow="${2:?pinned workflow directory required}"; shift 2 ;;
        --existing-export) existing_export="${2:?prepared export directory required}"; shift 2 ;;
        --existing-predictions) existing_predictions="${2:?completed prediction directory required}"; shift 2 ;;
        --dependency) dependency="${2:?Slurm dependency required}"; shift 2 ;;
        --controller-time) controller_time="${2:?Slurm walltime required}"; shift 2 ;;
        --blood-ch) blood_ch=true; shift ;;
        --resume) resume=true; shift ;;
        -h|--help)
            echo 'Usage: run_native_fir.sh --manifest samplesheet.csv --config site.env [--dependency TYPE:JOBID] [--controller-time D-HH:MM:SS] [--blood-ch] [--resume]'
            echo 'Or: run_native_fir.sh --entry neo-prepare|neo-predict|neo-review --run-dir ORIGINAL_RUN --workflow-dir PINNED_WORKFLOW [--existing-export DIR] [--existing-predictions DIR] [--dependency TYPE:JOBID] [--controller-time D-HH:MM:SS] [--resume]'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
case "$entry" in
    native)
        [[ -z "$original_run$recovery_workflow$existing_predictions$existing_export" ]] || { echo 'Recovery arguments require an optional entry' >&2; exit 2; }
        manifest="$(readlink -f "${manifest:?Use --manifest}")"
        config="$(readlink -f "${config:?Use --config}")" ;;
    neo-prepare|neo-predict|neo-review)
        [[ -z "$manifest$config" && "$blood_ch" == false ]] || { echo 'Optional entries reuse --run-dir package/site.env and options' >&2; exit 2; }
        original_run="$(readlink -f "${original_run:?Use --run-dir}")"
        recovery_workflow="$(readlink -f "${recovery_workflow:?Use --workflow-dir for the new pinned source}")"
        config="$original_run/package/site.env"
        if [[ "$entry" == neo-prepare ]]; then
            [[ -z "$existing_export" ]] || { echo '--existing-export is for prediction/review' >&2; exit 2; }
        elif [[ -n "$existing_export" ]]; then
            existing_export="$(readlink -f "$existing_export")"
        fi
        if [[ "$entry" == neo-review ]]; then
            existing_predictions="$(readlink -f "${existing_predictions:?Use --existing-predictions}")"
        else
            [[ -z "$existing_predictions" ]] || { echo '--existing-predictions is for neo-review' >&2; exit 2; }
        fi ;;
    *) echo "Unknown entry: $entry" >&2; exit 2 ;;
esac
source "$config"
# CLI overrides the prepared site value; old packages use the measured default.
if [[ "$entry" == native ]]; then
    controller_time="${controller_time:-${WIGITS_CONTROLLER_TIME:-4-00:00:00}}"
else
    controller_time="${controller_time:-1-00:00:00}"
fi
# Apply submission dependencies to this head only, never nested sbatch calls.
dependency="${dependency:-${SBATCH_DEPENDENCY:-}}"
unset SBATCH_DEPENDENCY
dependency_args=()
if [[ -n "$dependency" ]]; then
    dependency_args=(--dependency="$dependency")
fi
if [[ "$entry" != native ]]; then
    recovery="$original_run/extensions/recovery/native-intake-v1.3.0/$entry"
    mkdir -p "$original_run/logs" "$recovery/state"
    if [[ -f "$recovery/state/head_job.txt" ]]; then
        $resume || { echo 'Already submitted; use --resume after the previous recovery head ends' >&2; exit 2; }
        previous="$(cat "$recovery/state/head_job.txt")"
        [[ -z "$(squeue -h -j "$previous" -o '%i')" ]] || { echo "Recovery head $previous is still active" >&2; exit 2; }
    fi
    job="$(sbatch --parsable --account="$SLURM_ACCOUNT_CPU" --time="$controller_time" "${dependency_args[@]}" \
        --job-name="${original_run##*/}-$entry" --chdir="$original_run" \
        --output="$original_run/logs/$entry-%j.out" --error="$original_run/logs/$entry-%j.err" \
        --export="ALL,WIGITS_RECOVERY_WORKFLOW=$recovery_workflow" \
        "$recovery_workflow/run_optional_recovery.sbatch" "$original_run" "$entry" "$existing_predictions" "$existing_export")"
    printf '%s\n' "${job%%;*}" > "$recovery/state/head_job.txt"
    printf '%s\n' "$controller_time" > "$recovery/state/head_walltime.txt"
    printf '%s\n' "${dependency:-none}" > "$recovery/state/head_dependency.txt"
    printf '%s\n' "$entry" > "$recovery/state/entry.txt"
    echo "Submitted ${job%%;*}"
    echo "Entry: $entry"
    echo "Run: $recovery"
    exit 0
fi
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
for name in input_plan.json release.json case_manifest.json laboratory.json clinical_context.json native_runtime.config; do
    [[ ! -f "$(dirname "$manifest")/$name" ]] || copy_once "$(dirname "$manifest")/$name" "$package/$name"
done
if [[ -e "$package/workflow" || -L "$package/workflow" ]]; then
    [[ "$(readlink -f "$package/workflow")" == "$(readlink -f "$WIGITS_WORKFLOW_DIR")" ]] || {
        echo 'Existing package uses a different release' >&2; exit 2;
    }
else ln -s "$WIGITS_WORKFLOW_DIR" "$package/workflow"; fi
extensions="${WIGITS_EXTENSIONS_CONFIG:-}"
if $blood_ch; then
    extensions="${extensions:-${WIGITS_BLOOD_CH_CONFIG:?Prepared extension configuration required}}"
fi
if [[ -n "$extensions" ]]; then
    copy_once "$extensions" "$package/optional_extensions.json"
    python3 - "$package/optional_extensions.json" "$blood_ch" <<'PY'
import json, sys
d=json.load(open(sys.argv[1]))
if not d['cases'] or any((sys.argv[2]=='true' and c['options'].get('blood_ch') is not True) or
    (c['options'].get('blood_ch') and c['laboratory'].get('normal_specimen_type','').strip().casefold() != 'blood') for c in d['cases']):
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
native_args=()
if [[ -n "${WIGITS_NATIVE_EXTRA_CONFIG:-}" ]]; then
    native_args=(-c "$WIGITS_NATIVE_EXTRA_CONFIG")
fi
job="$(sbatch --parsable --account="$SLURM_ACCOUNT_CPU" --time="$controller_time" "${dependency_args[@]}" \
    --job-name="$WIGITS_RUN_ID" --chdir="$run_dir" \
    --output="$run_dir/logs/head-%j.out" --error="$run_dir/logs/head-%j.err" \
    --export="ALL,WIGITS_SITE_ENV=$package/site.env,WIGITS_WORKFLOW_DIR=$WIGITS_WORKFLOW_DIR" \
    "$WIGITS_WORKFLOW_DIR/run_full_native.sbatch" "$package" "$WIGITS_RUN_ID" "${native_args[@]}")"
printf '%s\n' "${job%%;*}" > "$run_dir/state/head_job.txt"
printf '%s\n' "$controller_time" > "$run_dir/state/head_walltime.txt"
printf '%s\n' "${dependency:-none}" > "$run_dir/state/head_dependency.txt"
echo "Submitted ${job%%;*}"
echo "Controller walltime: $controller_time"
echo "Dependency: ${dependency:-none}"
echo "Run: $run_dir"
