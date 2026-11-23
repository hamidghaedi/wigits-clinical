#!/usr/bin/env bash
#SBATCH --job-name=hartwig-metadata-check
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 java/21.0.1 apptainer
set -u
source_java=$1
fixture=$2
outdir=$3
neo_image=$4
orange_image=$5
mkdir -p "$outdir" "$SLURM_TMPDIR/classes"
javac --release 17 -d "$SLURM_TMPDIR/classes" "$source_java"
apptainer exec --bind "$SLURM_TMPDIR:$SLURM_TMPDIR" "$neo_image" cp /usr/local/share/hmftools-neo-1.3-0/neo.jar "$SLURM_TMPDIR/neo.jar"
java -Xmx1g -cp "$SLURM_TMPDIR/neo.jar:$SLURM_TMPDIR/classes" CheckCancerType "$fixture" \
    > "$outdir/neo_native_probe.log"
apptainer exec --bind "$SLURM_TMPDIR:$SLURM_TMPDIR" "$orange_image" cp /usr/local/share/hmftools-orange-5.0.1-0/orange.jar "$SLURM_TMPDIR/orange.jar"
javap -p -c -classpath "$SLURM_TMPDIR/orange.jar" com.hartwig.hmftools.orange.algo.OrangeAlgo \
    > "$outdir/orange_native_bytecode.txt"
sha256sum "$SLURM_TMPDIR/neo.jar" "$SLURM_TMPDIR/orange.jar" > "$outdir/native_jar_sha256.txt"
echo 'Native metadata checks completed; no genomic analysis or report rendering performed.'
