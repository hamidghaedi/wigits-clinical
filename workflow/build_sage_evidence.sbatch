#!/usr/bin/env bash
# Arguments: frozen source/patch directory, new build directory, SAGE image.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
package=$1
out=$2
image=$3
mkdir "$out"
cp -R "$package/source" "$out/source"
patch -d "$out/source" -p1 < "$package/sage-5.0.2-evidence-export.patch"
apptainer exec --cleanenv --bind /global "$image" bash -c '
    set -euo pipefail
    out=$1
    original=$(dirname "$(readlink -f "$(command -v sage)")")/sage.jar
    cp "$original" "$out/sage-5.0.2-evidence.jar"
    mkdir "$out/classes"
    javac -cp "$original" -d "$out/classes" \
        "$out/source/com/hartwig/hmftools/sage/evidence/ReadContextCounters.java" \
        "$out/source/com/hartwig/hmftools/sage/pipeline/RegionTask.java"
    jar uf "$out/sage-5.0.2-evidence.jar" -C "$out/classes" .
    sha256sum "$original" "$out/sage-5.0.2-evidence.jar" > "$out/jar-sha256.txt"
' bash "$out"
