#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -eo pipefail
# Arguments: frozen source directory, new build output, original GRIDSS image.
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
package=$1
out=$2
image=$3
mkdir "$out"
cp -R "$package/source" "$out/source"
patch -d "$out/source" -p1 < "$package/gridss-2.13.2-kraken-long-counts.patch"
apptainer exec --cleanenv --bind /global "$image" bash -c '
    set -euo pipefail
    out=$1; package=$2
    original=/env/share/gridss-2.13.2-3/gridss.jar
    cp "$original" "$out/gridss-2.13.2-long-counts.jar"
    mkdir "$out/classes" "$out/test"
    javac --release 8 -d "$out/classes" "$out/source/au/edu/wehi/idsv/kraken/KrakenReportLine.java"
    jar uf "$out/gridss-2.13.2-long-counts.jar" -C "$out/classes" au/edu/wehi/idsv/kraken/KrakenReportLine.class
    javac --release 8 -cp "$out/gridss-2.13.2-long-counts.jar" -d "$out/test" "$package/KrakenLongCountTest.java"
    if java -cp "$out/test:$original" KrakenLongCountTest > "$out/original-regression.log" 2>&1; then
        echo "Expected original overflow was not reproduced" >&2; exit 1
    fi
    grep -q NumberFormatException "$out/original-regression.log"
    java -cp "$out/test:$out/gridss-2.13.2-long-counts.jar" KrakenLongCountTest
    # Picard help exits 1 even when the command loads normally.
    help_status=0
    java -cp "$out/gridss-2.13.2-long-counts.jar" gridss.kraken.IdentifyViralTaxa --help > "$out/tool-help.log" 2>&1 || help_status=$?
    test "$help_status" -le 1
    grep -q INPUT_KRAKEN2_REPORT "$out/tool-help.log"
    sha256sum "$original" "$out/gridss-2.13.2-long-counts.jar" > "$out/jar-sha256.txt"
' bash "$out" "$package"
