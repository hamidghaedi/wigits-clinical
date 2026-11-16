#!/usr/bin/env python3
"""Build display-only oncoanalyser 3.0.0-qlmp.2.2 from pinned, unmodified upstream.

Native analytical code, report eligibility and chapter implementations are retained.
Only the ORANGE invocation gains source-value normalization and native JSON/receipts.
"""
import argparse
from pathlib import Path
import shutil
import subprocess

UPSTREAM = '7c74c87a43749952b38c9a18915947570f0595a0'
VERSION = '3.0.0-qlmp.2.2'


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Pinned upstream source differs at ' + repr(old))
    return text.replace(old, new)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('upstream_repo', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    if args.destination.exists():
        raise SystemExit('Destination already exists')
    subprocess.run(['git', 'clone', '--shared', '--no-checkout', str(args.upstream_repo), str(args.destination)], check=True)
    def git(*command):
        subprocess.run(['git', '-C', str(args.destination), *command], check=True)
    git('checkout', '--detach', UPSTREAM)
    overlay = Path(__file__).resolve().parent / 'overlays/orange_source_values'
    shutil.copytree(overlay, args.destination / 'assets/orange_source_values')
    module = args.destination / 'modules/local/orange/main.nf'
    text = module.read_text()
    text = replace_once(text, "    def log_level_arg", "    def source_values_java = file(\"${projectDir}/assets/orange_source_values/SourceValues.java\").text\n    def corrected_orange_java = file(\"${projectDir}/assets/orange_source_values/CorrectedOrangeApplication.java\").text\n\n    def log_level_arg")
    prefix = '''    cat > SourceValues.java <<'QLMP_SOURCE_VALUES'
${source_values_java}
QLMP_SOURCE_VALUES
    cat > CorrectedOrangeApplication.java <<'QLMP_CORRECTED_ORANGE'
${corrected_orange_java}
QLMP_CORRECTED_ORANGE
    mkdir -p source_values_classes
    javac --release 17 -cp /usr/local/share/hmftools-orange-5.0.1-0/orange.jar -d source_values_classes SourceValues.java CorrectedOrangeApplication.java
    java '''
    text = replace_once(text, '    orange ', prefix)
    text = replace_once(text, '        ${args}', '        -cp source_values_classes:/usr/local/share/hmftools-orange-5.0.1-0/orange.jar org.qlmp.orange.CorrectedOrangeApplication \\\\\n        ${args}')
    text = replace_once(text, "    path 'versions.yml'", "    tuple val(meta), path('output/*.source_values.json'), topic: orange_source_values\n    tuple val(meta), path('output/*.orange.source.json'), topic: orange_source_json\n    path 'versions.yml'")
    text = replace_once(text, '        orange: ', '        qlmp-orange-source-values: 0.1.0\n        orange: ')
    # Inserted Java has unindented lines. Keep the shell's versions terminator
    # unindented too, independent of Nextflow's multiline-string indentation.
    text = replace_once(text, '    END_VERSIONS\n', 'END_VERSIONS\n')
    module.write_text(text)
    outputs = args.destination / 'subworkflows/local/prepare_outputs/main.nf'
    text = outputs.read_text()
    old = "            channel.topic('orange_json').map { meta, d ->                      return [\"${meta.key}/orange/${d.name}\", d] },"
    if text.count(old) != 2:
        raise ValueError('Expected native WGTS and targeted ORANGE publication entries')
    text = text.replace(old, old + '\n'
        + "            channel.topic('orange_source_values').map { meta, d ->             return [\"${meta.key}/orange/${d.name}\", d] },\n"
        + "            channel.topic('orange_source_json').map { meta, d ->               return [\"${meta.key}/orange/${d.name}\", d] },")
    outputs.write_text(text)
    config = args.destination / 'nextflow.config'
    config.write_text(replace_once(config.read_text(), "version         = '3.0.0'", "version         = '" + VERSION + "'"))
    git('add', 'modules/local/orange/main.nf', 'subworkflows/local/prepare_outputs/main.nf', 'assets/orange_source_values', 'nextflow.config')
    git('diff', '--cached', '--check')
    git('commit', '-m', 'Correct ORANGE deletion extents and retain source-aware TINC availability')
    git('tag', 'v' + VERSION)
    git('rev-parse', 'HEAD')


if __name__ == '__main__':
    main()
