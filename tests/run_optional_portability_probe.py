#!/usr/bin/env python3
"""Compute-node probe of real container isolation and mixed-origin VCF indexing."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
from extension_tasks import Commands, sort_annotation_vcf


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pvactools-image', required=True)
    parser.add_argument('--samtools-image', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--vcf', type=Path, help='Optional small, uncompressed VEP VCF; source preserved')
    args = parser.parse_args()
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    resources = {'bind_roots':[str(root)], 'pvac':args.pvactools_image, 'samtools':args.samtools_image}
    calls = Commands(resources, root)
    version = subprocess.check_output(['apptainer','exec','--cleanenv',args.pvactools_image,'python','-c',
                                      'import sys; print("%d.%d" % sys.version_info[:2])'],text=True).strip()
    userbase = root/'userbase'
    site = userbase/'lib'/('python'+version)/'site-packages'
    site.mkdir(parents=True, exist_ok=True)
    (site/'numpy.py').write_text('raise RuntimeError("HOST_USER_SITE_WAS_IMPORTED")\n')
    previous = dict(os.environ)
    results = {}
    try:
        os.environ.update(APPTAINERENV_PYTHONUSERBASE=str(userbase), APPTAINERENV_PYTHONNOUSERSITE='')
        for mode in ('user_site','pythonpath'):
            if mode == 'pythonpath':
                os.environ['APPTAINERENV_PYTHONPATH'] = str(site)
            legacy = subprocess.run(['apptainer','exec','--cleanenv','-B',str(root),args.pvactools_image,
                                     'python','-c','import numpy'],text=True,capture_output=True)
            if legacy.returncode == 0 or 'HOST_USER_SITE_WAS_IMPORTED' not in legacy.stderr:
                raise AssertionError('Poisoned host import did not reproduce: '+mode)
            calls.container('pvac','python','-c',
                            'import numpy,pandas; print(numpy.__version__,numpy.__file__); '
                            'assert numpy.__file__.startswith("/usr/local/"); '
                            'assert pandas.__file__.startswith("/usr/local/")')
            calls.container('pvac','pvacseq','--help')
            results[mode] = {'legacy_exit':legacy.returncode,'fixed_import_and_cli':'PASS'}
    finally:
        os.environ.clear()
        os.environ.update(previous)
    path = root/'proximal.vep.vcf'
    if args.vcf:
        path.write_bytes(args.vcf.read_bytes())
    else:
        header = '##fileformat=VCFv4.2\n##contig=<ID=chr3,length=1000>\n##contig=<ID=chr6,length=1000>\n##contig=<ID=chr10,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        path.write_text(header+''.join(f'{c}\t{p}\t.\tA\tT\t.\tPASS\t.\n' for c,p in
                                     [('chr3',20),('chr6',80),('chr10',30),('chr6',15)]))
    before = path.read_text().splitlines()
    sort_annotation_vcf(path)
    after = path.read_text().splitlines()
    assert [l for l in before if l.startswith('#')] == [l for l in after if l.startswith('#')]
    assert Counter(l for l in before if not l.startswith('#')) == Counter(l for l in after if not l.startswith('#'))
    calls.container('samtools','bgzip','-kf',path)
    calls.container('samtools','tabix','-fp','vcf',str(path)+'.gz')
    results['vcf'] = {'index':'PASS','records':sum(not l.startswith('#') for l in after),'contents_preserved':True}
    (root/'completion.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__ == '__main__':
    main()
