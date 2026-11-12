#!/usr/bin/env python3
"""Diagnostic observer after native Google allele-counter normalization.

Run under the existing Google 1.10 image, under Slurm, with explicit bounded
native --regions. DS_CONTEXT_EXPORT_DIR names a new scratch directory. The
packaged Python runtime is unpacked there; nothing is installed or replaced.
All scientific options and computation are the native make_examples_somatic
entry. The observer writes existing in-memory reads after candidates_in_region;
it does not alter them, synthesize bases/CIGARs, or change returned candidates.
The exported BAMs are context diagnostics, not replacement counting inputs.
"""
import os
from pathlib import Path
import stat
import sys
import zipfile


def main():
    if not any(x == '--regions' or x.startswith('--regions=') for x in sys.argv[1:]):
        raise SystemExit('Explicit bounded --regions required for diagnostic export')
    out = Path(os.environ['DS_CONTEXT_EXPORT_DIR'])
    out.mkdir(exist_ok=False)
    runtime = out / 'runtime'
    with zipfile.ZipFile('/opt/deepvariant/bin/make_examples_somatic.zip') as z:
        z.extractall(runtime)
        # Match the native Bazel launcher: extension aliases are ZIP symlinks,
        # not ELF files containing the text of a relative target.
        for member in z.infolist():
            if stat.S_ISLNK(member.external_attr >> 16):
                path = runtime / member.filename
                path.unlink()
                path.symlink_to(z.read(member).decode())
    entry = list(runtime.rglob('deepvariant/make_examples_somatic.py'))
    if len(entry) != 1:
        raise RuntimeError('Expected one packaged native entry')
    root = entry[0].parents[1]
    # Match this archive's Bazel launcher, including its bundled C++ protobuf.
    # Appending these after site-packages can incorrectly select system upb.
    paths = [root.parent, root.parent / 'com_google_protobuf/python', root]
    paths.extend(p for p in root.parent.iterdir() if p.is_dir())
    sys.path[:0] = list(dict.fromkeys(map(str, paths)))
    from absl import app
    from deepvariant import make_examples_core, make_examples_somatic
    from third_party.nucleus.io import sam

    original = make_examples_core.RegionProcessor.candidates_in_region

    def observe(self, region, *args, **kwargs):
        result = original(self, region, *args, **kwargs)
        for sample in self.samples:
            if not sample.options.reads_filenames:
                continue
            path = out / (f'{sample.options.role}.{region.reference_name}.'
                          f'{region.start}-{region.end}.normalized.bam')
            with sam.SamWriter(str(path), header=sample.sam_readers[0].header) as writer:
                for read in sample.in_memory_sam_reader.query(region):
                    writer.write(read)
        return result

    make_examples_core.RegionProcessor.candidates_in_region = observe
    app.run(make_examples_somatic.main)


if __name__ == '__main__':
    main()
