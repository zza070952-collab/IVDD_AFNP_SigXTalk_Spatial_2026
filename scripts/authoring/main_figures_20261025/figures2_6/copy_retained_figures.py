#!/usr/bin/env python3
"""Copy old Figures 2–4 to new Figures 3–5 without changing artwork bytes."""
import argparse
import hashlib
import json
from pathlib import Path

from _io_guard import independent_output, reject_output_links


TARGETS = {
    3: ('Figure2', 'Figure2_Frozen_Programs_CrossCohort'),
    4: ('Figure3', 'Figure3_Predicted_Communication'),
    5: ('Figure4', 'Figure4_Anatomy_and_Depth'),
}
FORMATS = ('pdf', 'svg', 'png', 'tiff')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True,
                        help='External read-only V4 project with 04_figures/.')
    parser.add_argument('--output-root', type=Path, required=True,
                        help='Independent output project; writes figures/ and audit/.')
    args = parser.parse_args()
    source, output = independent_output(args.source_root, args.output_root)
    reject_output_links(output)
    jobs = []
    for number, (folder, stem) in TARGETS.items():
        for suffix in FORMATS:
            src = source / '04_figures' / folder / (stem + '.' + suffix)
            if not src.is_file():
                raise FileNotFoundError(src)
            jobs.append((number, suffix, src,
                         output / 'figures' / ('Figure' + str(number) + '.' + suffix), sha(src)))
    (output / 'figures').mkdir(parents=True, exist_ok=True)
    (output / 'audit').mkdir(parents=True, exist_ok=True)
    rows = []
    for number, suffix, src, dest, original in jobs:
        data = src.read_bytes()
        if hashlib.sha256(data).hexdigest() != original:
            raise RuntimeError('Input changed before copy: ' + str(src))
        dest.write_bytes(data)
        if sha(src) != original or sha(dest) != original:
            raise RuntimeError('Byte identity failed: ' + str(src))
        rows.append({'old_figure': number - 1, 'new_figure': number,
                     'format': suffix, 'input': src.relative_to(source).as_posix(),
                     'output': dest.relative_to(output).as_posix(), 'sha256': original,
                     'byte_identical': True, 'panel_content_changed': False})
    if any(sha(src) != original for _, _, src, _, original in jobs):
        raise RuntimeError('An input changed during the copy operation.')
    receipt = output / 'audit' / 'FIGURES3_5_BYTE_IDENTITY.json'
    receipt.write_text(json.dumps(rows, indent=2) + '\n', encoding='utf8')
    print(json.dumps({'unchanged_artwork_formats_copied': len(rows),
                      'old_inputs_unchanged': True}))


if __name__ == '__main__':
    main()
