"""Protect external inputs when producing authoring-only outputs."""
from pathlib import Path


def independent_output(source_root, output_root):
    """Require disjoint source/output trees, including resolved symlink targets."""
    source = Path(source_root).resolve()
    output = Path(output_root).resolve()
    if source == output or source in output.parents or output in source.parents:
        raise ValueError('Source and output must be disjoint directory trees.')
    return source, output


def reject_output_links(output_root):
    """Existing output links may point at read-only inputs; refuse them."""
    output = Path(output_root)
    if output.exists():
        for path in output.rglob('*'):
            if path.is_symlink() or getattr(path, 'is_junction', lambda: False)():
                raise ValueError('Output tree contains a link: ' + str(path))
