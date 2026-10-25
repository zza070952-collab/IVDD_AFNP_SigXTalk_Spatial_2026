"""Keep manuscript-bearing runtime artifacts outside this publication checkout."""
from pathlib import Path


def require_private_output(path):
    """Resolve an explicit output and reject paths inside this code checkout.

    The nearest Git root is used so a private workspace may itself have a
    separate parent repository. This check is a guard, not access control.
    """
    target = Path(path).expanduser().resolve()
    here = Path(__file__).resolve()
    root = next((p for p in here.parents if (p / '.git').exists()), None)
    if root is not None and (target == root or root in target.parents):
        raise ValueError('Runtime output may contain manuscript text; choose a private output outside this publication checkout.')
    return target
