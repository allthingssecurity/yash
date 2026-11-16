"""Dispatch public lead commands to the recorded Design #3 entry points."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile

from .io import child_env, require_assertions


def verify(repo):
    """Check recorded bytes without optimization, inference or network access."""
    require_assertions('Lead verification')
    result = subprocess.run(
        [sys.executable, str(Path(repo) / 'designs/design-3/verify_results.py')],
        cwd=repo, capture_output=True, text=True, env=child_env(),
    )
    if result.returncode:
        raise ValueError('Design #3 verification failed:\n' + result.stdout + result.stderr)
    return json.loads(result.stdout)


def reproduce(repo, out=None):
    """Replay a fixed design; independent designs belong in their own workspace."""
    require_assertions('Lead reproduction')
    repo = Path(repo).resolve()
    if out is None:
        with tempfile.TemporaryDirectory(prefix='yash-design3-replay-') as temp:
            result = reproduce(repo, Path(temp) / 'replay')
            result['output'] = None
            return result
    out = Path(out).resolve()
    if out.exists():
        raise ValueError('Choose a fresh output directory: ' + str(out))
    result = subprocess.run(
        [sys.executable, str(repo / 'designs/design-3/replay.py'), '--output', str(out)],
        cwd=repo, capture_output=True, text=True, env=child_env(),
    )
    if result.returncode:
        raise ValueError('Design #3 replay failed:\n' + result.stdout[-3000:] + result.stderr[-3000:])
    receipt = json.loads((out / 'replay.json').read_text())
    return dict(
        receipt, name='design-3', selected_regions=receipt['regions'],
        protein_length_aa=receipt['protein_length_aa'],
        all_solver_phases_optimal=receipt['all_solver_phases_optimal'],
        output=str(out), rna_length_nt=None,
    )
