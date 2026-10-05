"""Independent physical audit of terminal, incomplete contact checkpoints."""
import json
from pathlib import Path
import numpy as np
from .regular_plasma_reference import initialize
from .implicit_plasma_audit import measure
from .implicit_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    seeds = ['examples.implicit_contact_diagnosis']; identity = provenance(seeds)
    destination.mkdir(parents=True); freeze(destination, identity)
    entries = []
    for name, elements, steps in [('implicit_gauss_newton_contact', 4, 120),
                                 ('implicit_gauss_contact_time240', 4, 240),
                                 ('implicit_gauss_contact_ne8', 8, 120)]:
        directory = Path('outputs')/name
        terminal = json.loads((directory/'termination.json').read_text())
        if not terminal['terminal']: raise ValueError('Terminal checkpoint required')
        record = terminal['saved_checkpoint']; raw = directory/'progress.npz'
        if sha(raw) != record['arrays_sha256']: raise ValueError('Original saved arrays required')
        spec = dict(kind='contact', elements=elements, steps=steps, method='gauss')
        _, _, u, old, field = initialize(spec)
        saved = dict(np.load(raw)); arrays = dict(initial_state=u, final_state=saved['state'],
            initial_magnetic=field.packed_increment(old.geometry), final_magnetic=saved['magnetic'])
        entries.append(dict(input_directory=str(directory), arrays_sha256=sha(raw),
            input_source_sha256=record['source_sha256'], termination=terminal,
            complete_transit=False, order16=measure(spec, arrays, record['time'], 16),
            order20=measure(spec, arrays, record['time'], 20)))
    if provenance(seeds) != identity: raise RuntimeError('Audit source changed')
    (destination/'summary.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],
        terminal=True, cases=entries), indent=2)+'\n')
    print('incomplete contact checkpoints independently audited', flush=True)


if __name__ == '__main__': run('outputs/implicit_contact_failure_audit')
