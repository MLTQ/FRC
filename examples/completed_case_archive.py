"""Archive a fully completed signed case without completing its parent batch."""
import json
from pathlib import Path
import shutil
from .implicit_plasma_records import provenance, freeze, sha


def run(source, index, destination):
    source = Path(source); destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing overwrite')
    parent = json.loads((source/'summary.json').read_text())
    case = json.loads((source/f'case_{index:02}.json').read_text())
    if case['status'] != 'complete' or case != parent['cases'][index]:
        raise ValueError('Independently saved complete case and matching batch record required')
    raw = source/case['arrays_file']
    if sha(raw) != case['arrays_sha256']:
        raise ValueError('Original complete raw arrays required')
    original = json.loads((source/'provenance.json').read_text())
    if parent['source_sha256'] != original['source_sha256'] or any(sha(p) != h for p, h in original['source_files'].items()):
        raise ValueError('Unchanged original calculation source required')
    seeds = ['examples.completed_case_archive', *[p[:-3].replace('/', '.') for p in original['source_files']]]
    identity = provenance(seeds)
    destination.mkdir(parents=True); freeze(destination, identity)
    shutil.copy2(raw, destination/raw.name)
    (destination/'case_00.json').write_text(json.dumps(case, indent=2)+'\n')
    result = dict(terminal=True, all_cases_complete=True, selection_only=True,
        source_sha256=identity['source_sha256'], cases=[case],
        parent_directory=str(source), parent_index=index, parent_batch_terminal=parent['terminal'],
        parent_summary_sha256=sha(source/'summary.json'), parent_source_sha256=original['source_sha256'],
        parent_provenance_sha256=sha(source/'provenance.json'),
        scope='One fully completed signed workflow selected unchanged; no claim that parent batch completed')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
