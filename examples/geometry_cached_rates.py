"""Reuse immutable compatible geometry across same-grid signed workflows."""
import json
from pathlib import Path
from . import direct_metric_rates as direct
from . import regular_plasma_reference as reference
from .implicit_plasma_records import provenance


def run(destination, specs, **options):
    original = reference.CompatibleGeometry, direct.provenance
    cache = {}; hits = 0
    def geometry(mesh, **settings):
        nonlocal hits
        key = (mesh.radius, mesh.half_length, mesh.nr, mesh.nz,
               tuple((name, type(value), value) for name, value in sorted(settings.items())))
        if key not in cache:
            cache[key] = original[0](mesh, **settings)
        else:
            hits += 1
        return cache[key]
    try:
        reference.CompatibleGeometry = geometry
        direct.provenance = lambda seeds: provenance([*seeds, 'examples.geometry_cached_rates',
            'examples.direct_metric_rates', 'examples.regular_plasma_reference'])
        result = direct.run(destination, specs, **options)
        result['geometry_cache'] = dict(constructions=len(cache), hits=hits,
            scope='Only original field geometry within this run; every fluid/state/density and actual operator is fresh')
        (Path(destination)/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
        return result
    finally:
        reference.CompatibleGeometry, direct.provenance = original
