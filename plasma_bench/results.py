"""Numeric archives plus a complete, inspectable reproducibility manifest."""
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import platform
import numpy as np
from . import __version__


def provenance(config):
    digest = hashlib.sha256()
    root = Path(__file__).parent
    for p in sorted(root.glob("*.py")):
        digest.update(p.name.encode())
        digest.update(p.read_bytes())
    return {"schema_version": 1, "package_version": __version__,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "source_sha256": digest.hexdigest(), "python": platform.python_version(),
            "platform": platform.platform(),
            "dependencies": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib", "contourpy")},
            "configuration": config.to_dict(),
            "units": {"psi": "Wb (physical flux / 2pi)", "br": "T", "bz": "T", "j_theta": "A/m^2", "e_theta": "V/m", "dpsi_dt": "Wb/s", "time": "s", "energy": "J", "length": "m"},
            "method": "analytic elliptic-loop vacuum" if config.stage == 0 else "stationary scalar-conductivity flux diffusion; backward Euler; coupled elliptic vacuum",
            "boundary": "analytic free-space external field" if config.stage == 0 else "analytic free-space external field; induced flux zero on axis and finite outer boundary",
            "initial_state": "zero induced current; initial external field fully penetrates plasma",
            "energy_budget": "induced-field work = stored energy + Ohmic heat + backward-Euler loss; total chamber magnetic energy is separate",
            "topology_scope": "prescribed plasma cylinder; evaluated at saved snapshots",
            "current_time_convention": "instantaneous analytic E; J=0" if config.stage == 0 else "zero initial current; subsequent E and J average over previous integration interval"}


@dataclass
class Result:
    arrays: dict[str, np.ndarray]
    metadata: dict

    def save(self, directory):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(directory/"fields.npz", **self.arrays)
        (directory/"metadata.json").write_text(json.dumps(self.metadata, indent=2, allow_nan=False)+"\n")
        return directory

    @classmethod
    def load(cls, directory):
        directory = Path(directory)
        with np.load(directory/"fields.npz", allow_pickle=False) as data:
            arrays = {key: data[key] for key in data.files}
        return cls(arrays, json.loads((directory/"metadata.json").read_text()))
