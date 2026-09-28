"""Waveform-normalized comparisons and typed parameter sweeps."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
from .config import Experiment
from .waveforms import Unipolar, DampedSine, SineBurst
from .electromagnetics import loop_field
from .simulation import run


def comparison_config(kind, peak_field=.015, *, base=None):
    base = Experiment() if base is None else base
    if not np.isfinite(peak_field) or peak_field <= 0:
        raise ValueError("Peak perturbation must be finite and positive")
    waves = {"unipolar": Unipolar(1), "rlc": DampedSine(1), "burst": SineBurst(1)}
    if kind not in waves:
        raise ValueError(f"Choose {tuple(waves)}")
    wave = waves[kind]
    dense = np.linspace(0, base.duration, 20001)
    v = np.abs(wave.value(dense))
    idx = int(v.argmax())
    left, right = dense[max(idx-1, 0)], dense[min(idx+1, len(dense)-1)]
    optimum = minimize_scalar(lambda t: -abs(float(wave.value(t))), bounds=(left, right),
                              method="bounded", options={"xatol": 1e-16})
    peak = max(float(v.max()), -optimum.fun)
    field_per_amp = sum(c.turns*c.polarity*float(loop_field(0, 0, c.radius, c.z)[2]) for c in base.coils)
    if abs(field_per_amp) < 1e-20 or peak <= 0:
        raise ValueError("Cannot normalize zero center-field response")
    wave = replace(wave, amplitude=peak_field/(abs(field_per_amp)*peak))
    return replace(base, coils=tuple(replace(c, waveform=wave) for c in base.coils), label=kind)


def compare(directory, *, base=None, peak_field=.015):
    directory = Path(directory)
    results = {}
    for name in ("unipolar", "rlc", "burst"):
        config = comparison_config(name, peak_field, base=base)
        result = run(config)
        result.metadata["normalization"] = {"quantity": "peak absolute external perturbation Bz at r=z=0", "value_t": peak_field,
                                             "not_equalized": ["source energy", "Ohmic heat", "positive field excursion"]}
        result.save(directory/name)
        results[name] = result
    return results


def sweep(directory, *, base=None, conductivities=(5e3, 2e4, 8e4), peak_ratios=(1., 3., 6.)):
    base = Experiment() if base is None else base
    if base.bias == 0:
        raise ValueError("Peak/bias normalization requires nonzero bias")
    rows = []
    for kind in ("unipolar", "rlc", "burst"):
        for sigma in conductivities:
            for ratio in peak_ratios:
                config = comparison_config(kind, abs(base.bias)*ratio,
                                           base=replace(base, material=replace(base.material, conductivity=float(sigma))))
                result = run(config)
                rows.append({"waveform": kind, "conductivity_s_m": sigma, "peak_external_over_bias": ratio,
                             "summary": result.metadata["summary"], "configuration": config.to_dict(),
                             "warnings": result.metadata["warnings"]})
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    (directory/"sweep.json").write_text(json.dumps({"normalization": "equal peak absolute vacuum center perturbation within each amplitude column",
                                                    "provenance": {k: v for k, v in result.metadata.items() if k in ("source_sha256", "dependencies", "python", "platform", "method")},
                                                    "rows": rows}, indent=2, allow_nan=False)+"\n")
    return rows
