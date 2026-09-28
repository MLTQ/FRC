"""Orchestrate sources, independent solvers, diagnostics and reproducible results."""
import numpy as np
from .config import Experiment
from .electromagnetics import VacuumField
from .plasma import StationaryDiffusion, ResistiveMHD
from .numerics import field_from_flux, cell_divergence
from .diagnostics import characteristic_scales, magnetic_energy, longest_sampled_interval
from .topology import analyze_flux
from .results import Result, provenance


def run(config: Experiment, *, topology=True):
    if config.stage == 2:
        ResistiveMHD()
    metadata = provenance(config)  # Reject non-reproducible callbacks before running.
    g, m = config.grid, config.material
    r, z = g.mesh
    active = m.region.contains(r, z)
    axis_active = m.region.contains(0, g.z)
    sigma = np.where(active, m.conductivity, 0)
    sources = VacuumField(g, config.coils, config.bias)
    times = config.times
    solver = None if config.stage == 0 else StationaryDiffusion(g, sigma, sources.at(times[0]).psi)
    frame_indices = np.unique(np.r_[np.arange(0, len(times), config.output_stride), len(times)-1])
    frames = {k: [] for k in ("psi", "br", "bz", "j_theta", "e_theta", "dpsi_dt")}
    keys = ("min_axis_bz", "min_axis_bz_over_bias", "min_bz", "min_bz_over_bias", "vacuum_min_bz_over_bias",
            "center_bz", "vacuum_center_bz", "edge_bz", "end_bz", "peak_j_theta",
            "chamber_magnetic_energy", "induced_magnetic_energy", "ohmic_heat",
            "source_work", "numerical_loss", "energy_residual")
    history = {key: [] for key in keys}
    topologies = []
    heat = work = loss = residual = absolute_work = 0.0
    center = (0, int(np.argmin(abs(g.z))))
    edge = (int(np.argmin(abs(g.r-m.region.radius*.75))), center[1])
    end = (0, int(np.argmin(abs(g.z-m.region.half_length*.75))))
    maximum_divergence = 0.0
    for i, t in enumerate(times):
        external = sources.at(t)
        psi, br, bz = external.psi, external.br, external.bz
        electric, rate = external.e_theta, external.dpsi_dt
        current = np.zeros(g.shape)
        induced_energy = 0.0
        if solver is not None:
            if i:
                state = solver.step(external.psi, float(t-times[i-1]))
                ib_r, ib_z = field_from_flux(g, state.induced_psi)
                psi, br, bz = external.psi+state.induced_psi, br+ib_r, bz+ib_z
                electric, current = state.e_theta, state.j_theta
                rate = -r*electric
                induced_energy = state.induced_energy
                heat += state.ohmic_heat
                work += state.source_work
                absolute_work += abs(state.source_work)
                loss += state.numerical_loss
                residual += state.balance_residual
            else:
                # Initial state is field-penetrated with no current; interval E is undefined.
                electric, rate = np.zeros(g.shape), np.zeros(g.shape)
        values = dict(min_axis_bz=float(bz[0, axis_active].min()), min_bz=float(bz[active].min()),
                      min_axis_bz_over_bias=float(np.min(bz[0, axis_active]/config.bias)) if config.bias else np.nan,
                      min_bz_over_bias=float(np.min(bz[active]/config.bias)) if config.bias else np.nan,
                      vacuum_min_bz_over_bias=float(np.min(external.bz[active]/config.bias)) if config.bias else np.nan,
                      center_bz=float(bz[center]), vacuum_center_bz=float(external.bz[center]),
                      edge_bz=float(bz[edge]), end_bz=float(bz[end]), peak_j_theta=float(np.abs(current).max()),
                      chamber_magnetic_energy=magnetic_energy(g, br, bz, config.chamber),
                      induced_magnetic_energy=induced_energy, ohmic_heat=heat, source_work=work,
                      numerical_loss=loss, energy_residual=residual)
        for key in keys:
            history[key].append(values[key])
        if i in frame_indices:
            for key, arr in zip(frames, (psi, br, bz, current, electric, rate)):
                frames[key].append(arr.copy())
            topologies.append(analyze_flux(g, psi, m.region) if topology else {"closed_flux": False, "skipped": True})
            maximum_divergence = max(maximum_divergence, float(np.max(np.abs(cell_divergence(g, psi)))))
    arrays = {**{key: np.array(value) for key, value in frames.items()},
              **{key: np.array(value) for key, value in history.items()},
              "times": times, "frame_times": times[frame_indices], "frame_indices": frame_indices,
              "r": g.r, "z": g.z}
    scales, warnings = characteristic_scales(config)
    initial_internal = 1.5*scales["pressure_pa"]*np.pi*m.region.radius**2*2*m.region.half_length
    if heat > .1*initial_internal:
        warnings.append("Accumulated Ohmic heat exceeds 10% of initial reference thermal energy; fixed temperature/conductivity need review.")
    lifetime = longest_sampled_interval(times[frame_indices], [x["closed_flux"] for x in topologies])
    closed = any(x["closed_flux"] for x in topologies) if topology else None
    reversed_field = bool(np.min(arrays["min_bz_over_bias"]) < -1e-6) if config.bias else None
    vacuum_reversed = bool(np.min(arrays["vacuum_min_bz_over_bias"]) < -1e-6) if config.bias else None
    minimum_ratio = float(np.min(arrays["min_bz_over_bias"])) if config.bias else None
    classification = "no reversal"
    if reversed_field:
        classification = "transient reversal"
    if closed:
        classification = "closed flux"
        if lifetime+1e-12*config.duration >= config.persistence_threshold:
            classification = "persistent closed flux (sampled)"
    if config.bias == 0 and not closed:
        classification = "reversal undefined (zero bias)"
    if not topology:
        classification += "; topology not evaluated"
    metadata.update(scales=scales, warnings=warnings, topology=topologies,
                    probes={name: {"r_m": float(g.r[index[0]]), "z_m": float(g.z[index[1]])}
                            for name, index in (("center", center), ("edge", edge), ("end", end))},
                    summary={"classification": classification, "field_reversal": reversed_field,
                             "axis_field_reversal": bool(np.min(arrays["min_axis_bz_over_bias"]) < -1e-6) if config.bias else None,
                             "vacuum_field_reversal": vacuum_reversed, "closed_flux": closed,
                             "min_bz_over_bias": minimum_ratio,
                             "minimum_axis_bz_t": float(np.min(arrays["min_axis_bz"])),
                             "peak_j_theta_a_m2": float(np.max(arrays["peak_j_theta"])),
                             "ohmic_heat_j": heat, "numerical_loss_j": loss,
                             "source_work_j": work,
                             "closed_flux_lifetime_lower_bound_s": lifetime,
                             "energy_balance_relative_error": abs(residual)/max(absolute_work, 1e-30),
                             "max_face_divergence_t_per_m": maximum_divergence,
                             "ohmic_heat_over_initial_thermal_energy": heat/initial_internal,
                             "equilibrium": "not evaluated", "frc_formation": "not established",
                             "stability": "not evaluated", "particle_confinement": "not evaluated",
                             "excluded_flux_radius_m": None})
    return Result(arrays, metadata)
