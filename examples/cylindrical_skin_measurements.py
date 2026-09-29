"""One-dimensional cylindrical limit of the production induction flux machinery."""
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.vacuum_flux import flux_components
from .cylindrical_skin_reference import CylindricalSkin
from .conservative_measurements import provenance
from .verification_records import write_record


class RadialSkin:
    """Fixed gas used only to supply declared LLF speed; evolve magnetic flux."""
    def __init__(self, cells, reference=None, method="standard", cfl=.2):
        if method not in ("physical", "standard", "strip"):
            raise ValueError("Unknown skin method")
        self.reference = reference or CylindricalSkin()
        self.method, self.cfl = method, cfl
        self.mesh = MHDMesh(self.reference.radius, 4*self.reference.radius/cells, cells, 8)
        self.psi = self.reference.fields(self.mesh.rf, 0.)[0]
        self.time, self.steps = 0., 0
        self.gamma = 5/3
        self.forced = np.zeros((cells, 8), bool); self.forced[-1] = True

    def fields(self, psi, time):
        mesh, ref = self.mesh, self.reference
        b = np.diff(psi)/(mesh.r*mesh.dr)
        current = np.zeros(mesh.nr+1)
        current[1:-1] = -np.diff(b)/(MU0*mesh.dr)
        surface = ref.field_amplitude*np.cos(ref.omega*time)
        current[-1] = (b[-1]-surface)/(MU0*.5*mesh.dr)
        return b, current, surface

    def rate(self, psi, time):
        mesh, ref = self.mesh, self.reference
        b, current, surface = self.fields(psi, time)
        electric = ref.diffusivity*MU0*current
        if self.method == "physical":
            return -mesh.rf*electric
        flux = np.broadcast_to(psi[:, None], (mesh.nr+1, mesh.nz+1))
        q = np.zeros((mesh.nr, mesh.nz, 6))
        q[..., 0] = q[..., 3] = 1.; q[..., 5] = b[:, None]
        # Exactly the wall value of the production a*h*mu0/4 resistance.
        speed = np.sqrt(self.gamma+b[-1]**2/MU0)
        numerical = np.zeros_like(flux)
        numerical[-1] = MU0*.25*mesh.dr*speed*current[-1]
        physical = np.broadcast_to(electric[:, None], flux.shape)
        boundary = (np.full(mesh.nz, surface), np.zeros(mesh.nr), np.zeros(mesh.nr))
        high = flux_components(mesh, q, flux, self.gamma, physical, boundary, numerical, 2)
        if self.method == "strip":
            low = flux_components(mesh, q, flux, self.gamma, physical, boundary, numerical, 1)
            high, _ = high.local_fallback(low, self.forced)
        # Invariant-z reduction: the middle plane has no wall/end intersection.
        return -mesh.rf*high.electric[:, mesh.nz//2]

    def advance(self, end):
        mesh, ref = self.mesh, self.reference
        speed = np.sqrt(self.gamma+ref.field_amplitude**2/MU0)
        numerical_d = 0. if self.method == "physical" else .25*mesh.dr*speed
        limit = self.cfl/(2*speed/mesh.dr+4*(ref.diffusivity+numerical_d)/mesh.dr**2)
        while end-self.time > 8*np.finfo(float).eps*max(end, 1.):
            dt = min(limit, end-self.time)
            a = self.rate(self.psi, self.time)
            b = self.rate(self.psi+dt*a, self.time+dt)
            c = self.rate(self.psi+.25*dt*(a+b), self.time+.5*dt)
            self.psi += dt*((a+b)/6+2*c/3)
            self.time += dt; self.steps += 1


def measure(cells, method, cfl=.2):
    start = perf_counter(); solver = RadialSkin(cells, method=method, cfl=cfl)
    ref, mesh = solver.reference, solver.mesh
    rows = []
    for time in np.linspace(0., 1/ref.frequency, 9)[1:]:
        solver.advance(float(time))
        exact_psi, _, exact_j = ref.fields(mesh.rf, solver.time)
        exact_b = np.diff(exact_psi)/(mesh.r*mesh.dr)
        b, current, _ = solver.fields(solver.psi, solver.time)
        bscale = ref.field_amplitude
        jscale = ref.field_amplitude/(MU0*ref.skin_depth)
        rows.append({"time_s": solver.time,
                     "B_L2_over_applied": float(np.sqrt(np.sum((b-exact_b)**2*mesh.radial_measure)/mesh.radial_measure.sum())/bscale),
                     "J_L2_over_skin_scale": float(np.sqrt(np.sum((current[1:-1]-exact_j[1:-1])**2*mesh.rf[1:-1]*mesh.dr)/(.5*mesh.radius**2))/jscale),
                     "wall_J_error_over_skin_scale": float(abs(current[-1]-exact_j[-1])/jscale),
                     "flux_error_over_applied_flux": float(abs(solver.psi[-1]-exact_psi[-1])/(.5*mesh.radius**2*bscale)),
                     "Bz_t": b.tolist(), "exact_cell_Bz_t": exact_b.tolist(),
                     "current_a_m2": current.tolist(), "exact_current_a_m2": exact_j.tolist()})
    return {"completed": True, "radial_cells": cells, "method": method, "cfl": cfl,
            "reference": asdict(ref), "gamma": solver.gamma, "fixed_density_kg_m3": 1., "fixed_pressure_pa": 1.,
            "skin_depth_m": ref.skin_depth, "cells_per_skin_depth": ref.skin_depth/mesh.dr,
            "steps": solver.steps, "seconds": perf_counter()-start, "r_m": mesh.r.tolist(), "rf_m": mesh.rf.tolist(),
            "history": rows, "maximum_B_error": max(r["B_L2_over_applied"] for r in rows),
            "maximum_wall_J_error": max(r["wall_J_error_over_skin_scale"] for r in rows),
            "maximum_flux_error": max(r["flux_error_over_applied_flux"] for r in rows)}


def run(directory="outputs/boundary_verification/skin", *, reuse=False):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    for name in ("cylindrical_skin_reference", "cylindrical_skin_measurements"):
        source[name+"_sha256"] = hashlib.sha256(Path(f"examples/{name}.py").read_bytes()).hexdigest()
    rows = []
    cases = [(n, m, .2) for m in ("physical", "standard", "strip") for n in (16, 32, 64, 128)]
    cases += [(64, m, .1) for m in ("physical", "standard", "strip")]
    for n, method, cfl in cases:
        path = directory/f"{method}_{n}_cfl{cfl}.json"
        row = json.loads(path.read_text()) if reuse and path.exists() else None
        if not row or row["provenance"] != source:
            print("Starting skin", n, method, cfl, flush=True)
            row = measure(n, method, cfl); row["provenance"] = source
            write_record(path, row)
        print("Skin", n, method, cfl, "B error", row["maximum_B_error"], flush=True)
        rows.append(row)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true")
    run(reuse=parser.parse_args().reuse)
