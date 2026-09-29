"""Retain and independently account for the original curl-kernel rejection."""
from pathlib import Path
from dataclasses import asdict
import hashlib
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_mass_state import evaluate
from plasma_bench.smooth_finite_mass_flux import bulk_fluxes
from plasma_bench.primal_velocity import regular_velocity
from plasma_bench.moving_current import MovingCurrent
from plasma_bench.kernel_induction import KernelInduction,separate,difference_curl
from .finite_mass_contact_reference import FiniteMassContact
from .finite_mass_contact_case import initialize
from .moving_current_records import provenance as original_provenance
from .smooth_finite_mass_records import provenance
from .verification_records import write_record


def measure(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    ref=FiniteMassContact();mesh=MHDMesh(ref.radius,ref.length/2,8,16);state,par=initialize(mesh,ref)
    ev=evaluate(mesh,state,par);model=ev.model
    fr,fz,_=bulk_fluxes(mesh,state,ev,par)
    nt=model.mapping.density(-(np.diff(fr[...,0],axis=0)+np.diff(fz[...,0],axis=1))/mesh.volumes)/par.masses.total
    regular=regular_velocity(model.mapping,state.conserved)
    original=MovingCurrent(model,tolerance=par.linear_tolerance,max_iterations=par.max_iterations)
    terms,potential=original.terms(state.psi,state.g,state.electron_entropy,ev.ion_entropy,regular,nt)
    forcing=sum(terms.values());saved=forcing.copy()
    failure=None
    try:original.solver.solve(forcing)
    except RuntimeError as error:failure=str(error)
    assert failure is not None,'The exact historical failure must remain reproducible'
    np.testing.assert_array_equal(forcing,saved)
    solver=KernelInduction(model,tolerance=par.linear_tolerance,max_iterations=par.max_iterations)
    solved=solver.solve(forcing);strong=solver.mass_inverse(forcing);small,kernel=separate(model.geometry,strong)
    np.testing.assert_array_equal(difference_curl(model.geometry,kernel),0.)
    electric=solved.electric;curl=difference_curl(model.geometry,electric)
    ohm=model.mass@electric-model.inertia@solved.current-forcing
    ampere=model.mass@solved.current+model.geometry.faraday.T@(model.hessian@solved.magnetic)
    relative=lambda a,b:float(np.linalg.norm(a)/max(np.linalg.norm(b),np.finfo(float).tiny))
    replay=dict(ohm_backward=relative(ohm,abs(model.mass@electric)+abs(model.inertia@solved.current)+abs(forcing)),
        ampere_backward=relative(ampere,abs(model.mass@solved.current)+abs(model.geometry.faraday.T@(model.hessian@solved.magnetic))),
        faraday_backward=relative(curl-solved.magnetic,abs(model.geometry.faraday)@abs(electric)+abs(solved.magnetic)),
        full_faraday_rate_relative=relative(curl-solved.magnetic,solved.magnetic),
        electric_norm=float(np.linalg.norm(strong)),remainder_norm=float(np.linalg.norm(small)),
        magnetic_rate_energy_norm=float(np.sqrt(solved.magnetic@(model.hessian@solved.magnetic))))
    for key in ('ohm_backward','ampere_backward','faraday_backward'):assert replay[key]<par.linear_tolerance
    path=directory/'contact_kernel.npz'
    np.savez_compressed(path,conserved=state.conserved,psi=state.psi,g=state.g,electron_entropy=state.electron_entropy,
        density_rate=nt,regular_velocity=regular,forcing=forcing,potential=potential,
        strong_electric=strong,kernel=kernel,remainder=small,magnetic_rate=solved.magnetic,
        current_rate=solved.current,electric=solved.electric)
    record=dict(reference=asdict(ref),nr=mesh.nr,nz=mesh.nz,parameters=vars(par),
        original_source=original_provenance(),current_source=provenance(),original_failure=failure,
        new_residual=solved.residual,new_residual_components=solved.residual_components,replay=replay,
        fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()),
        scope='Exact separable Faraday-kernel representation, not tolerance relaxation. Full-array recombination remains rate-relative ill-conditioned; backward errors and weak-signal checks are separate.')
    write_record(directory/'contact_kernel.json',record);return record
