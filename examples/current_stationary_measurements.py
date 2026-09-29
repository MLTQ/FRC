"""Repeated uniform-field controls with independently reconstructed physical B."""
import hashlib
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.current_thermal_midpoint import CurrentThermalMidpoint
from .current_thermal_norms import physical_energy


def measure(directory,nr,bias,dt,steps=10):
    mesh=MHDMesh(.035,.12,nr,2*nr);n=np.full((nr+1,2*nr),1e20)
    model=CurrentThermal(mesh,n,39.948*atomic_mass);geo=model.geometry
    solver=CurrentThermalMidpoint(model);psi=.5*bias*mesh.corners[0]**2;g=np.zeros((nr,2*nr))
    ue=model.electron.from_pressure(n*5*E);ui=model.ion.from_pressure(n*E)
    arrays=dict(number_density=n,radial_faces=mesh.rf,axial_faces=mesh.zf);history=[]
    def state(prefix):
        current=model.current(geo.pack(psi,g))
        for name,value in dict(psi=psi,g=g,electron_entropy=ue,ion_entropy=ui,current=current).items():
            arrays[prefix+'_'+name]=value.copy()
        return physical_energy(mesh,psi,g,current,n,ue,ui,reduced_mass=model.reduced_mass,order=10)
    before=state('initial');initial_energy=float(sum(np.sum(mesh.volumes*v) for v in before.values()))
    flux_energy=0.
    for index in range(steps):
        inputs=[a.copy() for a in (psi,g,ue,ui)]
        step=solver.step(psi,g,ue,ui,dt)
        assert all(np.array_equal(a,b) for a,b in zip(inputs,(psi,g,ue,ui)))
        psi,g,ue,ui=step.psi,step.g,step.electron_entropy,step.ion_entropy
        # No production curl/field interpolation is used for the physical error.
        bz=np.diff(psi,axis=0)/mesh.radial_measure[:,None]
        br=-np.diff(psi,axis=1)[1:]/(mesh.rf[1:,None]*mesh.dz)
        maximum=max(float(abs(bz-bias).max()),float(abs(br).max()),float(abs(mesh.rf[1:,None]*g).max()))
        flux_energy+=float(np.sum(mesh.volumes*step.total_energy_change))
        history.append(dict(time_s=(index+1)*dt,magnetic_max_error_t=maximum,
            iterations=step.iterations,nonlinear_residual=step.nonlinear_residual,
            linear_residual=step.linear_residual,stored_state_residual=step.stored_state_residual,
            conservative_energy_change_j=flux_energy))
        assert maximum<2e-13*abs(bias)
    after=state('final');final_energy=float(sum(np.sum(mesh.volumes*v) for v in after.values()))
    for key in ('electron_entropy','ion_entropy'):
        np.testing.assert_array_equal(arrays['initial_'+key],arrays['final_'+key])
    relative=abs(final_energy-initial_energy)/initial_energy
    assert relative<2e-13 and abs(flux_energy)/initial_energy<2e-13
    name=f'stationary_{nr}_{bias:g}_{dt:g}';path=directory/(name+'.npz')
    np.savez_compressed(path,**arrays)
    with np.load(path,allow_pickle=False) as saved:
        for key in arrays:np.testing.assert_array_equal(saved[key],arrays[key])
    return dict(name=name,completed=True,nr=nr,nz=2*nr,bias_t=bias,dt_s=dt,steps=steps,
        radius_m=.035,length_m=.24,density_m3=1e20,ion_mass=model.ion_mass,electron_ev=5.,ion_ev=1.,
        history=history,independent_energy_relative_drift=relative,
        initial_energy_j=initial_energy,final_energy_j=final_energy,
        fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
