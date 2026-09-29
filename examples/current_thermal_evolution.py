"""Actual finite-mass steps with separate conservative and component energy."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.current_thermal_midpoint import CurrentThermalMidpoint
from plasma_bench.hall_midpoint import HallConvergenceError
from .current_beltrami_reference import CurrentBeltrami
from .electron_reference import ElectronReference
from .current_thermal_norms import wave_norms,wave_face_power


def initialize(case):
    ref=(CurrentBeltrami if case['kind']=='wave' else ElectronReference)(**case['reference'])
    mesh=MHDMesh(ref.radius,ref.length/2,case['nr'],case['nz'])
    if case['kind']=='wave':
        n=np.full((mesh.nr+1,mesh.nz),ref.number_density)
        pe,pi=n*E*ref.electron_ev,n*E*ref.ion_ev
        psi,g=ref.initial_fields(mesh);end=case['phase']/abs(ref.omega)
    else:
        f=ref.fields(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'))
        n,pe=f['density'],f['pressure'];pi=np.full_like(pe,case['ion_pressure'])
        psi=np.zeros((mesh.nr+1,mesh.nz+1));g=ref.cell_average(mesh,'g');end=case['end_time']
    model=CurrentThermal(mesh,n,case['ion_mass'],model=case['model'],gamma_e=case['gamma_e'],
                         gamma_i=case['gamma_i'],density_scale=case['density_scale'],quadrature_order=case['quadrature'])
    return ref,model,psi,g,model.electron.from_pressure(pe),model.ion.from_pressure(pi),end


def measure(case):
    ref,model,psi,g,ue,ui,end=initialize(case);mesh=model.mesh;geo=model.geometry
    solver=CurrentThermalMidpoint(model,tolerance=case['tolerance'],max_iterations=case['max_iterations'])
    total=sum(model.energy_density(geo.pack(psi,g),ue,ui).values())
    arrays=dict(number_density=model.electron.density,radial_faces=mesh.rf,axial_faces=mesh.zf,
                primal_volumes=mesh.volumes,dual_volumes=model.electron.geometry.volumes)
    def state(prefix,p,t,e,i,et):
        arrays.update({prefix+'_'+k:v.copy() for k,v in dict(psi=p,g=t,electron_entropy=e,ion_entropy=i,
            total_energy=et,current=model.current(geo.pack(p,t))).items()})
    def diagnostics(p,t,e,i,et,time,iterations=0,residual=0.,linear=0.,stored=0.):
        packed=geo.pack(p,t);d=model.diagnostics(packed,e,i);energies=model.energy_density(packed,e,i)
        defect=(case['gamma_i']-1)*(et-sum(energies.values()))
        inferred=(case['gamma_i']-1)*(et-energies['magnetic']-energies['current_kinetic']-energies['electron'])
        d.update(time_s=time,conservative_energy_j=float(np.sum(mesh.volumes*et)),
            min_inferred_ion_pressure_pa=float(inferred.min()),
            local_pressure_defect_rms_pa=float(np.sqrt(np.sum(mesh.volumes*defect*defect)/np.sum(mesh.volumes))),
            local_pressure_defect_max_pa=float(abs(defect).max()),
            axis_pressure_defect_max_pa=float(abs(defect[0]).max()),
            nonlinear_iterations=iterations,nonlinear_residual=residual,linear_residual=linear,stored_state_residual=stored)
        return d
    initial=diagnostics(psi,g,ue,ui,total,0.);history={k:[v] for k,v in initial.items()};last=initial
    state('initial',psi,g,ue,ui,total);status='succeeded';failure=None;accepted=0;dt=end/case['steps']
    for step in range(case['steps']):
        unchanged=[v.copy() for v in (psi,g,ue,ui)]
        try:
            trial=solver.step(psi,g,ue,ui,dt)
        except HallConvergenceError as error:
            assert all(np.array_equal(a,b) for a,b in zip(unchanged,(psi,g,ue,ui)))
            status='rejected';failure=dict(reason=str(error),trial_time_s=(step+1)*dt,input_unchanged=True)
            break
        next_total=total+trial.total_energy_change
        diagnostic=diagnostics(trial.psi,trial.g,trial.electron_entropy,trial.ion_entropy,next_total,
            (step+1)*dt,trial.iterations,trial.nonlinear_residual,trial.linear_residual,trial.stored_state_residual)
        if diagnostic['min_inferred_ion_pressure_pa'] <= 0:
            state('rejected',trial.psi,trial.g,trial.electron_entropy,trial.ion_entropy,next_total)
            status='rejected';failure=dict(reason='Conservative ion partition became nonpositive; no repair',
                trial_time_s=(step+1)*dt,diagnostic=diagnostic,input_unchanged=True)
            break
        arrays.update(last_mid_psi=.5*(psi+trial.psi),last_mid_g=.5*(g+trial.g),
            last_electric=trial.electric_coefficients,last_radial_power=trial.radial_energy_flux,last_axial_power=trial.axial_energy_flux)
        psi,g,ue,ui,total=trial.psi,trial.g,trial.electron_entropy,trial.ion_entropy,next_total
        last=diagnostic;accepted+=1
        for k,v in last.items():history[k].append(v)
    state('final',psi,g,ue,ui,total)
    metrics={k:last[k] for k in ('local_pressure_defect_rms_pa','local_pressure_defect_max_pa',
        'axis_pressure_defect_max_pa','min_inferred_ion_pressure_pa')}
    if case['kind']=='wave':
        metrics.update(wave_norms(mesh,ref,last['time_s'],psi,g,arrays['final_current']))
        if accepted:
            exact=wave_face_power(mesh,ref,last['time_s']-.5*dt)
            for label,actual,target in zip(('radial','axial'),(arrays['last_radial_power'],arrays['last_axial_power']),exact):
                metrics[label+'_power_relative_L2']=float(np.linalg.norm(actual-target)/np.linalg.norm(target))
                metrics[label+'_power_rms_error_w']=float(np.sqrt(np.mean((actual-target)**2)))
    return dict(status=status,completed=status=='succeeded',initial=initial,final=last,metrics=metrics,
        accepted_steps=accepted,requested_end_time_s=end,history_columns=history,failure=failure,
        factorizations=solver.factorizations),arrays
