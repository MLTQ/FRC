"""Kinetic heat-wave evolutions and independent local/ballistic comparisons."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from plasma_bench.kinetic_heat import KineticHeat
from .implicit_plasma_records import provenance,freeze


def characteristic_reference(duration,intervals):
    u=np.linspace(0,duration,intervals+1);du=u[1]
    gaussian=np.exp(-u*u/2);kernel=(u**3-3*u)*gaussian
    rhs=(1-2.5*u*u+.5*u**4)*gaussian;field=np.empty_like(u);field[0]=1.
    heat_kernel=(5*u-u**3)*gaussian/3
    temperature=(1-7*u*u/6+u**4/6)*gaussian
    for i in range(1,len(u)):
        field[i]=rhs[i]-du*(.5*field[0]*kernel[i]+field[1:i]@kernel[i-1:0:-1])
        temperature[i]+=du*(.5*field[0]*heat_kernel[i]+field[1:i]@heat_kernel[i-1:0:-1])
    return u,temperature,field


def run(destination, *, degrees=(8,12,16), reference_intervals=4000):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh kinetic heat output required')
    if not degrees or any(d<3 for d in degrees) or reference_intervals<400 or reference_intervals%400:
        raise ValueError('Degrees >=3 and reference intervals divisible by400 required')
    out.mkdir(parents=True);identity=provenance(['examples.kinetic_heat_relaxation']);freeze(out,identity)
    started=time.monotonic();report=dict(terminal=False,source_sha256=identity['source_sha256'],
        scope='Linear uniform periodic electron heat transport, stationary infinite-mass ions, zero current. No vacuum edge, finite drift, nonlinear heat fronts or full formation coupling.',
        normalization='s=v/sqrt(T/m), tau=nu0*t, K=k*vT/nu0; nu0=n*4*pi*(e²/(4*pi*eps0))²*lnLambda/(m²*vT³)',
        degrees=list(degrees),reference_intervals=reference_intervals,local=[],waves=[],ballistic=[])
    def save():
        report['elapsed_seconds']=time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();braginskii=3.16*3*np.sqrt(np.pi/2);lorentz=64*np.sqrt(2/np.pi)
    report['independent_local_coefficients']={'approximate_braginskii_z1':braginskii,'analytic_lorentz':lorentz}
    fig,axes=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for degree in degrees:
        model=KineticHeat(degree);local=model.local_conductivity();lorentz_numeric=model.local_conductivity(lorentz_only=True)
        np.savez_compressed(out/f'operators_degree{degree}.npz',streaming=model.streaming,
                            electron_collision=model.self_collision,ion_collision=model.ion_collision,
                            indices=np.array(model.indices),active=model.active,heat=model.heat)
        report['local'].append(dict(degree=degree,moments=model.count,conductivity=local,
            conductivity_in_braginskii_units=local/(3*np.sqrt(np.pi/2)),
            lorentz_conductivity=lorentz_numeric,lorentz_relative_difference=lorentz_numeric/lorentz-1))
        print('degree',degree,'local',report['local'][-1],flush=True)
        for axis,k in zip(axes.flat,(.001,.03,.3,1.)):
            duration=5/((2/3)*braginskii*k*k) if k==.001 else 4/k
            times=np.linspace(0,duration,401);result=model.evolve(k,times)
            filename=f'wave_degree{degree}_K{k:g}.npz';np.savez_compressed(out/filename,**result)
            classical=np.exp(-(2/3)*braginskii*k*k*times)
            report['waves'].append(dict(degree=degree,K=k,filename=filename,duration_nu0_t=duration,
                final_temperature_amplitude=float(result['temperature'][-1].real),
                maximum_difference_from_approximate_classical=float(abs(result['temperature']-classical).max()),
                maximum_imaginary_temperature=float(abs(result['temperature'].imag).max()),
                maximum_positive_free_energy_increment=float(np.diff(result['quadratic_free_energy']).max()),
                minimum_collision_dissipation=float(result['collision_dissipation'].min()),
                maximum_eigenvalue_real=result['maximum_eigenvalue_real'],eigenvector_condition=result['eigenvector_condition']))
            axis.plot(times*k,result['temperature'].real,label=f'degree {degree}')
            if degree==degrees[-1]:axis.plot(times*k,classical,'k--',label='Local Braginskii approximation')
            axis.set(title=f'K = {k:g}',xlabel='k vT t',ylabel='Temperature amplitude / initial');axis.grid(alpha=.2);axis.legend(fontsize=8)
        save()
    # Independent collisionless characteristic integrals, not the moment matrix.
    u,reference,field=characteristic_reference(4.,reference_intervals)
    uc,coarse,ec=characteristic_reference(4.,reference_intervals//2)
    report['ballistic_reference_maximum_change_on_halving_interval']=float(abs(reference[::2]-coarse).max())
    np.savez_compressed(out/'characteristic_reference.npz',u=u,temperature=reference,field=field,
                        coarse_u=uc,coarse_temperature=coarse)
    for degree in (12,20,28):
        stride=reference_intervals//400
        model=KineticHeat(degree,collisions=False);result=model.evolve(1.,u[::stride])
        np.savez_compressed(out/f'ballistic_degree{degree}.npz',**result)
        error=float(abs(result['temperature']-reference[::stride]).max())
        report['ballistic'].append(dict(degree=degree,moments=model.count,maximum_temperature_error=error,
            maximum_quadratic_free_energy_change=float(abs(result['quadratic_free_energy']/result['quadratic_free_energy'][0]-1).max())))
        axes[1,1].plot(result['time'],result['temperature'].real,label=f'degree {degree}')
        print('ballistic',report['ballistic'][-1],flush=True)
    axes[1,1].plot(u,reference,'k--',label='Independent characteristics')
    axes[1,1].set(title='Collisionless zero-current heat relaxation',xlabel='k vT t',ylabel='Temperature amplitude / initial')
    axes[1,1].legend(fontsize=8);axes[1,1].grid(alpha=.2)
    degrees=[r['degree'] for r in report['local']]
    axes[1,2].plot(degrees,[r['conductivity_in_braginskii_units'] for r in report['local']],'o-',label='Full linearized Landau, Z=1')
    axes[1,2].axhline(3.16,color='black',ls='--',label='Braginskii approximation')
    axes[1,2].set(title='Local conductivity from kinetic equation',xlabel='Maximum velocity-polynomial degree',ylabel='Conductivity / (n vT² tau_e)')
    axes[1,2].legend(fontsize=8);axes[1,2].grid(alpha=.2)
    report['terminal']=True;save()
    fig.suptitle('Linear kinetic heat transport\nFull electron thermalization and velocity-dependent ion scattering; periodic fixed background, no pulse coupling')
    fig.savefig(out/'heat_relaxation.png',dpi=170);fig.savefig(out/'heat_relaxation.pdf');plt.close(fig)
    print('Kinetic heat calculation complete',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination')
    parser.add_argument('--degrees',type=int,nargs='+',default=[8,12,16])
    parser.add_argument('--reference-intervals',type=int,default=4000)
    args=parser.parse_args();run(args.destination,degrees=tuple(args.degrees),reference_intervals=args.reference_intervals)
