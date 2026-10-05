"""Measured field-direction and full signed-trajectory verification figures."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from .implicit_plasma_records import provenance, freeze, sha


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    root=Path('outputs');inputs={}
    def read(name):
        path=root/name/'summary.json';inputs[str(path)]=sha(path)
        value=json.loads(path.read_text())
        if not value['terminal']:raise ValueError('Complete measured inputs required')
        return value
    rates=read('matched_mode_directional_rates')['cases']
    fine=read('matched_mode_fine_controls')['cases']
    base='matched_mode_guarded_ne4_z32_'
    audits=[read(base+suffix+'_audit') for suffix in ('t30','t60','t30_eps2e4')]
    comparisons=[read('matched_mode_time_comparison'),read('matched_mode_amplitude_comparison')]
    identity=provenance(['examples.matched_mode_report']);destination.mkdir(parents=True);freeze(destination,identity)
    names=['density','radial_velocity','azimuthal_velocity','axial_velocity','electron_entropy','ion_entropy',
        'radial_B','azimuthal_B','axial_B','current','electric','electron_pressure','ion_pressure']
    labels=['$n$','$v_r$',r'$v_\theta$','$v_z$','$s_e$','$s_i$','$B_r$',r'$B_\theta$','$B_z$','$J$','$E$','$p_e$','$p_i$']
    fig,axes=plt.subplots(2,2,figsize=(13,8.5),layout='constrained')
    colors=['#176B93','#C26420','#4A7A46']; ax=axes[0,0]
    selected=[v for v in rates if v['spec']['field_nr']==32]
    for name,label in [('radial_B','Radial B rate'),('current_rate','Current rate')]:
        ax.loglog([v['spec']['field_nz'] for v in selected],
            [v['quadrature'][1]['rate_errors'][name]['relative_L2'] for v in selected],'-o',label=label)
    refined=fine[0]['quadrature'][1]['rate_errors']
    ax.scatter([128],[refined['current_rate']['relative_L2']],marker='D',color=colors[2],label='Current: radial 64')
    ax.set(xlabel='Axial field cells (radial 32 except diamond)',ylabel='Relative L2 error',title='A  Independent instantaneous mode rates')
    ax.set_xticks([32,64,128],labels=['32','64','128']);ax.xaxis.set_minor_locator(NullLocator());ax.legend(fontsize=8);ax.grid(True,which='both',alpha=.2)
    ax=axes[0,1];x=np.arange(len(names));legend=['30 steps, ε=10⁻⁴','60 steps, ε=10⁻⁴','30 steps, ε=2×10⁻⁴']
    for index,audit in enumerate(audits):
        values=audit['quadrature'][1]['physical']['final']['errors']
        ax.semilogy(x,[values[name]['relative_L2'] for name in names],'-o',color=colors[index],label=legend[index],ms=3)
    ax.set(xticks=x,xticklabels=labels,ylabel='Relative L2 error',title='B  Full 100 μs centered physical response')
    ax.legend(fontsize=8);ax.grid(True,which='both',alpha=.2)
    ax=axes[1,0]
    for value,label,color in zip(comparisons,['60 versus 30 steps','2ε versus ε'],colors):
        errors=value['quadrature'][1]['differences']
        ax.semilogy(x,[errors[name]['relative_L2'] for name in names],'-o',label=label,color=color,ms=3)
    ax.set(xticks=x,xticklabels=labels,ylabel='Difference / first response norm',title='C  Separate time and amplitude controls')
    ax.legend(fontsize=8);ax.grid(True,which='both',alpha=.2)
    ax=axes[1,1]
    for index,audit in enumerate(audits):
        energy=[abs(audit['quadrature'][1]['conservation'][sign]['relative_energy_drift']) for sign in ('minus','plus')]
        ax.semilogy([index-.12,index+.12],energy,'o',color=colors[index],label=legend[index])
    ax.set(xticks=[0,1,2],xticklabels=['30 steps','60 steps','Double amplitude'],ylabel='|Δ physical energy| / initial total energy',
        title='D  Independent continuous energy (both signs)')
    ax.grid(True,which='both',alpha=.2)
    fig.suptitle('Growing-mode verification: filled periodic cylinder, physical masses\nSeparate field refinement; coarse trajectory mesh nr=8, nz=32, fluid elements=4',fontsize=13)
    fig.savefig(destination/'convergence.png',dpi=170);fig.savefig(destination/'convergence.pdf');plt.close(fig)
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_sha256=inputs,
        independent_quadrature=20,rate_cases=8,paired_trajectories=3,
        limitations=['Coarse trajectory mesh; no full trajectory mesh convergence yet',
            'Two time resolutions do not determine a convergence order',
            'Energy normalized by full background energy; not a small-signal accuracy certificate',
            'No formation, exterior vacuum, end loss or laboratory confinement result'],
        files={name:sha(destination/name) for name in ('convergence.png','convergence.pdf')})
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    import sys
    run(sys.argv[1])
