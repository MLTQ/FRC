"""Independent saved-state continuous audits and scientific convergence figure."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from .cylindrical_entropy_study import initialize
from .cylindrical_entropy_reference import regular_sample,physical
from .radial_work_records import sha


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
    root=args.directory;out=args.output;out.mkdir(parents=True,exist_ok=True)
    rates=[json.loads(p.read_text()) for p in sorted(root.glob('rates_*.json'))]
    spectra=[json.loads(p.read_text()) for p in sorted(root.glob('spectra_*.json'))]
    fluids=[json.loads(p.read_text()) for p in sorted(root.glob('fluid_[0-9][0-9].json'))]
    audit=json.loads((root/'audit.json').read_text());continuous=[]
    for record in fluids:
        ref,d,_=initialize(record['spec']);arrays=np.load(root/record['arrays_file']);rows=[]
        for state,step in zip(arrays['snapshots'],arrays['snapshot_steps'],strict=True):
            energy=[];moment=[]
            for order,factor in ((32,4),(48,8)):
                u,r,z,volume=regular_sample(d,state,order,factor);f=physical(u,r,ref)
                if not all(np.isfinite(v).all() for v in f.values()) or np.min(f['pressure'])<=0:raise ValueError('Invalid continuous physical sample')
                energy.append(float(np.sum(volume*f['energy'])));moment.append(np.sum(volume[...,None]*u,axis=(0,1,2)))
            nodal_energy=float(np.sum(d.volume*physical(state,d.r,ref)['energy']))
            nodal_moment=np.sum(d.volume[...,None]*state,axis=(0,1,2))
            rows.append(dict(step=int(step),continuous_energy=energy[1],nodal_energy=nodal_energy,
                        quadrature_difference=abs(energy[0]-energy[1])/abs(energy[1]),nodal_continuous_energy_difference=(nodal_energy-energy[1])/energy[1],
                        continuous_moments=moment[1].tolist(),nodal_moments=nodal_moment.tolist()))
        continuous.append(dict(file=record['arrays_file'],spec=record['spec'],states=rows,
            continuous_energy_drift=(rows[-1]['continuous_energy']-rows[0]['continuous_energy'])/rows[0]['continuous_energy']))
    summary=dict(source_sha256=audit['source_sha256'],rate_cases=len(rates),spectral_cases=len(spectra),trajectories=len(fluids),
                 maximum_energy_rate_defect=max(r['global_physical_energy_rate_relative'] for r in rates),
                 maximum_conserved_rate_defect=max(max(r['conserved_rate_relative']) for r in rates),
                 maximum_history_conserved_drift=max(max(r['history_conserved_drift']) for r in fluids),
                 maximum_independent_local_energy_difference=max(r.get('independent_local_energy_rate_relative',0) for r in audit['records']),
                 maximum_saved_state_quadrature_difference=max(s['quadrature_difference'] for r in continuous for s in r['states']),
                 maximum_saved_state_nodal_continuous_energy_difference=max(abs(s['nodal_continuous_energy_difference']) for r in continuous for s in r['states']),
                 continuous=continuous)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,3,figsize=(14,8.4),layout='constrained');colors={3:'#cc6e30',5:'#187e99',7:'#6b56a5'}
    for degree in (3,5,7):
        chosen=[r for r in rates if r['spec']['degree']==degree and r['spec']['points']==128];x=[r['spec']['elements'] for r in chosen]
        for ax,key in ((axes[0,0],'velocity'),(axes[0,1],'axis_radial')):
            ax.loglog(x,[r['errors'][key]['relative_L2'] for r in chosen],'-o',color=colors[degree],label=f'degree {degree}')
    axes[0,0].set(title='Original 2-D acceleration',xlabel='Radial elements (128 axial points)',ylabel='Whole-domain relative L2 error')
    axes[0,1].set(title='Central-axis accuracy',xlabel='Radial elements (128 axial points)',ylabel='First-element radial relative L2 error')
    for degree in (3,5):
        chosen=[r for r in spectra if r['spec']['degree']==degree]
        axes[0,2].loglog([r['spec']['elements'] for r in chosen],[abs(r['first_three_frequency_relative_errors'][0]) for r in chosen],'-o',color=colors[degree],label=f'degree {degree}')
    axes[0,2].set(title='Cylindrical sound frequency',xlabel='Radial elements',ylabel='Relative error vs first J1 zero')
    for key,label in (('n','Density'),('temperature','Species temperature')):
        axes[1,0].loglog([r['spec']['points'] for r in fluids[:3]],[r['final_metrics'][key] for r in fluids[:3]],'-o',label=label)
    axes[1,0].set(title='One complete contact transit',xlabel='Periodic axial points',ylabel='Continuous-field relative L2 error')
    chosen=[fluids[i] for i in (1,10,15)];x=np.array([r['dt'] for r in chosen]);y=np.array([abs(r['physical_energy_drift']) for r in chosen])
    axes[1,1].loglog(x/x[0],y,'-o',color=colors[5],label='Actual physical energy');axes[1,1].loglog(x/x[0],y[0]*(x/x[0])**4,'--',color='#777777',label='Fourth-order reference')
    axes[1,1].set(title='Time-step control at 32 axial points',xlabel='Time step / base time step',ylabel='Absolute relative energy drift')
    chosen=[fluids[i] for i in (12,13,14)]
    axes[1,2].loglog([r['spec']['amplitude'] for r in chosen],[r['density_perturbation_error'] for r in chosen],'-o',color=colors[7],label='Density error / initial perturbation')
    axes[1,2].set(title='One radial sound-wave period',xlabel='Initial fractional amplitude',ylabel='Departure from linearized Bessel solution')
    for ax in axes.flat:ax.grid(True,which='major',alpha=.2);ax.legend(fontsize=8,frameon=False)
    for ax,ticks in ((axes[0,0],[2,4,8,16,32]),(axes[0,1],[2,4,8,16,32]),(axes[0,2],[2,4,8]),
                     (axes[1,0],[16,32,64]),(axes[1,1],[.25,.5,1])):
        ax.set_xticks(ticks,[str(t) for t in ticks]);ax.xaxis.set_minor_locator(NullLocator())
    fig.suptitle('Cylindrical transport: physical energy and angular momentum',fontsize=18,fontweight='bold')
    fig.supxlabel('Unforced fluid verification • regular axis • impermeable wall • periodic z • no magnetic coupling yet',fontsize=10)
    fig.savefig(out/'report.png',dpi=180);fig.savefig(out/'report.pdf');plt.close(fig)
    (out/'report_source.json').write_text(json.dumps(dict(script_sha256=sha(__file__),audit_sha256=sha(root/'audit.json')),indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='continuous'},indent=2))


if __name__=='__main__':main()
