"""Scientific figure for fixed-wavelength accuracy and strict solve recovery."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from .implicit_plasma_records import provenance,freeze,sha


def run(destination='outputs/radial16_report'):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    names=['radial16_reference_audit','radial16_matched_rates','radial16_fine_rates','radial16_balanced_rates','radial16_balanced_probe']
    sources={n:Path('outputs')/n/'summary.json' for n in names}
    data={n:json.loads(p.read_text()) for n,p in sources.items()}
    if not all(s['terminal'] for s in data.values()):raise ValueError('Terminal inputs required')
    cases=data[names[1]]['cases'][:6]+data[names[2]]['cases'][:1]+data[names[3]]['cases'][1:2]
    if any(c['status']!='complete' for c in cases):raise ValueError('Passing rates required')
    q=[c['quadrature'][-1] for c in cases];x=np.arange(len(q))
    labels=[f"F{c['spec']['elements']} B{c['spec']['field_nr']}\nz{c['spec']['points']}/{c['spec']['field_nz']}" for c in cases]
    identity=provenance(['examples.radial16_report']);destination.mkdir(parents=True);freeze(destination,identity)
    plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(3,2,figsize=(14,13),layout='constrained')
    colors=['#176B87','#C85C32','#617A35']
    def series(ax,values,label,color):ax.semilogy(x,values,'o-',label=label,color=color,lw=1.8,ms=4)
    for key,label,col in zip(['radial_velocity','azimuthal_B','current_rate'],['Radial velocity rate','Toroidal B rate','Current rate'],colors):
        series(axes[0,0],[r['rate_errors'][key]['relative_L2'] for r in q],label,col)
    axes[0,0].axhline(1,color='#777777',ls=':',lw=1)
    axes[0,0].set_title('A. Original-equation rate error remains substantial')
    for key,label,col in zip(['density','axial_B','initial_current'],['Density','Axial B','Weak current'],colors):
        values=[r['initial_errors'][key]['relative_L2'] if key!='initial_current' else r['representation_and_faces']['volume'][key]['relative_L2'] for r in q]
        series(axes[0,1],values,label,col)
    axes[0,1].set_title('B. Initial representation converges faster than rates')
    for direction,key,label,col in [('radial','normal_current','Radial normal-J jump',colors[0]),('radial','tangential_B','Radial tangential-B jump',colors[1]),('axial','normal_current','Axial normal-J jump',colors[2])]:
        series(axes[1,0],[r['representation_and_faces']['internal_faces'][direction][key]['jump_over_reference_L2'] for r in q],label,col)
    axes[1,0].set_title('C. Both one-sided internal traces are retained')
    for key,label,col in zip(['radial_B','axial_B','electric'],['Radial B rate','Axial B rate','Full electric field'],colors):
        series(axes[1,1],[r['rate_errors'][key]['relative_L2'] for r in q],label,col)
    axes[1,1].set_title('D. Separate fluid / field / axial refinements matter')
    for ax in axes[:2].ravel():
        ax.set_xticks(x,labels,fontsize=8);ax.set_ylabel('Relative physical L2 error');ax.legend(fontsize=8)
        ax.grid(alpha=.2);ax.set_xlabel('F = fluid radial elements; B = field radial cells; z = fluid / field axial points')
    axes[1,0].set_ylabel('Jump / analytic face L2 norm')
    ref=data[names[0]]['cases'];a=axes[2,0]
    a.semilogy([r['points'] for r in ref],[r['curl_controls'][1]['component_relative_errors'][0] for r in ref],'o-',color=colors[0],label='Independent curl(B)/mu0 vs Jr')
    a.set_xticks([64,128,256]);a.set_xlabel('Independent axial reference points');a.set_ylabel('Relative current mismatch');a.grid(alpha=.2);a.legend(fontsize=8)
    a.set_title('E. N64 growth convergence hid a Nyquist derivative defect')
    trace=data[names[4]]['trace'][0]['trace'];a=axes[2,1]
    for key,label,col in zip(['operator','faraday_dual','faraday_physical'],['Original operator','Dual Faraday','Physical Faraday'],colors):
        a.semilogy(range(len(trace)),[r[key] for r in trace],'o-',color=col,label=label)
    a.axhline(1e-10,color='#222222',ls='--',label='Unchanged acceptance limit');a.set_xticks(range(len(trace)))
    a.set_xlabel('Inner refinement');a.set_ylabel('Original residual');a.set_title('F. Saved failed input recovered with all bounds enforced');a.grid(alpha=.2);a.legend(fontsize=8)
    fig.suptitle('Fixed radial16 verification: better resolution, accuracy gate still open',fontsize=17)
    fig.supxlabel('Filled periodic physical-mass fixture; no added diffusion, atmosphere, filtering or energy repair.\nInstantaneous consistency and strict solve recovery do not establish FRC formation or confinement.',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(destination/f'verification.{ext}',dpi=150)
    plt.close(fig)
    record=dict(terminal=True,source_sha256=identity['source_sha256'],inputs={str(p):sha(p) for p in sources.values()},
        selected_specs=[c['spec'] for c in cases],figures={f'verification.{e}':sha(destination/f'verification.{e}') for e in ('png','pdf')})
    (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
