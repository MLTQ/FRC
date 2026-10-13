"""Scientific report of physical spatial contributions and current charge defects."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .implicit_plasma_records import provenance,freeze,sha


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh report destination required')
    paths={key:Path('outputs')/name/'summary.json' for key,name in dict(
        parts='trajectory_spatial_parts',audit='trajectory_spatial_parts_audit',charge='trajectory_current_charge',
        original='trajectory_tangent_audit',parts_replay='trajectory_spatial_parts_replay_verification',
        audit_replay='trajectory_spatial_parts_audit_replay_verification',charge_replay='trajectory_current_charge_replay_verification').items()}
    records={key:json.loads(path.read_text()) for key,path in paths.items()}
    if any(not value['terminal'] for value in records.values()):raise ValueError('Complete scientific inputs required')
    for key in ('parts','audit','charge'):
        replay=records[key+'_replay']
        if not replay['bitwise_metrics'] or sha(paths[key])!=replay['original_summary_sha256']:
            raise ValueError('Original full replay proof required')
        for side in ('original','repeated'):
            for name,row in replay['arrays'].items():
                if sha(Path(replay[side+'_directory'])/name)!=row[side+'_sha256']:raise ValueError('Original replay arrays required')
    inputs={str(path):sha(path) for path in paths.values()}
    identity=provenance(['examples.spatial_parts_report']);destination.mkdir(parents=True);freeze(destination,identity)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'semibold'})
    fig,axes=plt.subplots(2,2,figsize=(12.6,9.2));fig.subplots_adjust(left=.09,right=.97,top=.87,bottom=.14,hspace=.52,wspace=.30)
    fig.suptitle('Where the rapid spatial response enters',fontsize=20,x=.06,ha='left',y=.972)
    fig.text(.06,.918,'Original F8/64, B32/64 equations · h = 10⁻⁵ · physical quadrature 20 · complete fresh replays',color='#526371',fontsize=11)
    names=['ideal_lorentz','hall','electron_species','ion_species'];labels=['Ideal / Lorentz','Hall','Electron species','Ion species']
    colors=['#3b8e9c','#5477b6','#cf6944','#b79b4d']
    q=next(q for q in records['audit']['quadrature'] if q['order']==20)
    for i,direction in enumerate(('saved_error','independent_mode')):
        case=next(c for c in q['cases'] if c['direction']==direction and c['amplitude']==1e-5)
        groups={g['name']:g for g in case['groups']}
        values=[groups[n]['metrics']['azimuthal_B']['projection_per_second']/1000 for n in names]
        ax=axes[0,i];ax.barh(labels,values,color=colors,height=.65);ax.axvline(0,color='#80909a',linewidth=.8);ax.invert_yaxis()
        for j,v in enumerate(values):ax.annotate(f'{v:,.2f}',(v,j),xytext=(5 if v>=0 else -5,0),textcoords='offset points',ha='left' if v>=0 else 'right',va='center',fontsize=9)
        ax.set_xlim(min(values)*1.45,max(values)*1.26);ax.set_xlabel('Bθ directional contribution (10³ / s)')
        ax.set_title(['Observed error direction','Independent continuum-mode direction'][i],loc='left',pad=11)
        ax.text(.01,-.28,'Signed contributions can cancel. These are not eigenvalues.',transform=ax.transAxes,fontsize=9,color='#5b6772')
    ax=axes[1,0];keys=['radius_times_divergence_over_current','normalized_radial_jump','normalized_axial_jump']
    x=np.arange(3);width=.35
    for i,(direction,color,label) in enumerate([('saved_error','#cf6944','Saved error'),('independent_mode','#3b8e9c','Independent mode')]):
        case=next(c for c in records['charge']['cases'] if c['order']==20 and c['direction']==direction)
        values=[case[k] for k in keys];ax.bar(x+(i-.5)*width,values,width,color=color,label=label)
    ax.set_yscale('log');ax.set_xticks(x,['Volume div J','Radial Jn jumps','Axial Jn jumps']);ax.set_ylim(1e-9,2)
    ax.set_ylabel('Dimensionless physical norm');ax.set_title('Local current continuity is unresolved',loc='left',pad=12)
    ax.legend(frameon=False,fontsize=9);ax.grid(axis='y',alpha=.15)
    ax.text(0,-.29,'Volume: R‖div J‖ / ‖J‖.  Surface: √R ‖[Jn]‖ / ‖J‖.\nCellwise Gauss-law relative defects ≤ 5.39 × 10⁻¹⁶.',transform=ax.transAxes,fontsize=9,color='#5b6772')
    ax=axes[1,1];case=next(c for c in q['cases'] if c['direction']=='independent_mode' and c['amplitude']==1e-5)
    groups={g['name']:g for g in case['groups']}
    old=next(c for q0 in records['original']['quadrature'] if q0['order']==20 for c in q0['cases'] if c['direction']=='independent_mode' and c['amplitude']==1e-5)
    values=[groups[n]['independent_continuum_errors']['current_axial']['relative_L2'] for n in names[:3]]
    values.append(old['continuum_rate_errors']['current_axial']['relative_L2'])
    ax.bar(['Ideal','Hall','Electron','Full response'],values,color=[*colors[:3],'#303f4b'])
    ax.set_yscale('log');ax.set_ylim(5e-5,.02);ax.set_ylabel('Relative L² error in Jz response')
    for k,v in enumerate(values):ax.text(k,v*1.22,f'{v:.2e}',ha='center',fontsize=9)
    ax.set_title('Known-mode cancellation and accuracy',loc='left',pad=12);ax.grid(axis='y',alpha=.15)
    ax.text(0,-.29,'Each bar uses its own continuum target norm.\nSeparated solves retain measured superposition error.',transform=ax.transAxes,fontsize=9,color='#5b6772')
    fig.text(.06,.045,'The mechanism is still under investigation. No term has been removed; no diffusion, filtering, floor or mass change.\nStrong interior and trace measurements do not yet establish a complete weak entropy-flux budget or FRC formation.',fontsize=10,color='#465763')
    fig.savefig(destination/'verification.png',dpi=170);fig.savefig(destination/'verification.pdf');plt.close(fig)
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Original report inputs changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_sha256=inputs,
                figures={p.name:sha(p) for p in destination.glob('verification.*')},scope='Scientific attribution report; interpretation limits remain explicit.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(destination/'verification.png',flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
