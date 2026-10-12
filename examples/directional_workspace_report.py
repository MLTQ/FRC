"""Render replayed directional physical errors and complete-workflow memory."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, cases, resources):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Fresh report path required')
    inputs = {}
    def read(path):
        path = Path(path); result = json.loads(path.read_text()); inputs[str(path)] = sha(path)
        if not result.get('terminal'): raise ValueError('Terminal input required: '+str(path))
        return result
    def replay(path):
        value = read(path)
        if not value['bitwise_metrics'] or not all(a['bitwise_equal'] for a in value['arrays'].values()):
            raise ValueError('Complete scientific replay required')
        for prefix in ('original', 'repeated'):
            summary = Path(value[prefix+'_directory'])/'summary.json'
            if sha(summary) != value[prefix+'_summary_sha256']: raise ValueError('Replay input changed')
            inputs[str(summary)] = sha(summary)
        return value
    physical = []
    for item in cases:
        calculation = read(Path(item['calculation'])/'summary.json')
        audit = read(Path(item['audit'])/'summary.json')
        replay(Path(item['audit']+'_verification')/'summary.json')
        replay(item['verification'])
        if not calculation['all_cases_complete'] or len(calculation['cases']) != 1:
            raise ValueError('One complete signed pair required')
        for path, expected in audit['input_hashes'].items():
            if sha(path) != expected: raise ValueError('Physical audit input changed')
            inputs[path] = expected
        if audit['inputs'] != [[item['calculation'], 0]]: raise ValueError('Matching physical audit required')
        q = audit['cases'][0]['quadrature'][-1]; record = calculation['cases'][0]
        field = record['quadrature'][-1]
        if q['order'] != 20 or field['order'] != 20 or not q['original_fingerprints_reproduced']:
            raise ValueError('Full independent q20 reconstruction required')
        physical.append(dict(label=item['label'], spec=record['spec'],
            current=[100*q['metrics'][a]['relative_L2']['current_error'] for a in ('radial','azimuthal','axial')],
            strong=[100*q['metrics'][a]['relative_L2']['strong_error'] for a in ('radial','azimuthal','axial')],
            magnetic=[100*field['rate_errors'][a]['relative_L2'] for a in ('radial_B','azimuthal_B','axial_B')]))
    memory = []
    for item in resources:
        rows = [read(path) for path in item['controls']]
        if len(rows) != 2 or not rows[1].get('complete_replay_verified'):
            raise ValueError('Primary and full fresh replay resource records required')
        for row in rows:
            if (not row['all_cases_complete'] or row['resource_guard_triggered']
                    or row['peak_rss_bytes'] > row['maximum_rss_bytes']
                    or row['minimum_free_bytes'] < row['threshold_free_bytes']):
                raise ValueError('Resource guard and final bounds must pass')
        replay(rows[1]['destination']+'_verification/summary.json')
        memory.append(dict(label=item['label'], GiB=[r['peak_rss_bytes']/2**30 for r in rows]))
    identity = provenance(['examples.directional_workspace_report'])
    destination.mkdir(parents=True); freeze(destination, identity)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 6.5))
    fig.subplots_adjust(left=.060, right=.985, bottom=.22, top=.72, wspace=.34)
    fig.suptitle('Directional field refinement with physical replay controls', x=.06, ha='left', y=.97, fontsize=19)
    fig.text(.06,.89,'Fixed fluid 64 / 256  |  Independent radial-16, N128 reference  |  Original solver tolerance 1e-10',fontsize=11)
    fig.text(.06,.845,'Field-grid labels give radial / axial cell counts. Plotted errors use independent order-20 volume quadrature.',fontsize=10,color='#444444')
    colors = ('#355c7d','#437f97','#27806a','#bf7634','#8a6388'); width=.76/max(len(physical),1); x=np.arange(3)
    for k, row in enumerate(physical):
        pos=x+(k-(len(physical)-1)/2)*width
        for axis,key in zip(axes[:2],('current','magnetic')):
            axis.bar(pos,row[key],width,color=colors[k%len(colors)],label=row['label'])
        axes[0].scatter(pos,row['strong'],s=25,facecolors='white',edgecolors='black',zorder=3)
    for axis,labels,title in zip(axes[:2],((r'$J_r$',r'$J_\theta$',r'$J_z$'),(r'$B_r$',r'$B_\theta$',r'$B_z$')),
                                  ('Current-rate component errors','Magnetic-rate component errors')):
        axis.set(xticks=x,xticklabels=labels,yscale='log',ylabel='Relative physical L2 error (%)',title=title)
        axis.grid(axis='y',which='both',alpha=.15);axis.set_axisbelow(True)
    axes[0].legend(frameon=False,fontsize=8,ncol=2,loc='lower left',bbox_to_anchor=(-.02,1.07))
    x=np.arange(len(memory)); width=.34
    for k,label in enumerate(('Primary','Full replay')):
        bars=axes[2].bar(x+(k-.5)*width,[m['GiB'][k] for m in memory],width,color=colors[k*2],label=label)
        axes[2].bar_label(bars,fmt='%.2f',padding=3,fontsize=8)
    axes[2].axhline(8,color='#a64635',ls='--',lw=1,label='8 GiB guard')
    axes[2].set(xticks=x,xticklabels=[m['label'] for m in memory],ylim=(0,9.4),ylabel='Own-process peak RSS (GiB)',title='Complete workflow memory')
    axes[2].legend(frameon=False,fontsize=8,ncol=3,loc='lower center',bbox_to_anchor=(.5,1.07))
    axes[2].grid(axis='y',alpha=.15);axes[2].set_axisbelow(True)
    fig.text(.06,.115,'Hollow circles: independent strong curl of the actual magnetic rate, excluding face distributions. All runs and audits replay completely.',fontsize=10)
    fig.text(.06,.065,'Modest directional increments do not establish full convergence, time accuracy or formation. Memory observations vary between processes.',fontsize=10,color='#444444')
    for suffix in ('png','pdf'): fig.savefig(destination/('verification.'+suffix),dpi=160)
    plt.close(fig)
    if provenance(['examples.directional_workspace_report']) != identity: raise ValueError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_hashes=inputs,physical=physical,memory=memory,
        artifact_hashes={name:sha(destination/name) for name in ('verification.png','verification.pdf')})
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result
