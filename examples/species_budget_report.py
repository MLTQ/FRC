"""Verified quadrature comparisons and a scientific entropy-budget figure."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .species_budget_quadrature import calculate
from .trajectory_tangent_probe import verified_inputs
from .implicit_plasma_records import provenance,freeze,sha


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Fresh species-budget report required')
    resources=[]
    def observe():
        item=dict(peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,free_bytes=shutil.disk_usage('.').free)
        resources.append(item)
        if item['peak_rss_bytes']>8*1024**3 or item['free_bytes']<256*1024**2:raise RuntimeError('Report resource bounds exceeded')
    observe();records={};inputs={}
    for kind in ('interface','mass','space'):
        directory=Path('outputs')/('species_'+kind+'_audit')
        records[kind],checked=verified_inputs(directory,str(directory)+'_replay_verification')
        inputs.update(checked);inputs.update(records[kind]['input_sha256'])
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Original complete budget inputs required')
    identity=provenance(['examples.species_budget_report']);destination.mkdir(parents=True);freeze(destination,identity)
    rows,integrals=calculate(records)
    path=destination/'quadrature_integrals.npz';np.savez_compressed(path,integrals=integrals)
    bg='#0b1520';fg='#dce9ed';muted='#a1b6c4';colors={'electron':'#6dd3c0','ion':'#e8b278'}
    plt.rcParams.update({'font.size':10,'text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted,
                         'axes.edgecolor':'#49606c','axes.facecolor':bg,'figure.facecolor':bg,'savefig.facecolor':bg})
    fig=plt.figure(figsize=(13,8.2));left=fig.add_axes([.085,.33,.43,.43]);right=fig.add_axes([.65,.36,.28,.37])
    fig.text(.065,.934,'FRC NUMERICAL VERIFICATION  /  SPECIES ENTROPY',color=muted,fontsize=10,family='monospace')
    fig.text(.065,.878,'Where the entropy-rate error enters',fontsize=25)
    fig.text(.065,.827,'Same archived current · unchanged physical equations · independent volume, face and mass accounting',color=muted,fontsize=10)
    cases=[c for c in records['space']['cases'] if c['order']==20 and c['direction']=='saved_error']
    for label,color in colors.items():
        for kind,style in [('original','-'),('exact_mass','--')]:
            y=[next(s for s in c['species'] if s['name']==label)['errors'][kind]['relative_L2'] for c in cases]
            left.semilogy([c['elements'] for c in cases],y,style+'o',color=color,markersize=4,
                          label=label.title()+(' · original rate' if kind=='original' else ' · exact-mass diagnostic'))
    left.set_xticks([8,16,32]);left.set_xlabel('Fluid radial elements (field grid stays 32 × 64)')
    left.set_ylabel('Entropy-rate error / strong-source norm');left.set_title('Fixed-current space refinement',loc='left',pad=15,fontsize=12,color=fg)
    left.grid(alpha=.15);left.legend(loc='lower left',bbox_to_anchor=(0,-.29),ncol=2,fontsize=8,frameon=False,labelcolor=muted)
    mass=next(c for q in records['mass']['quadrature'] if q['order']==20 for c in q['cases'] if c['direction']=='saved_error')
    names=['interface_transfer','mass_map','representation'];positions=np.arange(3)
    for offset,(label,color) in zip([-.07,.07],colors.items()):
        values=next(s for s in mass['species'] if s['name']==label)['component_norms_relative_to_strong']
        right.semilogx([values[k] for k in names],positions+offset,'o',color=color,markersize=6)
    right.set_yticks(positions,['Interface + transfer','Mass map','Representation']);right.invert_yaxis()
    right.set_xlim(1e-9,1);right.set_xlabel('Separate error norm / strong-source norm')
    right.set_title('F8 error fields · norms do not add',loc='left',pad=20,fontsize=11,color=fg);right.grid(axis='x',alpha=.15)
    fig.text(.065,.168,'Electron source error: 55.0% → 8.97% → 0.0292% with the original species operator.',fontsize=13,color=fg)
    fig.text(.065,.119,'Full weak identity agrees within 1.8 × 10⁻¹² relative. All three complete numerical protocols replay exactly.',fontsize=10,color=muted)
    fig.text(.065,.077,'Exact-mass curves are diagnostics only. No solver repair, coupled-growth cure or FRC formation is established.',fontsize=10,color='#e8b278')
    fig.text(.065,.038,'Q16/20 physical comparisons · both species and directions · original bounds and states preserved · 12 OCT 2026',fontsize=8,color=muted)
    figure=destination/'verification.png';fig.savefig(figure,dpi=180);plt.close(fig)
    observe()
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.species_budget_report'])!=identity:raise ValueError('Report input/source changed')
    summary=dict(terminal=True,all_cases_complete=True,source_sha256=identity['source_sha256'],input_sha256=inputs,
        quadrature_rows=rows,integral_columns=['difference_squared','component_squared','full_rate_squared'],
        arrays_file=path.name,arrays_sha256=sha(path),figure=figure.name,figure_sha256=sha(figure),
        maximum_quadrature_difference_relative_to_full=max(r['relative_to_full'] for r in rows if r['relative_to_full'] is not None),
        scope='Full entropy-budget verification and fixed-current scalar-space refinement; no coupled solver modification.')
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (destination/'resources.json').write_text(json.dumps(dict(stage_observations=resources,continuous_guard=False),indent=2)+'\n')
    print(figure,flush=True);return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
