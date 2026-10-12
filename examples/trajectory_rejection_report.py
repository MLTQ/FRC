"""Scientific figure for a failed refined trajectory and replayed density audit."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from .implicit_plasma_records import provenance,freeze,sha


def run(localization,verification,destination):
    localization=Path(localization);verification=Path(verification);destination=Path(destination)
    if destination.exists():raise ValueError('Fresh scientific-report destination required')
    observed=json.loads((localization/'summary.json').read_text())
    replay=json.loads((verification/'summary.json').read_text())
    if not observed['terminal'] or not replay['terminal'] or not replay['bitwise_metrics']:
        raise ValueError('Complete replayed density localization required')
    if replay['original_directory']!=str(localization):raise ValueError('Matching verification required')
    inputs={str(localization/'summary.json'):sha(localization/'summary.json'),
        str(verification/'summary.json'):sha(verification/'summary.json')}
    for prefix in ('original','repeated'):
        path=Path(replay[prefix+'_directory'])/'summary.json'
        if sha(path)!=replay[prefix+'_summary_sha256']:raise ValueError('Unchanged replay summary required')
        inputs[str(path)]=sha(path)
    for filename,row in replay['arrays'].items():
        if row['original_sha256']!=row['repeated_sha256']:raise ValueError('Byte-identical localization replay required')
        for prefix in ('original','repeated'):
            path=Path(replay[prefix+'_directory'])/filename
            if sha(path)!=row[prefix+'_sha256']:raise ValueError('Unchanged replay arrays required')
            inputs[str(path)]=sha(path)
    arrays=[]
    for case in observed['cases']:
        path=localization/case['arrays_file']
        if sha(path)!=case['arrays_sha256']:raise ValueError('Original localization arrays required')
        arrays.append(dict(np.load(path)))
    failed=Path(observed['cases'][0]['directory'])/'summary.json'
    physical=json.loads(failed.read_text());inputs[str(failed)]=sha(failed)
    if not physical['terminal'] or physical['success'] or physical.get('failed_step')!=18 or physical['completed_steps']!=17:
        raise ValueError('This report describes the recorded17-step/step18 rejection only')
    identity=provenance(['examples.trajectory_rejection_report']);destination.mkdir(parents=True);freeze(destination,identity)
    plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.spines.top':False,
        'axes.spines.right':False,'figure.facecolor':'#fbfcfe','axes.facecolor':'#fbfcfe'})
    fig,axes=plt.subplots(2,2,figsize=(13,9.8));colors=['#ba3933','#24738a','#8c719e']
    labels=['Fine split: F8/64, B32/64','Older medium: F4/64, B16/64','Coarse split: F4/32, B8/32']
    metrics=[next(q for q in case['quadrature'] if q['order']==20)['metrics'] for case in observed['cases']]
    for case,data,label,color in zip(metrics,arrays,labels,colors):
        t=1e6*data['time']
        axes[0,0].semilogy(t,[row['relative_L2_error'] for row in case],label=label,color=color,lw=2)
        axes[1,1].semilogy(t[1:],[100*row['sampled_max_interface_jump_over_background'] for row in case[1:]],color=color,lw=2)
    axes[0,0].set(title='Density departure from the independent linear target',
        xlabel='Physical time (µs)',ylabel='Density L2 error / growing linear signal')
    axes[0,0].legend(fontsize=8.5,loc='upper left')
    edges=arrays[0]['fluid_radial_edges'];centers=(edges[:-1]+edges[1:])/(2*edges[-1])
    fractions=np.array(metrics[0][-1]['element_error_fractions'])
    axes[0,1].bar(centers,100*fractions,width=.85*np.diff(edges)/edges[-1],color=colors[0])
    axes[0,1].set(title='Final refined error occupies several radial elements',
        xlabel='Radius / wall radius',ylabel='Share of volume squared density error (%)',xlim=(0,1))
    t=1e6*arrays[0]['time']
    axes[1,0].semilogy(t,[row['low_mode_fraction'] for row in metrics[0]],color='#276b48',lw=2,label='Low axial cell modes: |k| ≤ 4')
    axes[1,0].semilogy(t,[row['upper_quarter_mode_fraction'] for row in metrics[0]],color='#b48b34',lw=2,label='Upper quarter: |k| ≥ 24')
    axes[1,0].set(title='The refined error is dominated by low axial cell modes',
        xlabel='Physical time (µs)',ylabel='Fraction of volume squared density error')
    axes[1,0].legend(fontsize=8.5,loc='lower left')
    axes[1,1].set(title='Internal density traces separate as the error grows',
        xlabel='Physical time (µs)',ylabel='Maximum sampled jump / background density (%)')
    for ax in axes.flat:ax.grid(alpha=.18,which='major',axis='both');ax.set_axisbelow(True)
    fig.suptitle('Refined contact trajectory: rapid growth and step-18 rejection',x=.055,y=.974,ha='left',fontsize=17,weight='bold')
    fig.text(.055,.932,'100 µs requested; 17 accepted steps through 28.33 µs. Original tolerances retained. Quadrature 20 shown.',fontsize=10.5)
    fig.subplots_adjust(left=.075,right=.975,top=.855,bottom=.16,hspace=.38,wspace=.27)
    fig.text(.055,.105,'Independent density observations and localization replay bitwise; the complete failed-trajectory replay is a separate gate.',fontsize=10)
    fig.text(.055,.079,'Single sign and a linear reference. Cell spectra fold aliases. Initial tiny-error fractions are quadrature-sensitive.',fontsize=9.5)
    fig.text(.055,.053,'The older medium control uses different induction arithmetic. Location and trace jumps do not identify the cause or prove a flux-budget defect.',fontsize=9.5)
    fig.text(.055,.027,'Filled periodic verification fixture; no formation, laboratory stability or particle-confinement claim.',fontsize=9.5,color='#4c5665')
    for suffix in ('png','pdf'):fig.savefig(destination/f'verification.{suffix}',dpi=170)
    plt.close(fig)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.trajectory_rejection_report'])!=identity:
        raise ValueError('Report inputs or source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_sha256=inputs,
        figure_sha256={name:sha(destination/name) for name in ('verification.png','verification.pdf')},
        scope='Rendering of complete replayed density localization and the terminal primary numerical rejection; no assertion that the separate full rejection replay is complete.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('localization');parser.add_argument('verification');parser.add_argument('destination')
    args=parser.parse_args();run(args.localization,args.verification,args.destination)
