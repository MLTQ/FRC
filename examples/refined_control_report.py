"""Consolidate fully replayed refined controls and render their physical errors."""
import argparse
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, audit_directory):
    destination = Path(destination); audit_directory = Path(audit_directory); inputs = {}
    if destination.exists(): raise ValueError('Fresh report destination required')
    def check(path, expected=None):
        path = Path(path); actual = sha(path)
        if expected is not None and actual != expected: raise ValueError('Changed evidence: '+str(path))
        inputs[str(path)] = actual
    def read(path):
        check(path); value = json.loads(Path(path).read_text())
        if not value.get('terminal'): raise ValueError('Terminal evidence required: '+str(path))
        for name, expected in value.get('input_hashes', {}).items(): check(name, expected)
        if 'arrays_file' in value: check(Path(path).parent/value['arrays_file'], value['arrays_sha256'])
        return value
    def replay(path):
        value = read(path)
        if not value['bitwise_metrics']: raise ValueError('Full scientific replay required')
        for prefix in ('original','repeated'):
            check(Path(value[prefix+'_directory'])/'summary.json', value[prefix+'_summary_sha256'])
        for name, row in value['arrays'].items():
            if not row['bitwise_equal']: raise ValueError('Every replay array must match')
            for prefix in ('original','repeated'):
                check(Path(value[prefix+'_directory'])/name, row[prefix+'_sha256'])
        return value
    def independent(path):
        value = read(Path(path)/'summary.json'); replay(str(path)+'_verification/summary.json')
        for q in value.get('quadrature',[]):
            if not q['original_fingerprints_reproduced']: raise ValueError('Original reconstruction required')
        return value
    audits = read(audit_directory/'summary.json'); check(audits['campaign'],audits['campaign_sha256'])
    campaign = read(audits['campaign'])
    if not campaign['all_jobs_terminal']: raise ValueError('All physical jobs must be terminal')
    for row in campaign['jobs']:
        if row['returncode'] != 0: raise ValueError('Uncollected or interrupted physical protocol')
        check(Path(row['job']['control'])/'summary.json',row['control_sha256']); check(row['log'],row['log_sha256'])
    identity = provenance(['examples.refined_control_report'])
    result = dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],baseline={},
        maximum_component_quadrature_difference=0.,maximum_rate_quadrature_difference=0.,
        component_quadrature_difference_by_region=dict.fromkeys(('volume','axis','wall','interior'),0.),
        maximum_spectral_fraction_quadrature_difference=0.,maximum_curl_quadrature_difference=0.,
        maximum_parseval_relative_defect=0.,maximum_curl_decomposition_relative_defect=0.)
    curls = {}
    for entry in audits['curls']:
        value = independent(entry['directory'])
        for key, case in zip(value['inputs'],value['cases']):
            if not all(q['original_fingerprints_reproduced'] for q in case['quadrature']):
                raise ValueError('Complete current reconstruction required')
            curls[key[0]] = case
            for axis,metric in case['quadrature'][1]['metrics'].items():
                for name,value in metric['relative_L2'].items():
                    result['maximum_curl_quadrature_difference']=max(result['maximum_curl_quadrature_difference'],
                        abs(value-case['quadrature'][0]['metrics'][axis]['relative_L2'][name]))
            for q in case['quadrature']:
                for metric in q['metrics'].values():
                    result['maximum_curl_decomposition_relative_defect'] = max(
                        result['maximum_curl_decomposition_relative_defect'],metric['energy_decomposition_relative_defect'])
    for grid in ('axial','radial'):
        base = 'outputs/component_'+grid+'320'
        simulation = read(base+'/summary.json')
        replay(base+('_replay_serial' if grid == 'axial' else '_replay')+'_verification/summary.json')
        result['baseline'][grid] = dict(directory=base,case=simulation['cases'][0],
            curl=independent(base+'_curl')['cases'][0])
    for entry in audits['cases']:
        label = entry['label']; grid = label.split('_')[0]; directory = entry['input']
        check(entry['physical_replay'],entry['physical_replay_sha256']); replay(entry['physical_replay'])
        calculation = read(directory+'/summary.json'); case = calculation['cases'][0]
        check(Path(directory)/case['arrays_file'],case['arrays_sha256'])
        controls = [read('outputs/refined_control_diagnostics/'+label+suffix+'_control/summary.json')
            for suffix in ('','_replay')]
        for control in controls:
            if (not control['resources_within_bounds'] or control['resource_guard_triggered']
                    or control['peak_rss_bytes'] > control['maximum_rss_bytes']
                    or control['minimum_free_bytes'] < control['threshold_free_bytes']):
                raise ValueError('Physical guards and final observations must pass')
        if not controls[1]['complete_protocol_replay_verified']: raise ValueError('Full protocol replay required')
        if calculation['all_cases_complete'] != entry['accepted']: raise ValueError('Acceptance mismatch')
        row = dict(label=label,directory=directory,accepted=entry['accepted'],case=case,
            resource_controls=controls); result['cases'].append(row)
        if not entry['accepted']:
            baseline = result['baseline'][grid]['directory']
            with np.load(Path(directory)/case['arrays_file']) as failed, np.load(baseline+'/case_00.npz') as original:
                same={}
                for key in failed.files:
                    if key not in original: same[key]=False;continue
                    a,b=failed[key],original[key]
                    same[key]=a.shape==b.shape and a.dtype==b.dtype and a.tobytes()==b.tobytes()
            if not all(same.values()): raise ValueError('Rejected saved inputs differ from baseline')
            row['saved_input_arrays_equal_baseline'] = same
            continue
        if case['status'] != 'complete' or case['max_induction_residual'] > case['spec']['tolerance']:
            raise ValueError('Accepted original equations required')
        component = independent(entry['components']); row['components'] = component['quadrature']
        if component['inputs'] != [[result['baseline'][grid]['directory'],0],[directory,0]]:
            raise ValueError('Same-grid baseline and control required')
        row['curl'] = curls[directory]
        checked_current=0
        for physical,current in zip(component['quadrature'],row['curl']['quadrature']):
            if physical['order'] != current['order']: raise ValueError('Matching quadrature orders required')
            for axis,metric in current['metrics'].items():
                values=physical['metrics']['current_'+axis]['regions']['volume']['squared_integrals']
                for first,second in (('second_error','current_error'),('reference','reference')):
                    if values[first] != metric['squared_integrals'][second]:
                        raise ValueError('Separate component and curl physical integrals must agree')
                    checked_current+=1
        row['current_component_curl_integrals_equal']=checked_current
        row['relative_instantaneous_work_defects'] = [abs(d['total_power'])/sum(abs(p) for p in d['powers'].values())
            for d in case['rate_diagnostics']]
        q16,q20 = component['quadrature']
        for name,metric in q20['metrics'].items():
            for region,values in metric['regions'].items():
                for key,value in values['relative_L2'].items():
                    earlier=q16['metrics'][name]['regions'][region]['relative_L2'][key]
                    if value is not None and earlier is not None:
                        result['maximum_component_quadrature_difference']=max(
                            result['maximum_component_quadrature_difference'],abs(value-earlier))
                        result['component_quadrature_difference_by_region'][region]=max(
                            result['component_quadrature_difference_by_region'][region],abs(value-earlier))
        for name,metric in case['quadrature'][1]['rate_errors'].items():
            result['maximum_rate_quadrature_difference']=max(result['maximum_rate_quadrature_difference'],
                abs(metric['relative_L2']-case['quadrature'][0]['rate_errors'][name]['relative_L2']))
        if 'spectrum' in entry:
            modes=independent(entry['spectrum']); row['spectrum']=modes['quadrature']; checked=0
            for physical,spectral in zip(component['quadrature'],modes['quadrature']):
                for name,terms in spectral['metrics'].items():
                    for key,value in terms.items():
                        if value['direct_squared_integral'] != physical['metrics'][name]['regions']['volume']['squared_integrals'][key]:
                            raise ValueError('Separate component/spectrum integrals must agree')
                        checked+=1
                        result['maximum_parseval_relative_defect']=max(result['maximum_parseval_relative_defect'],value['parseval_relative_defect'])
            row['independent_integrals_equal']=checked; row['above_reference_error_RMS']={}
            for name,terms in modes['quadrature'][1]['metrics'].items():
                for key,metric in terms.items():
                    for group,value in metric['fractions'].items():
                        result['maximum_spectral_fraction_quadrature_difference']=max(
                            result['maximum_spectral_fraction_quadrature_difference'],
                            abs(value-modes['quadrature'][0]['metrics'][name][key]['fractions'][group]))
                norm=terms['reference']['direct_squared_integral']
                values={key: math.sqrt(terms[key]['direct_squared_integral']*terms[key]['fractions']['above_reference_band']/norm)
                    for key in ('first_error','second_error','difference')}
                values['small_over_baseline']=values['second_error']/values['first_error'] if values['first_error'] else None
                row['above_reference_error_RMS'][name]=values
    destination.mkdir(parents=True); freeze(destination,identity)
    render(destination,result)
    for path,expected in inputs.items():
        if sha(path) != expected: raise ValueError('Evidence changed during reporting')
    if provenance(['examples.refined_control_report']) != identity: raise ValueError('Report source changed')
    result.update(terminal=True,input_hashes=inputs,artifact_hashes={name:sha(destination/name)
        for name in ('verification.png','verification.pdf')},
        scope='Replayed signed-rate controls, not time evolution or formation. Cell modes fold aliases; strong curl excludes face distributions. Rejections have no accepted rate. No unit tests.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result


def render(destination,result):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(left=.085,right=.98,bottom=.19,top=.81,hspace=.55,wspace=.30)
    fig.suptitle('Refined-grid amplitude controls',x=.085,ha='left',y=.97,fontsize=20)
    fig.text(.085,.925,'Full signed rates and independent component / magnetic-curl audits, each freshly replayed',fontsize=11)
    fig.text(.085,.89,'Fluid 64 / 256  |  radial-16, N128 continuum reference  |  original tolerance 1e-10  |  order-20 norms',fontsize=10)
    colors=('#31688e','#a26725','#22856a'); axes_names=('radial','azimuthal','axial'); names=(r'$J_r$',r'$J_\theta$',r'$J_z$')
    statuses=[]
    for col,grid in enumerate(('axial','radial')):
        rows=[r for r in result['cases'] if r['label'].startswith(grid+'_') and not r['label'].endswith('_tight') and r['accepted']]
        baseline=result['baseline'][grid]; records=[(1e-4,baseline['curl'])]+[(r['case']['epsilon'],r['curl']) for r in rows]
        records.sort(key=lambda r:r[0]);rows.sort(key=lambda r:r['case']['epsilon'])
        for axis,color,name in zip(axes_names,colors,names):
            for metric,marker,fill in (('current_error','o',color),('strong_error','o','white')):
                axes[0,col].loglog([r[0] for r in records],[100*r[1]['quadrature'][1]['metrics'][axis]['relative_L2'][metric]
                    for r in records],color=color,marker=marker,markerfacecolor=fill,lw=1 if fill != 'white' else 0,
                    label=name if fill != 'white' else None)
            axes[1,col].loglog([r['case']['epsilon'] for r in rows],[100*r['components'][1]['metrics']['current_'+axis]['regions']['volume']['relative_L2']['difference']
                for r in rows],color=color,marker='o',linestyle='none',label=name)
        title='B256 / 320: axial refinement' if grid=='axial' else 'B320 / 256: radial refinement'
        axes[0,col].set(title=title,ylabel='Component-relative L2 error (%)')
        axes[1,col].set(title='Actual change from epsilon = 1e-4',ylabel='Change / continuum component norm (%)')
        for ax in axes[:,col]:
            ax.set(xlabel='Perturbation amplitude epsilon',xticks=[1e-5,1e-4,1e-3]);ax.grid(which='both',alpha=.15)
        axes[0,col].legend(frameon=False,ncol=3,fontsize=9)
        tight=next(r for r in result['cases'] if r['label']==grid+'_tight')
        statuses.append(grid+': '+('accepted' if tight['accepted'] else 'strictly rejected'))
    fig.text(.085,.12,'Tighter 1e-11 controls at epsilon = 1e-4: '+ '; '.join(statuses)+'. Both protocols replay.',fontsize=10)
    fig.text(.085,.079,'Hollow markers: independent inside-cell curl of the magnetic rate; face distributions are excluded.',fontsize=10)
    fig.text(.085,.040,'Amplitude sensitivity does not establish time accuracy, physical formation, stability or particle confinement.',fontsize=10,color='#444444')
    for suffix in ('png','pdf'):fig.savefig(destination/('verification.'+suffix),dpi=160)
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('audits')
    args=parser.parse_args();run(args.destination,args.audits)
