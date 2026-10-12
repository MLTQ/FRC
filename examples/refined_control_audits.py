"""Complete independent component/curl/spectral audits of replayed controls."""
import argparse
import json
from pathlib import Path
import resource
import shutil
from .component_mode_audit import run as components
from .component_axial_spectrum import run as spectrum
from .current_rate_curl_audit import run as curl
from .wall_diagnostic_replay import verify
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, campaign_directory):
    destination = Path(destination); campaign_directory = Path(campaign_directory)
    if destination.exists(): raise ValueError('Fresh independent-audit destination required')
    campaign_path = campaign_directory/'summary.json'; campaign = json.loads(campaign_path.read_text())
    if not campaign['terminal'] or not campaign.get('all_jobs_terminal'):
        raise ValueError('Complete serial physical campaign required')
    source = json.loads((campaign_directory/'provenance.json').read_text())
    if any(sha(p) != h for p,h in source['source_files'].items()): raise ValueError('Unchanged campaign source required')
    identity = provenance(['examples.refined_control_audits'])
    destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=False, source_sha256=identity['source_sha256'], campaign=str(campaign_path),
        campaign_sha256=sha(campaign_path), cases=[], curls=[], active_stage=None,
        minimum_free_bytes=shutil.disk_usage('.').free, peak_rss_bytes=0)
    def save():
        result['minimum_free_bytes'] = min(result['minimum_free_bytes'],shutil.disk_usage('.').free)
        result['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    def stage(name, function, *args):
        result['active_stage'] = name; save(); print(name,flush=True)
        if result['minimum_free_bytes'] < 256*1024**2: raise RuntimeError('Audit disk preflight below256MiB')
        value = function(*args); save(); return value
    save()
    for grid in ('axial','radial'):
        baseline = 'outputs/component_'+grid+'320'; accepted = []
        for size in ('small','large','tight'):
            label = grid+'_'+size
            job = next(row for row in campaign['jobs'] if row['job']['name'] == label+'_replay')
            control_path = Path(job['job']['control'])/'summary.json'
            if sha(control_path) != job['control_sha256']: raise ValueError('Original control changed')
            control = json.loads(control_path.read_text())
            if not control.get('complete_protocol_replay_verified'): raise ValueError('Independent complete replay required')
            original = job['job']['original']; verification = Path(job['job']['destination']+'_verification')/'summary.json'
            replay = json.loads(verification.read_text())
            if not replay['terminal'] or not replay['bitwise_metrics']: raise ValueError('Full scientific replay required')
            for prefix in ('original','repeated'):
                if sha(Path(replay[prefix+'_directory'])/'summary.json') != replay[prefix+'_summary_sha256']:
                    raise ValueError('Scientific replay summary changed')
            row = dict(label=label, input=original, physical_replay=str(verification),
                physical_replay_sha256=sha(verification), accepted=control['accepted_physical_replay'])
            result['cases'].append(row); save()
            if not row['accepted']:
                row['scope']='Replayed numerical rejection has no accepted rate to audit';save();continue
            accepted.append((original,0)); path=destination/(label+'_components')
            first=(baseline,0); second=(original,0)
            measured=stage(label+' components',components,path,first,second)
            stage(label+' component replay',components,str(path)+'_replay',first,second)
            stage(label+' component verification',verify,path,str(path)+'_replay',str(path)+'_verification')
            row['components']=str(path)
            if size == 'small':
                path=destination/(label+'_spectrum')
                modes=stage(label+' spectrum',spectrum,path,first,second)
                stage(label+' spectrum replay',spectrum,str(path)+'_replay',first,second)
                stage(label+' spectrum verification',verify,path,str(path)+'_replay',str(path)+'_verification')
                checked=0
                for physical,spectral in zip(measured['quadrature'],modes['quadrature']):
                    if physical['order'] != spectral['order']: raise ValueError('Same quadrature required')
                    for name,terms in spectral['metrics'].items():
                        for key,value in terms.items():
                            expected=physical['metrics'][name]['regions']['volume']['squared_integrals'][key]
                            if expected != value['direct_squared_integral']: raise ValueError('Independent physical integrals differ')
                            checked+=1
                row.update(spectrum=str(path),direct_component_integrals_bitwise_equal=checked)
            save()
        if accepted:
            path=destination/(grid+'_curl')
            stage(grid+' magnetic/current curl',curl,path,accepted)
            stage(grid+' magnetic/current curl replay',curl,str(path)+'_replay',accepted)
            stage(grid+' magnetic/current curl verification',verify,path,str(path)+'_replay',str(path)+'_verification')
            result['curls'].append(dict(directory=str(path),inputs=accepted));save()
    if provenance(['examples.refined_control_audits']) != identity or sha(campaign_path) != result['campaign_sha256']:
        raise ValueError('Source or campaign changed')
    result.update(terminal=True,active_stage=None,
        scope='Complete independent audits and replays of accepted signed rates; rejected protocols are not assigned physical errors. Spectra are cell-translation modes with aliases, and strong curl excludes face distributions. Resources are observed at stage boundaries, not continuously guarded. No unit tests.')
    save();return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('campaign')
    args=parser.parse_args();run(args.destination,args.campaign)
