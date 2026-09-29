"""Audited physical thermal convergence, axis correction and transport report."""
import argparse
import json
from pathlib import Path
import tarfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.kernel_induction import difference_curl
from .physical_scalar_sampler import PhysicalScalarSampler
from .log_entropy_reference import LogEntropyReference
from .log_thermal_records import sha


def checked(directory,filename,count):
    records=json.loads((directory/filename).read_text());audit=json.loads((directory/'audit.json').read_text())
    if len(records)!=count or len(audit['cases'])!=count or sha(directory/filename)!=audit['records_sha256']:
        raise RuntimeError('Complete unchanged audited protocol required')
    for record in records:
        if record['status']!='complete' or sha(directory/record['arrays_file'])!=record['arrays_sha256']:
            raise RuntimeError('Incomplete or changed raw data')
    return records,audit


def full_field_conditioning(directory,records):
    results=[]
    for record in records:
        spec=record['spec'];nr=spec['nr'];ref=record['reference']
        mesh=MHDMesh(ref['radius'],ref['half_length'],nr,2*nr)
        geo=CompatibleGeometry(mesh,degree=spec.get('degree',5),quadrature_order=spec.get('quadrature_order',10))
        with np.load(directory/record['arrays_file']) as arrays:
            bt=arrays['joint_magnetic_rate'];difference=difference_curl(geo,arrays['joint_electric'])-bt
        squared=max(float(difference@(geo.hessian@difference)),0.)
        scale=max(float(bt@(geo.hessian@bt)),0.)
        results.append(dict(spec=spec,physical_faraday_rate_relative=float(np.sqrt(squared/scale)) if scale>0 else None,
                            physical_faraday_rms_t_s=float(np.sqrt(MU0*squared/mesh.volumes.sum()))))
    return results


def time_control(directory,records):
    base=records[1];fine=records[3];ref=LogEntropyReference(**base['reference'])
    nr=base['spec']['nr'];mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    grid=PhysicalScalarSampler(mesh,5,14)
    with np.load(directory/base['arrays_file']) as a,np.load(directory/fine['arrays_file']) as b:
        difference=grid.scalar(a['final']-b['final'])
    target=ref.fields(grid.r,grid.z,ref.transit_time)['entropy']-ref.background_entropy
    return dict(nr=nr,relative_to_entropy_wave=float(np.sqrt(np.sum(grid.weight*difference**2)/np.sum(grid.weight*target**2))),
                rms_difference=float(np.sqrt(np.sum(grid.weight*difference**2)/mesh.volumes.sum())))


def report(rate_directory,transport_directory):
    rates,rate_audit=checked(rate_directory,'rates.json',18)
    transports,transport_audit=checked(transport_directory,'transports.json',9)
    if rate_audit['source_sha256']!=transport_audit['source_sha256']:raise RuntimeError('Rate and transport source identities differ')
    conditioning=full_field_conditioning(rate_directory,rates);time=time_control(transport_directory,transports)
    maxima=dict(pressure_work=max(s['pressure_work_backward'] for r in rates for s in r['work'].values()),
        physical_entropy_rate=max(s['physical_entropy_backward'] for r in rates for s in r['work'].values()),
        total_energy_exchange=max(r['energy_exchange_backward'] for r in rates),
        scalar_residual=max(s['scalar_residuals'][key] for r in rates for s in r['work'].values() for key in ('entropy_mass','thermal_projection')),
        induction_residual=max(r['induction_residual'] for r in rates),
        original_backward=max(a['original_equations'][key] for a in rate_audit['cases'] for key in ('ohm','ampere','faraday_backward')),
        independent_local_energy=max(v for a in rate_audit['cases'] for values in a['independent_energy'].values() for v in values.values()),
        independent_quadrature=max(v for r in rates for values in r['independent_quadrature_agreement'].values() for v in values.values()),
        transport_scalar_residual=max(r['max_scalar_residual'] for r in transports),
        transport_entropy_drift=max(r['conservation']['entropy']['max_relative_drift'] for r in transports),
        transport_integral_agreement=max(v for a in transport_audit['cases'] for row in a['independent_integral_agreement'] for v in row.values()),
        full_faraday_physical_rate_relative=max(r['physical_faraday_rate_relative'] for r in conditioning if r['physical_faraday_rate_relative'] is not None),
        full_faraday_physical_rms_t_s=max(r['physical_faraday_rms_t_s'] for r in conditioning))
    nr=[8,16,32,64];main=rates[4:8];fig,axes=plt.subplots(3,2,figsize=(12,12),constrained_layout=True)
    for scheme,color,style in (('legacy','#b34732','--'),('physical','#007d8c','-')):
        axes[0,0].loglog(nr,[r['norms']['nodal']['electron'][scheme]['axis_four']['relative_L2'] for r in main],style+'o',color=color,label=scheme)
    for key,color,style in (('legacy_proxy_rate','#b34732','--'),('proxy_rate','#007d8c','-')):
        axes[0,1].loglog(nr,[r['norms']['metrics']['electron'][key]['relative_L2'] for r in main],style+'o',color=color,label=key.replace('_',' '))
    for species,color in (('electron','#007d8c'),('ion','#987123')):
        axes[1,0].loglog(nr,[r['norms']['metrics'][species]['electric']['relative_L2'] for r in main],'-o',color=color,label=species)
    for key,style in (('relative_L2','-o'),('axis_relative_L2','--s')):
        axes[1,1].loglog([8,16,32],[r['norms'][-1]['entropy_wave'][key] for r in transports[:3]],style,label='whole domain' if key=='relative_L2' else 'first annulus')
    for key,style in (('rms_error','-o'),('sampled_maximum','--s')):
        axes[2,0].loglog(nr,[r['norms']['metrics']['joint_magnetic_rate_absolute'][key] for r in rates[14:18]],style,label='RMS' if key=='rms_error' else 'sampled maximum')
    for key in ('entropy','proxy','internal'):
        axes[2,1].loglog([8,16,32],[r['conservation'][key]['max_relative_drift'] for r in transports[:3]],'-o',label=key)
    titles=('Near-axis nodal log-entropy rate','Continuous physical proxy rate',
            'Pressure electric field','Full-transit sheared entropy wave',
            'Residual isothermal magnetic battery','Full-transit conservation (absolute drift)')
    for index,(ax,title) in enumerate(zip(axes.flat,titles)):
        ax.set_title(title,fontsize=11);ax.set_xlabel('Radial cells; axial cells = 2 × radial')
        ax.set_ylabel('Relative error' if index!=4 else 'Spurious B rate (T/s)')
        ticks=[8,16,32] if index in (3,5) else nr
        ax.set_xticks(ticks,[str(n) for n in ticks]);ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(True,which='both',alpha=.18);ax.legend(fontsize=9)
    fig.suptitle('Physical entropy component: accuracy, work pairing and finite-time transport',fontsize=14)
    fig.savefig(rate_directory/'verification.png',dpi=160);plt.close(fig)
    helper=Path(__file__)
    with tarfile.open(rate_directory/'report_sources.tar.gz','w:gz') as archive:
        for path in (helper,helper.with_suffix('.md')):archive.add(path,arcname=path.name)
    output=dict(source_sha256=rate_audit['source_sha256'],rate_cases=18,transport_cases=9,maxima=maxima,
        rate_summary=[dict(spec=r['spec'],norms=r['norms'],work=r['work'],energy_exchange_backward=r['energy_exchange_backward']) for r in rates],
        transport_summary=[{k:r[k] for k in ('spec','norms','conservation','dt','steps','end_time')} for r in transports],
        full_field_conditioning=conditioning,time_control=time,
        rate_audit_sha256=sha(rate_directory/'audit.json'),transport_audit_sha256=sha(transport_directory/'audit.json'),
        report_source_sha256=sha(helper),report_companion_sha256=sha(helper.with_suffix('.md')),
        figure_sha256=sha(rate_directory/'verification.png'),
        scope='Physical current/thermal component and prescribed-current transport; full conservative moving state and vacuum remain uncoupled')
    (rate_directory/'report.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(dict(maxima=maxima,time_control=time),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('rates',type=Path);parser.add_argument('transports',type=Path)
    args=parser.parse_args();report(args.rates,args.transports)
