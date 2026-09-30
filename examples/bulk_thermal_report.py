"""Audited bulk continuity/thermal convergence and mechanical work limitations."""
import argparse
import json
from pathlib import Path
import tarfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .bulk_thermal_records import sha
from .bulk_compression_reference import BulkCompressionReference
from .physical_scalar_sampler import PhysicalScalarSampler
from plasma_bench.mhd_mesh import MHDMesh


def checked(directory,name,count):
    records=json.loads((directory/name).read_text());audit=json.loads((directory/'audit.json').read_text())
    if len(records)!=count or len(audit['cases'])!=count or sha(directory/name)!=audit['records_sha256']:
        raise RuntimeError('Complete unchanged audited protocol required')
    for r in records:
        if r['status']!='complete' or sha(directory/r['arrays_file'])!=r['arrays_sha256']:raise RuntimeError('Raw record changed')
    return records,audit


def time_difference(directory,records):
    base=records[1];half=records[6];ref=BulkCompressionReference(**base['reference']);nr=base['spec']['nr']
    grid=PhysicalScalarSampler(MHDMesh(ref.radius,ref.half_length,nr,2*nr),5,18)
    with np.load(directory/base['arrays_file']) as a,np.load(directory/half['arrays_file']) as b:
        na=grid.density(a['density_final']);nb=grid.density(b['density_final'])
        result=dict(density=grid.metric(na,nb))
        for label,gamma in (('electron',ref.gamma_e),('ion',ref.gamma_i)):
            sa=grid.scalar(a[label+'_final']);sb=grid.scalar(b[label+'_final'])
            background=3+.1*(grid.r/ref.radius)**2 if label=='electron' else 2.5-.15*(grid.r/ref.radius)**2
            result[label]=dict(entropy_wave=grid.metric(sa-background,sb-background),
                temperature=grid.metric(np.exp(sa)*(na/1e20)**(gamma-1),np.exp(sb)*(nb/1e20)**(gamma-1)),
                work_defect_ratio=base['conservation'][label]['work_max_relative_defect']/half['conservation'][label]['work_max_relative_defect'])
    return result


def report(rate_directory,evolution_directory,probe_directory):
    rates,ra=checked(rate_directory,'rates.json',14);runs,ta=checked(evolution_directory,'transports.json',10)
    probe=json.loads((probe_directory/'probe.json').read_text())
    pa=json.loads((probe_directory/'audit.json').read_text())
    if (len(probe['cases'])!=6 or pa['cases']!=6 or not pa['all_arrays_bitwise_replayed']
            or pa['record_sha256']!=sha(probe_directory/'probe.json')
            or pa['source_sha256']!=ra['source_sha256']):
        raise RuntimeError('Complete unchanged audited pressure probes required')
    if ra['source_sha256']!=ta['source_sha256'] or probe['calculation_source_sha256']!=ra['source_sha256']:
        raise RuntimeError('Mixed calculation sources')
    for row in probe['cases']:
        if sha(probe_directory/row['arrays_file'])!=row['arrays_sha256']:raise RuntimeError('Changed pressure-work probe')
    maxima=dict(scalar_residual=max(w['residuals'][k] for r in rates for w in r['work'].values() for k in ('entropy_mass','thermal_projection')),
        pressure_work_backward=max(w['backward'] for r in rates for w in r['work'].values()),
        independent_local_energy=max(v for r in ra['cases'] for v in r['independent_local_energy'].values()),
        physical_entropy_rate=max(b['entropy_backward'] for r in ra['cases'] for b in r['norms']['balances'].values()),
        weak_continuity=max(b['weak_continuity_backward'] for r in ra['cases'] for b in r['norms']['balances'].values()),
        transport_particle_drift=max(r['conservation']['particles'] for r in runs),
        transport_entropy_drift=max(r['conservation'][label]['entropy_max_relative_drift'] for r in runs for label in ('electron','ion')),
        transport_work_defect=max(r['conservation'][label]['work_max_relative_defect'] for r in runs for label in ('electron','ion')),
        independent_snapshot_integrals=max(v for r in ta['cases'] for row in r['independent_integral_agreement'] for v in row.values()))
    time=time_difference(evolution_directory,runs);main=rates[4:8];nr=[8,16,32,64]
    fig,axes=plt.subplots(3,2,figsize=(12,12),constrained_layout=True)
    for key,label in (('density_rate','density rate'),('particle_flux','particle flux')):
        axes[0,0].loglog(nr,[r['norms']['metrics'][key]['relative_L2'] for r in main],'-o',label=label)
    for species,style in (('electron','-o'),('ion','--s')):
        axes[0,1].loglog(nr,[r['norms']['metrics'][species]['entropy_rate']['axis_relative_L2'] for r in main],style,label=species)
    for selected,label,style in ((runs[:3],'compression','-o'),(runs[3:6],'full transit','--s')):
        axes[1,0].loglog([8,16,32],[r['norms'][-1]['density']['relative_L2'] for r in selected],style,label=label)
        for species,color in (('electron','#007d8c'),('ion','#b34732')):
            axes[1,1].loglog([8,16,32],[r['norms'][-1][species]['temperature']['relative_L2'] for r in selected],style,color=color,label=species+' '+label)
    for species,style in (('electron','-o'),('ion','--s')):
        axes[2,0].loglog([8,16,32],[r['conservation'][species]['work_max_relative_defect'] for r in runs[:3]],style,label=species)
    for key,label,style in (('conservative_force_relative_L2','pressure-force error','-o'),('work_backward','mechanical/thermal work defect','--s')):
        axes[2,1].loglog(nr,[r[key] for r in probe['cases'][:4]],style,label=label)
    titles=('Shared physical particle transport','Entropy rate in first physical annulus',
            'Evolved density','Evolved species temperature','Compression pressure-work ledger error',
            'Candidate mechanical pressure force (diagnostic)')
    for index,(ax,title) in enumerate(zip(axes.flat,titles)):
        ax.set_title(title,fontsize=11);ax.set_xlabel('Radial cells; axial cells = 2 × radial');ax.set_ylabel('Relative error / stated normalization')
        ticks=nr if index in (0,1,5) else [8,16,32]
        ax.set_xticks(ticks,[str(n) for n in ticks]);ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(True,which='both',alpha=.18);ax.legend(fontsize=8)
    fig.suptitle('Shared continuity and bulk entropy: physical accuracy, evolving density and work',fontsize=14)
    fig.savefig(rate_directory/'verification.png',dpi=160);plt.close(fig)
    helper=Path(__file__)
    with tarfile.open(rate_directory/'report_sources.tar.gz','w:gz') as archive:
        archive.add(helper,arcname=helper.name);archive.add(helper.with_suffix('.md'),arcname=helper.with_suffix('.md').name)
    output=dict(source_sha256=ra['source_sha256'],rate_cases=14,transport_cases=10,pressure_probe_cases=6,maxima=maxima,time_control=time,
        rate_summary=[{k:r[k] for k in ('spec','norms','work','independent_local_energy')} for r in rates],
        transport_summary=[{k:r[k] for k in ('spec','norms','conservation','dt','steps','end_time')} for r in runs],pressure_probe=probe,
        rate_audit_sha256=sha(rate_directory/'audit.json'),transport_audit_sha256=sha(evolution_directory/'audit.json'),probe_sha256=sha(probe_directory/'probe.json'),
        probe_audit_sha256=sha(probe_directory/'audit.json'),
        report_source_sha256=sha(helper),report_companion_sha256=sha(helper.with_suffix('.md')),figure_sha256=sha(rate_directory/'verification.png'),
        scope='Actual changing-density/two-entropy transport with prescribed bulk velocity; mechanical/current/field coupling remains incomplete')
    (rate_directory/'report.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(dict(maxima=maxima,time_control=time),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('rates',type=Path);parser.add_argument('evolution',type=Path);parser.add_argument('probe',type=Path)
    args=parser.parse_args();report(args.rates,args.evolution,args.probe)
