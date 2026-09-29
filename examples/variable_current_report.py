"""Assemble audited variable-density measurements and a scientific comparison plot."""
import argparse
import json
from pathlib import Path
import tarfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .variable_current_records import sha
from .variable_current_reference import VariableCurrentReference
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_current import DualCurrent


def thermal_regions(record,arrays):
    ref=VariableCurrentReference(**record['reference']);nr=record['spec']['nr']
    mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr);dual=DualCurrent(mesh);result={}
    for label,gamma in (('electron',5/3),('ion',1.4)):
        exact=ref.thermal(dual.r,dual.z,gamma=gamma)['entropy_rate']
        error=arrays[label+'_entropy_rate']-exact
        weighted=dual.volumes*error**2;total=weighted.sum()
        regions={}
        for name,rows in (('axis_four',slice(0,4)),('interior',slice(4,-4)),
                          ('wall_four',slice(-4,None))):
            w=dual.volumes[rows];ee=error[rows];target=exact[rows]
            regions[name]=dict(squared_error_fraction=float(weighted[rows].sum()/total),
                relative_L2=float(np.sqrt(np.sum(w*ee**2)/np.sum(w*target**2))),
                rms_error=float(np.sqrt(np.sum(w*ee**2)/np.sum(w))),
                sampled_maximum=float(abs(ee).max()))
        result[label]=dict(regions=regions,radial_squared_error=weighted.sum(axis=1).tolist(),
                           peak_radial_row=int(np.argmax(weighted.sum(axis=1))))
    return result


def report(directory):
    records=json.loads((directory/'records.json').read_text())
    audit=json.loads((directory/'audit.json').read_text())
    if len(records)!=14 or len(audit['cases'])!=14:raise RuntimeError('Complete 14-case protocol required')
    if sha(directory/'records.json')!=audit['records_sha256']:raise RuntimeError('Audited records changed')
    for record in records:
        if record['status']!='complete' or sha(directory/record['arrays_file'])!=record['arrays_sha256']:
            raise RuntimeError('Incomplete or changed raw arrays')
    def metric(row,key):return row['norms']['metrics'][key]['relative_L2']
    rows=[]
    for record in records:
        with np.load(directory/record['arrays_file']) as arrays:regions=thermal_regions(record,arrays)
        rows.append(dict(spec=record['spec'],metrics=record['norms']['metrics'],pressure=record['pressure'],
            residual=record['residual'],iterations=record['iterations'],energy_exchange_backward=record['energy_exchange_backward'],
            thermal_error_regions=regions))
    maxima=dict(original_backward=max(v for a in audit['cases'] for k,v in a['original_equations'].items() if k!='faraday_rate_relative'),
        full_faraday_rate_relative=max(a['original_equations']['faraday_rate_relative'] for a in audit['cases']),
        independent_local_energy=max(v for a in audit['cases'] for v in a['independent_energy_relative_difference'].values()),
        density_pullback=max(a['density_pullback_backward'] for a in audit['cases']),
        production_quadrature=max(v for r in records for v in r['production_energy_relative_difference'].values()),
        independent_quadrature=max(v for r in records for v in r['independent_quadrature_relative_difference'].values()),
        solver_residual=max(r['residual'] for r in records),
        energy_exchange_backward=max(r['energy_exchange_backward'] for r in records),
        continuity_relative=max(r['continuity_relative_L2'] for r in records),
        continuity_integral_backward=max(r['continuity_integral_backward'] for r in records),
        pressure_work_backward=max(s['work_backward'] for r in records for s in r['pressure'].values()),
        log_entropy_backward=max(s['entropy_backward'] for r in records for s in r['pressure'].values()))
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    nr=[8,16,32,64]
    for degree,color,style in ((3,'#8f6b32','--'),(5,'#007d8c','-')):
        group=records[:4] if degree==3 else records[4:8]
        axes[0,0].loglog(nr,[metric(r,'magnetic_rate') for r in group],style+'o',color=color,label=f'B rate, degree {degree}')
        axes[0,0].loglog(nr,[metric(r,'current_rate') for r in group],style+'s',color=color,alpha=.7,label=f'J rate, degree {degree}')
    group=records[4:8]
    for key,label in (('current_energy_rate','Current-energy rate'),('density_energy_rate','Density contribution'),('magnetic_energy_rate','Magnetic-energy rate')):
        axes[0,1].loglog(nr,[metric(r,key) for r in group],'-o',label=label)
    for key,label in (('electron_electric','Electron pressure E'),('ion_electric','Ion pressure E')):
        axes[1,0].loglog(nr,[metric(r,key) for r in group],'-o',label=label)
    axes[1,0].loglog(nr,[r['pressure']['electron']['nodal_entropy_rate_relative_L2'] for r in group],'-s',label='Electron entropy rate (nodal)')
    axes[1,0].loglog(nr,[r['pressure']['ion']['nodal_entropy_rate_relative_L2'] for r in group],':^',label='Ion entropy rate (nodal)')
    axes[1,1].loglog(nr,[r['pressure']['electron']['wrong_old_stream_work_relative'] for r in group],'-o',color='#b24032',label='Old stream: relative net-work mismatch')
    axes[1,1].loglog(nr,[r['pressure']['electron']['work_backward'] for r in group],'-s',color='#007d8c',label='New pairing: work backward error')
    titles=('Variable-density original induction','Independent local energy rates (degree 5)',
            'Pressure bridge retains accuracy limits','Different work-error scales; neither proves accuracy')
    for ax,title in zip(axes.flat,titles):
        ax.set_title(title,fontsize=11);ax.set_xlabel('Radial cells (axial cells = 2 × radial)')
        ax.set_ylabel('Measured dimensionless error');ax.set_xticks(nr,[str(n) for n in nr])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(True,which='both',alpha=.18);ax.legend(fontsize=8)
    fig.suptitle('Changing-density component verification — not a formation trajectory',fontsize=14)
    fig.savefig(directory/'verification.png',dpi=170);plt.close(fig)
    helper=Path(__file__);source_hashes={helper.name:sha(helper),helper.with_suffix('.md').name:sha(helper.with_suffix('.md'))}
    with tarfile.open(directory/'report_sources.tar.gz','w:gz') as archive:
        for path in (helper,helper.with_suffix('.md')):archive.add(path,arcname=path.name)
    result=dict(source_sha256=audit['source_sha256'],cases=14,rows=rows,maxima=maxima,
        records_sha256=sha(directory/'records.json'),audit_sha256=sha(directory/'audit.json'),
        figure_sha256=sha(directory/'verification.png'),report_source_files=source_hashes,
        scope='Instantaneous manufactured component and explicit lumped pressure bridge; no full moving evolution or formation adoption')
    (directory/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(maxima,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    report(parser.parse_args().directory)
