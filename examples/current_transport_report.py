"""Compact measured report from immutable transport study and audit records."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def main():
    root=Path(__file__).resolve().parents[1];src=root/'outputs/current_transport_verification';out=root/'outputs/current_transport_diagnostics'
    records=[json.loads((src/f'case_{i:02d}.json').read_text()) for i in range(16)]
    audit=json.loads((src/'audit.json').read_text());assert audit['terminal'] and len(audit['records'])==16
    cases=audit['records'];rows=[]
    for i in range(4,8):
        r=records[i];a=cases[i]['independent'];rows.append(dict(elements=r['spec']['elements'],field_cells=2*r['spec']['elements'],electric=r['electric_relative_L2'],force={k:a['norms']['stress_rate_'+k]['relative_L2'] for k in ('radial','azimuthal','axial')},axis_force={k:a['norms']['stress_rate_'+k]['first_fluid_element'] for k in ('radial','azimuthal','axial')},wall_force={k:a['norms']['stress_rate_'+k]['last_fluid_element'] for k in ('radial','azimuthal','axial')},grid_power=r['total_energy_rate'],continuous_power=a['powers']['continuous_total_rate']))
    summary=dict(source_sha256=audit['source_sha256'],cases=len(cases),degree5=rows,
        max_grid_energy_backward=max(r['total_energy_backward'] for r in records),max_continuous_energy_backward=max(a['independent']['powers']['continuous_total_backward'] for a in cases),
        max_original_residual=max(r['residual'] for r in records),
        max_work_backward={k:max(r[k] for r in records) for k in ('advection_work_backward','stress_work_backward','self_drift_work_backward','axial_force_scaled','angular_torque_scaled')},
        max_power_difference_scaled=max(v['scaled_to_physical_rate'] for a in cases for v in a['production_power_difference'].values()),
        max_power_quadrature_scaled=max(v['absolute']/a['independent']['powers']['continuous_rate_scale'] for a in cases for k,v in a['quadrature_power'].items() if k not in ('continuous_total_backward','continuous_rate_scale')))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'figure.dpi':150})
    fig,axs=plt.subplots(2,3,figsize=(15,8.3),layout='constrained');x=np.array([r['elements'] for r in rows]);colors=['#146c94','#cf4d36','#26846b']
    for k,c in zip(('advection','stretching','self_drift'),colors):axs[0,0].loglog(x,[r['electric'][k] for r in rows],'o-',label=k.replace('_',' '),color=c)
    for k,c in zip(('radial','azimuthal','axial'),colors):
        axs[0,1].loglog(x,[r['force'][k] for r in rows],'o-',label=k,color=c)
        axs[0,2].loglog(x,[r['axis_force'][k] for r in rows],'o-',label=k+' axis',color=c)
        axs[0,2].loglog(x,[r['wall_force'][k] for r in rows],'s--',label=k+' wall',color=c,alpha=.7)
    for k,label in [('continuity_point_defect','continuity defect'),('charge_divergence_scaled','scaled div J')]:axs[1,0].loglog(x,[records[i][k] for i in range(4,8)],'o-',label=label)
    axs[1,0].loglog(x,[records[i]['normal_current_jumps'][1] for i in range(4,8)],'o-',label='axial current jumps')
    axs[1,1].loglog(x,[cases[i]['independent']['norms']['omitted_stress_rate']['relative_L2']/cases[i]['independent']['norms']['stress_rate']['relative_L2'] for i in range(4,8)],'o-',label='drift force')
    axs[1,1].loglog(x,[records[i]['advection_without_faces']/records[i]['electric_relative_L2']['advection'] for i in range(4,8)],'o-',label='current advection')
    axs[1,2].semilogy(x,[abs(r['grid_power']) for r in rows],'o-',label='degree 5 grid total')
    axs[1,2].semilogy(x,[abs(r['continuous_power']) for r in rows],'s--',label='degree 5 continuous total')
    axs[1,2].semilogy(x,[abs(cases[i]['independent']['powers']['continuous_total_rate']) for i in range(4)],'d-',label='degree 3 continuous total',color='#94603e')
    titles=['Original electric transport · relative L2','Physical drift force · relative L2','Force accuracy at axis and wall','Reconstructed continuity defects','Error increase when faces are omitted','Total energy rate magnitude (W)']
    for ax,title in zip(axs.ravel(),titles):ax.set_title(title);ax.set_xlabel('Radial fluid elements (degree 5)');ax.set_xticks(x,labels=x);ax.xaxis.set_minor_formatter(NullFormatter());ax.legend(fontsize=8)
    axs[1,2].set_xlabel('Radial fluid elements')
    fig.suptitle('Current transport verification | 16 frozen cases + independent physical audit\nOne instantaneous filled-cylinder fixture; species pressure and coupled evolution remain open',fontsize=14)
    fig.savefig(out/'convergence.png');fig.savefig(out/'convergence.pdf');plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
