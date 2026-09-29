"""Audited thermal improvement with the unresolved magnetic error kept visible."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def plot(report,directory):
    rates={r['configuration']['name']:r for r in report['rates']}
    runs={r['configuration']['name']:r for r in report['trajectories']}
    old_rates={r['configuration']['name']:r for r in report['baseline_rates']}
    old_runs={r['configuration']['name']:r for r in report['baseline_trajectories']}
    fig,axes=plt.subplots(2,2,figsize=(11.7,8.5),layout='constrained')
    ax=axes[0,0];nr=(8,16,32,64)
    for rows,label,style in ((old_rates,'Previous limited transport','o--'),(rates,'Shared smooth transport','o-')):
        ax.loglog(nr,[rows[f'wave_0_{n}']['metrics']['ion_pressure']['relative_L2'] for n in nr],style,label=label)
    ax.axhline(1,color='.5',ls=':');ax.set_title('Slow-wave ion-pressure rate improves')
    ax.set_ylabel('Error / physical pressure-rate signal');ax.legend(fontsize=8)
    ax=axes[0,1];controls=report['slow_mode_diagnosis']['cases']
    for key,label,style in (('both_numerical','Actual velocity/current representation','o-'),
                             ('exact_velocity','Exact velocity control','s--'),('both_exact','Both exact: oracle control','^--')):
        ax.loglog([r['nr'] for r in controls],[r['magnetic_rate_relative_L2'][key] for r in controls],style,label=label)
    ax.axhline(1,color='.5',ls=':');ax.set_title('Slow-wave magnetic rate still fails')
    ax.set_ylabel('Error / physical magnetic-rate signal');ax.legend(fontsize=8)
    ax=axes[1,0];nr=(8,16,32)
    for rows,label,style in ((old_runs,'Previous ion pressure','o--'),(runs,'Smooth ion pressure','o-')):
        selected=[n for n in nr if rows[f'wave_{n}']['completed']]
        ax.loglog(selected,[rows[f'wave_{n}']['metrics']['ion_pressure']['relative_L2'] for n in selected],style,label=label)
    ax.loglog(nr,[runs[f'wave_{n}']['metrics']['magnetic']['relative_L2'] for n in nr],'s-',label='Magnetic field (unchanged)')
    ax.set_title('Fast wave after one quarter period');ax.set_ylabel('Error / oscillation amplitude');ax.legend(fontsize=8)
    ax=axes[1,1];nz=[n for n in (16,32,64) if runs[f'contact_{n}']['completed']]
    for key,label in (('electron_pressure','Nodal electron pressure'),('ion_pressure','Cell-mean ion pressure')):
        ax.loglog(nz,[runs[f'contact_{n}']['metrics'][key]['relative_L2'] for n in nz],'o-',label=label)
    ax.loglog(nz,[runs[f'contact_{n}']['initial_metrics']['ion_pressure']['relative_L2'] for n in nz],'k:',label='Initial ion projection error')
    limited=runs['contact_limited_32']
    if limited['completed']:ax.plot([32],[limited['metrics']['electron_pressure']['relative_L2']],'rx',ms=9,label='Limited entropy comparison')
    failed=[name for name in runs if name.startswith('contact_') and not runs[name]['completed']]
    if failed:ax.text(.03,.03,'Partial runs excluded: '+', '.join(failed),transform=ax.transAxes,fontsize=7,wrap=True)
    ax.set_title('Unmagnetized contact after one transit');ax.set_ylabel('Error / full species pressure');ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(True,which='both',alpha=.2);ax.set_xlabel('Radial cells (nz = 2 nr)')
        ax.set_xticks([8,16,32,64],labels=['8','16','32','64'])
        ax.xaxis.set_minor_formatter(NullFormatter())
    axes[1,1].set_xlabel('Axial cells (nr = 4)')
    fig.suptitle('Shared transport improves temperature accuracy; magnetic cancellation remains unresolved',fontsize=12)
    path=directory/'accuracy.png';fig.savefig(path,dpi=170);plt.close(fig);return path
