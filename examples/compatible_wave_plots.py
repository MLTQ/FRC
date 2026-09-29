"""Scientific comparison of resolved spatial error and retained bias roundoff."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def plot(report,directory):
    rows={r['configuration']['name']:r for r in report['cases']}
    raw={r['configuration']['name']:r for r in report['rounded_cases']}
    nr=(8,16,32,64);fig,axes=plt.subplots(2,2,figsize=(11.5,8),layout='constrained')
    ax=axes[0,0]
    old=report['previous_linear_diagnosis']['cases']
    ax.loglog(nr,[r['magnetic_rate_relative_L2']['both_numerical'] for r in old],'o--',label='Previous linear representation')
    ax.loglog(nr,[rows[f'matched_linear_{n}']['metrics']['magnetic_rate']['relative_L2'] for n in nr],'o-',label='Compatible degree 5, matched input')
    ax.set_title('Slow mode: matched linear comparison');ax.set_ylabel('Error / magnetic-rate signal');ax.legend(fontsize=8)
    ax=axes[0,1]
    for collection,label,style in ((raw,'Rounded total field','s--'),(rows,'Explicit harmonic + increment','o-')):
        ax.loglog(nr,[collection[f'branch0_{n}']['metrics']['magnetic_rate']['relative_L2'] for n in nr],style,label=label)
    ax.plot([32],[raw['weak_full_bias']['metrics']['magnetic_rate']['relative_L2']],'rx',ms=9,label='Rounded: 100× weaker signal')
    ax.set_title('Same physical bias: precision matters');ax.set_ylabel('Error / magnetic-rate signal');ax.legend(fontsize=8)
    ax=axes[1,0]
    for branch,label in ((0,'Slow'),(1,'Intermediate'),(2,'Fast')):
        ax.loglog(nr,[rows[f'branch{branch}_{n}']['metrics']['magnetic_rate']['relative_L2'] for n in nr],'o-',label=label)
    ax.set_title('All branches: full physical bias');ax.set_ylabel('Error / magnetic-rate signal');ax.legend(fontsize=8)
    ax=axes[1,1]
    for key,label in (('magnetic','Initial magnetic field'),('current','Initial Ampere current'),('current_rate','Current rate')):
        ax.loglog(nr,[rows[f'branch0_{n}']['metrics'][key]['relative_L2'] for n in nr],'o-',label=label)
    ax.set_title('Slow mode: separate physical field norms');ax.set_ylabel('Error / physical signal');ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(True,which='both',alpha=.2);ax.set_xlabel('Radial cells (nz = 2 nr)')
        ax.set_xticks(nr,labels=[str(n) for n in nr]);ax.xaxis.set_minor_formatter(NullFormatter())
    fig.suptitle('Compatible field verification — linear components, not nonlinear formation',fontsize=12)
    path=directory/'accuracy.png';fig.savefig(path,dpi=170);plt.close(fig);return path
