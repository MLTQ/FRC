"""Scientific physical-rate refinement figure, including axis and density controls."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(records,path):
    lookup={r['configuration']['name']:r for r in records if r['completed']}
    fig,axes=plt.subplots(2,3,figsize=(12,7),sharex=True)
    for column,(key,title) in enumerate((('magnetic_rate','Magnetic rate'),('current_rate','Current rate'),('electric','Electric field'))):
        for row,(metric,region) in enumerate((('relative_L2','Whole volume'),('axis_relative_L2','First annulus'))):
            ax=axes[row,column]
            for prefix,label,color in (('moving','Ar+, n₀ = 10¹⁸ m⁻³','#2067a5'),('target','Ar+, n₀ = 10²⁰ m⁻³','#cc6b20'),('massless','Massless control, n₀ = 10¹⁸','#6b8e38')):
                group=[lookup[f'{prefix}_{n}'] for n in (8,16,32,64,128) if f'{prefix}_{n}' in lookup]
                ax.loglog([r['configuration']['nr'] for r in group],[r['metrics'][key][metric] for r in group],'-o',label=label,color=color,markersize=4)
            ax.set_title(f'{title} · {region}',fontsize=11);ax.grid(True,which='both',alpha=.2)
            ax.set_xticks([8,16,32,64,128],labels=['8','16','32','64','128'])
            if column==0:ax.set_ylabel('Relative L² error')
            if row==1:ax.set_xlabel('Radial cells (axial cells = 2×)')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Moving finite-mass induction: independently forced continuum reference',fontsize=14)
    fig.text(.5,.94,'Instantaneous rates; no time evolution, free plasma edge or formation claim',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.06,1,.92));fig.savefig(path,dpi=170);plt.close(fig)
