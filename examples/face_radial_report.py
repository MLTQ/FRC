"""Scientific convergence and conservation figure from preserved numerical records."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root=Path(__file__).resolve().parents[1];out=root/'outputs';target=out/'face_radial_diagnostics';target.mkdir(exist_ok=True)
    rates=json.loads((out/'face_radial_rates/rates.json').read_text())
    old=json.loads((out/'radial_work_rates/rate_audit.json').read_text())['records']
    spectrum=json.loads((out/'face_radial_spectrum/records.json').read_text())['records']
    groups={p:[r for r in rates if r['spec'].get('degree')==p] for p in (3,5)}
    fine=groups[5];nr=[r['spec']['nr'] for r in fine]
    plt.rcParams.update({'font.size':11,'axes.titlesize':13,'axes.labelsize':11,'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(17,10));fig.subplots_adjust(left=.065,right=.985,top=.84,bottom=.095,wspace=.30,hspace=.48)
    fig.suptitle('A common radial momentum field fixes the near-axis force',x=.065,y=.965,ha='left',fontsize=23,weight='bold')
    fig.text(.065,.905,'Original fluid equations • independent physical quadrature • preserved failed alternatives',fontsize=13,color='#44546a')
    a=axes[0,0]
    a.loglog(nr,[r['velocity_components']['radial']['axis_relative_L2'] for r in old[4:8]],'s--',color='#bd4936',label='Previous cell-radial work pairing')
    a.loglog(nr,[r['norms']['velocity_components']['radial']['axis_relative_L2'] for r in fine],'o-',color='#126d9c',label='Common radial-face field')
    a.set(title='Radial acceleration in the first annulus',ylabel='Relative L2 error');a.legend()
    a=axes[0,1]
    for key,label,color in [('pressure','Pressure response','#126d9c'),('centrifugal','Centrifugal response','#be8530'),('advection','Advection','#2b8d74')]:
        a.loglog(nr,[r['radial_components'][key]['axis_relative_L2'] for r in fine],'o-',label=label,color=color)
    a.set(title='Individual near-axis force contributions',ylabel='Relative L2 error');a.legend()
    a=axes[0,2]
    for p,color in [(3,'#8a69a1'),(5,'#126d9c')]:
        a.loglog(nr,[r['norms']['velocity']['relative_L2'] for r in groups[p]],'o-',label=f'Degree {p}',color=color)
    a.set(title='Full two-dimensional acceleration',ylabel='Relative L2 error');a.legend()
    a=axes[1,0]
    for p,color in [(3,'#8a69a1'),(5,'#126d9c')]:
        a.loglog(nr,[abs(r['scaled_global_energy_defect']) for r in groups[p]],'o-',label=f'Degree {p}',color=color)
    a.set(title='Physical energy-rate defect is retained',ylabel='Absolute defect / (initial energy × c/R)');a.legend()
    a=axes[1,1];selected=[r for r in spectrum if r['spec'].get('degree')==5 and 'quadrature_order' not in r['spec']]
    for k,color in enumerate(('#126d9c','#2b8d74','#be8530')):
        a.loglog([r['spec']['nr'] for r in selected],[max(abs(r['first_three_frequency_relative_errors'][k]),1e-16) for r in selected],'o-',label=f'Radial mode {k+1}',color=color)
    a.set(title='Sound frequencies vs analytic Bessel modes',ylabel='Absolute relative frequency error');a.legend()
    a=axes[1,2]
    with np.load(out/'face_radial_pilot/fluid_00.npz') as data:h=data['history']
    x=h[:,0]/h[-1,0];energy0=h[0,5]
    for y,label,color,style in [((h[:,5]-energy0)/energy0,'Physical energy','#126d9c','-'),(h[:,12]/energy0,'Integrated spatial defect','#be8530','--'),((h[:,4]-h[0,4])/energy0,'Energy account','#2b8d74',':')]:
        a.plot(x,y*1e9,label=label,color=color,linestyle=style,linewidth=2)
    a.set(title='Full contact transit: separate energy budgets',xlabel='Fraction of one transit',ylabel='Relative change × 10⁹');a.legend()
    for i,a in enumerate(axes.flat):
        a.grid(True,which='major',alpha=.18)
        if i<5:
            a.set_xlabel('Radial cells');a.set_xticks(nr if i!=4 else [8,16,32])
            a.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    fig.text(.065,.025,'Smooth unmagnetized verification only. Convective energy is not exact; moving magnetic/current coupling, heated edge and FRC formation remain open.',fontsize=10,color='#44546a')
    fig.savefig(target/'report.png',dpi=190);plt.close(fig)


if __name__=='__main__':main()
