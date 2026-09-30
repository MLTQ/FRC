"""Scientific report from saved physical coupling evidence."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from .radial_work_records import sha


def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    r=[json.loads(f.read_text()) for f in sorted(a.study.glob('case_*.json'))];audit=json.loads((a.study/'audit.json').read_text())
    if len(r)!=16 or len(audit['records'])!=16:raise ValueError('Complete study and audit required')
    colors=['#157f91','#db6944','#7b5bb7'];names=['Radial','Azimuthal','Axial'];ne=[2,4,8,16]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,3,figsize=(15,8.8),layout='constrained')
    for k,(color,name) in enumerate(zip(colors,names)):
        ax[0,0].loglog(ne,[r[i]['force_relative_L2'][k] for i in range(4,8)],'o-',color=color,label=name)
        ax[0,1].loglog(ne,[r[i]['axis_force_relative_L2'][k] for i in range(4,8)],'o-',color=color,label=name)
        ax[0,2].loglog(ne,[r[i]['lorentz_force_error'][k] for i in range(4,8)],'o-',color=color,label=name)
    titles=['Independent smooth force, degree 5','First element next to the axis','Actual Ampère current → Lorentz force']
    for axis,title in zip(ax[0],titles):
        axis.set(title=title,xlabel='Fluid radial elements',ylabel='Relative L2 force-rate error');axis.set_xticks(ne,labels=ne)
    ax[0,0].legend(frameon=False)
    cells=[8,16,32]
    for key,color,label in [('current_error',colors[0],'Current'),('magnetic_rate_error',colors[1],'Magnetic rate'),('current_rate_error',colors[2],'Current rate')]:
        ax[1,0].loglog(cells,[r[i][key] for i in (5,6,7)],'o-',color=color,label=label)
    ax[1,0].set(title='Original finite-inertia induction',xlabel='Field radial cells',ylabel='Physical relative L2 error');ax[1,0].set_xticks(cells,labels=cells);ax[1,0].legend(frameon=False)
    ax[1,1].loglog(cells,[r[i]['magnetic_rate_error'] for i in (13,6,14)],'o-',color=colors[1],label='Magnetic rate: 8 fluid elements')
    ax[1,1].loglog(cells,[r[i]['continuity_error'] for i in (13,6,14)],'s--',color=colors[0],label='Same actual continuity rate')
    ax[1,1].set(title='Independent field-grid control',xlabel='Field radial cells',ylabel='Physical relative L2 error');ax[1,1].set_xticks(cells,labels=cells);ax[1,1].legend(frameon=False,fontsize=8)
    ax[1,2].plot(ne,[1000*r[i]['incomplete_total_energy_rate'] for i in range(4,8)],'o-',color=colors[1],label='Fluid + magnetic + current rate')
    ax[1,2].plot(ne,[1000*r[i]['unresolved_density_exchange'] for i in range(4,8)],'x--',color=colors[0],label='Unclosed density exchange')
    ax[1,2].set(title='Missing transport remains explicit',xlabel='Fluid radial elements',ylabel='Energy rate (mW)');ax[1,2].set_xticks(ne);ax[1,2].ticklabel_format(axis='y',style='plain',useOffset=False);ax[1,2].legend(frameon=False,fontsize=8)
    for axis in ax.ravel():
        axis.xaxis.set_minor_locator(NullLocator());axis.grid(alpha=.18)
    fig.suptitle('Connecting the regular fluid to physical current and magnetic work',fontsize=18,weight='bold')
    fig.supxlabel('Instantaneous component verification • all 16 cases independently audited • full moving-species transport still open',fontsize=10)
    fig.savefig(a.output/'report.png',dpi=180);fig.savefig(a.output/'report.pdf');plt.close(fig)
    (a.output/'report_source.json').write_text(json.dumps(dict(source_sha256=audit['source_sha256'],study=str(a.study),audit_sha256=sha(a.study/'audit.json'),script_sha256=sha(__file__)),indent=2)+'\n')


if __name__=='__main__':main()
