"""Independent physical kinetic metric and analytic-force obstruction at the axis."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh
from plasma_bench.lobatto_sbp import LobattoSBP
from .regular_euler_sampling import basis
from .radial_work_records import sha


def calculate(degree):
    q,w=np.polynomial.legendre.leggauss(degree+1);nodes=.5*(q+1);weights=.5*w*nodes
    qh,wh=np.polynomial.legendre.leggauss(2*degree+4);r=.5*(qh+1);volume=.5*wh*r
    old=r[:,None]*basis(nodes**2,r**2)/nodes[None,:]
    metric=old.T@(volume[:,None]*old);lumped=np.diag(weights);eigen=eigvalsh(metric,lumped)
    force=old.T@(volume*r)/weights
    old_error=float(np.sqrt(np.sum(weights*(force-nodes)**2)/np.sum(weights*nodes**2)))
    b=LobattoSBP(degree);x=.5*(b.nodes+1);wx=.5*b.weights
    # In x=r² coordinates the common constant factors cancel. The physical
    # transverse kinetic metric is integral x*u² dx over active positive nodes.
    xh=.5*(qh+1);vx=basis(x[1:],xh);new=vx.T@((.5*wh*xh)[:,None]*vx);new_lumped=np.diag((wx*x)[1:]);new_eigen=eigvalsh(new,new_lumped)
    acceleration=(vx.T@(.5*wh*xh))/(wx*x)[1:]
    new_error=float(np.sqrt(np.sum((wx*x)[1:]*(acceleration-1)**2)/np.sum((wx*x)[1:])))
    record=dict(degree=degree,gauss_r_regular_interpolation_energy_eigenvalue_range=[float(eigen[0]),float(eigen[-1])],
        gauss_r_consistent_force_lumped_inverse_error=old_error,gauss_r_point_force_ratios=(force/nodes).tolist(),
        regular_x_active_velocity_energy_eigenvalue_range=[float(new_eigen[0]),float(new_eigen[-1])],
        regular_x_linear_force_error=new_error)
    return record,dict(old_metric=metric,old_lumped_metric=lumped,old_force=force,old_nodes=nodes,new_metric=new,new_lumped_metric=new_lumped,new_regular_acceleration=acceleration,new_nodes=x)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory;d.mkdir(parents=True,exist_ok=False);records=[]
    for p in (3,5,7):
        record,arrays=calculate(p);path=d/f'degree_{p}.npz';np.savez_compressed(path,**arrays);record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record)
    (d/'records.json').write_text(json.dumps(dict(script_sha256=sha(__file__),records=records),indent=2)+'\n')
    (d/'probe_source.py').write_bytes(Path(__file__).read_bytes());(d/'probe_source.md').write_bytes(Path(__file__).with_suffix('.md').read_bytes())
    print(json.dumps(records,indent=2))


if __name__=='__main__':main()
