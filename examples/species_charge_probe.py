"""Measure species pressure/entropy compatibility on saved verified currents."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E
from .regular_coupling_reference import initialize
from .current_transport_sampling import TransportSampler


def entropy(r,z,ref,k):
    y=(r/ref.radius)**2;phase=np.pi*z/ref.half_length
    return np.log(2. if k==0 else 1.5)+(0.13 if k==0 else -.09)*y*np.sin(phase)+.07*np.cos(2*phase)


def measure(spec,arrays,order):
    ref,d,u,m,_=initialize(spec);u=arrays['state'].copy()
    for k in range(2):u[...,k+4]=u[...,0]*entropy(d.r,d.z,ref,k)
    s=TransportSampler(d,m.mesh,order);n=s.scalar(u[...,0]);j=s.electric(arrays['current']);weight=s.weight
    rn,zn=m.mesh.corners;rn=rn[:,:-1];zn=zn[:,:-1];nn=ref.fields(rn,zn)['n'];result=[]
    for k,(gamma,a) in enumerate(zip(d.thermal.gammas,(-m.beta,m.alpha))):
        c=s.scalar(u[...,k+4]);se=c/n;h=E*np.exp(se)*(n/d.thermal.n0)**(gamma-1)/(gamma-1)
        gs=np.zeros_like(j);gh=np.zeros_like(j)
        for axis in (0,2):
            dn=s.scalar(u[...,0],axis=axis);ds=(s.scalar(u[...,k+4],axis=axis)-se*dn)/n
            gs[...,axis]=ds;gh[...,axis]=h*((gamma-1)*dn/n+ds)
        hn=E*np.exp(entropy(rn,zn,ref,k))*(nn/d.thermal.n0)**(gamma-1)/(gamma-1)
        gradient=s.electric(m.geometry.gradient(gamma*hn));entropy_gradient=s.electric(m.geometry.gradient(entropy(rn,zn,ref,k)))
        pressure=a/E*(gamma*gh-h[...,None]*gs);drift_entropy_rate=-a/E*np.sum(j*gs,axis=-1)
        thermal_work=float(np.sum(weight*h*drift_entropy_rate));pressure_work=float(np.sum(weight*np.sum(j*pressure,axis=-1)))
        compatible_pressure=a/E*(gradient-h[...,None]*gs)
        physical_gradient_work=float(a/E*np.sum(weight*np.sum(j*(gamma*gh),axis=-1)))
        compatible_gradient_work=float(a/E*np.sum(weight*np.sum(j*gradient,axis=-1)))
        entropy_total=float(np.sum(weight*drift_entropy_rate));compatible_entropy_total=float(-a/E*np.sum(weight*np.sum(j*entropy_gradient,axis=-1)))
        denominator=float(np.sum(weight*np.sum(pressure*pressure,axis=-1)))
        result.append(dict(species=('electron','ion')[k],mass_fraction=a,pressure_work=pressure_work,advective_thermal_work=thermal_work,
            unpaired_pressure_work=pressure_work-thermal_work,physical_gradient_work=physical_gradient_work,
            compatible_gradient_work=compatible_gradient_work,compatible_pressure_minus_thermal_work=float(np.sum(weight*np.sum(j*compatible_pressure,axis=-1)))-thermal_work,
            potential_reconstruction_relative=float(np.sqrt(np.sum(weight*np.sum((compatible_pressure-pressure)**2,axis=-1))/denominator)),
            advective_entropy_total=entropy_total,compatible_gradient_entropy_total=compatible_entropy_total,
            entropy_total_scale=float(np.sum(weight*abs(drift_entropy_rate))),work_scale=float(np.sum(weight*(abs(h*drift_entropy_rate)+abs(np.sum(j*pressure,axis=-1)))))))
    return result


def main():
    root=Path(__file__).resolve().parents[1];source=root/'outputs/current_transport_verification';out=root/'outputs/species_charge_probe';out.mkdir(exist_ok=False)
    identity=json.loads((source/'provenance.json').read_text());files=dict(identity['source_files']);name='examples/species_charge_probe.py';files[name]=hashlib.sha256((root/name).read_bytes()).hexdigest()
    for name,digest in files.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
    from .radial_work_records import freeze
    provenance=dict(source_files=files,source_sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),base_source_sha256=identity['source_sha256'],python=identity['python'],numpy=identity['numpy'],scipy=identity['scipy'])
    freeze(out,provenance);records=[]
    for i in range(4,8):
        path=source/f'case_{i:02d}.json';record=json.loads(path.read_text());raw=source/record['arrays_file']
        assert hashlib.sha256(raw.read_bytes()).hexdigest()==record['arrays_sha256'];arrays=dict(np.load(raw))
        low=measure(record['spec'],arrays,16);high=measure(record['spec'],arrays,20)
        row=dict(base_record=path.name,arrays_sha256=record['arrays_sha256'],spec=record['spec'],order16=low,order20=high);records.append(row)
        (out/f'case_{i:02d}.json').write_text(json.dumps(row,indent=2)+'\n');print(record['spec'],json.dumps(high),flush=True)
    (out/'summary.json').write_text(json.dumps(dict(source_sha256=provenance['source_sha256'],base_source_sha256=identity['source_sha256'],records=records),indent=2)+'\n')


if __name__=='__main__':main()
