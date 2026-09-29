"""Explicit conditional argon regimes and a vanishing-density edge diagnostic."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import atomic_mass,proton_mass,mu_0 as MU0
from plasma_bench.coulomb_transport import regime
from .verification_records import write_record


def study():
    sources=[Path('plasma_bench/coulomb_transport.py'),Path(__file__),Path('examples/verification_records.py')]
    source_hashes={str(p.relative_to(Path.cwd()) if p.is_absolute() else p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    cases=[];mi=39.948*atomic_mass
    for ne in np.logspace(16,24,5):
        for te in (1.,5.,20.,100.):
            for length in (.035,.001):
                r=regime(ne,te,1.,mi,field=.005,current_density=.005/(MU0*length),length=length,timescale=5e-6)
                cases.append({'name':f'Ar_ne{ne:g}_Te{te:g}_L{length:g}',
                              'interpretation':'Assumed Ar+ with prescribed state; J=B/(mu0 L) is an order-of-magnitude gradient diagnostic, not evolved current',**r})
    for name,mass,x in [('default_argon',mi,1.),('hydrogen_comparison',proton_mass,1.),('partial_argon',mi,.1)]:
        cases.append({'name':name,**regime(1e20,5.,1.,mass,field=.005,length=.035,timescale=5e-6,ionization=x)})
    edge=[]
    for xi in 1-np.logspace(-1,-5,41):
        ne=1e20*(1-xi*xi);length=.035*(1-xi*xi)/(2*xi)
        edge.append({'initial_radius_fraction':float(xi),'distance_to_edge_m':float(.035*(1-xi)),
                     **regime(ne,5.,1.,mi,field=.005,current_density=.005/(MU0*.035),length=length,timescale=5e-6)})
    return {'completed':True,'sources':source_hashes,'units':'SI except explicitly eV temperatures',
            'interpretation':'Conditional scans, not chosen lab parameters or evolved formation. Fixed charge/ionization is an input, not an equilibrium prediction.',
            'cases':cases,'edge':{'density':'ne=1e20*(1-(s/R)^2), R=.035 m; exact exterior vacuum excluded from plasma samples',
              'temperature':'Te=5 eV, Ti=1 eV fixed; not a thermal solution',
              'length':'ne/abs(dne/ds)','current':'fixed gradient scale .005/(mu0*.035) A/m², not a self-consistent edge current',
              'purpose':'Expose loss of local closure; no claim that the prescribed edge is physically attainable','samples':edge}}


def run(directory='outputs/current_verification'):
    result=study();directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    write_record(directory/'transport.json',result)
    print('Recorded',len(result['cases']),'conditional states and',len(result['edge']['samples']),'edge samples',flush=True)
    return result


if __name__=='__main__':run()
