from decimal import Decimal as D, localcontext
import json
import argparse
from .implicit_plasma_records import provenance, freeze
from .precise_material import save_states
from pathlib import Path
import numpy as np
from scipy.constants import atomic_mass
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_vector_moments import VectorizedAtomicMoments
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination
from plasma_bench.precise_geometry import coordinates
from plasma_bench.material_reactions import ReactionOverlapState
import plasma_bench.precise_reactions as pr

parser=argparse.ArgumentParser()
parser.add_argument('destination')
args=parser.parse_args()
out=Path(args.destination)
if out.exists(): raise ValueError('Fresh diagnostic output required')
out.mkdir(parents=True)
identity=provenance(['examples.reaction_projection_diagnosis'])
freeze(out,identity)
report=dict(scope='Exact saved source and knots: independent 100-digit projection, local errors, and explicit candidate comparisons. No formation claim.', arithmetic=[])
source=VectorizedAtomicMoments(AtomicNetwork(ImpactIonization.from_file('data/atomic/argon/argon.json'),RadiativeRecombination.from_file('data/atomic/recombination/argon_rr.json')),39.948*atomic_mass,[0.,0.,0.])
overlap=pr.PreciseReactionOverlap(source,1e14)
p=Path('outputs/reactive_material_argon_retained_knots_10ns')
with np.load(p/'rejected_source.npz') as z:
 state=ReactionOverlapState(coordinates(z['edges_decimal_m']),coordinates(z['points_decimal_m']),z['weights_m'],z['moments'],0.)
 knots={int(k.split('species')[1].split('_')[0]):coordinates(z[k]) for k in z.files if k.startswith('velocity_knots')}
loads=[]
original=pr.solve_mass_banded
def capture(matrix,load):
 loads.append((matrix.copy(),load.copy()))
 return original(matrix,load)
pr.solve_mass_banded=capture
model,result,record=pr.project_precise_reactions(overlap,state,velocity_cells=32,velocity_knots=knots,energy_projection='projected_frame')

def dec(x): return D(float(x))
def solve(diag,off,rhs):
 d=diag.copy(); b=rhs.copy()
 for i in range(1,len(d)):
  ratio=off[i-1]/d[i-1];d[i]-=ratio*off[i-1];b[i]-=ratio*b[i-1]
 x=[D(0)]*len(d);x[-1]=b[-1]/d[-1]
 for i in range(len(d)-2,-1,-1):x[i]=(b[i]-off[i]*x[i+1])/d[i]
 return x
with localcontext() as ctx:
 ctx.prec=100
 for index in (0,9):
  bidx=list(model.species_indices).index(index);m=model.models[bidx];s=result.blocks[bidx]
  matrix,load=loads[bidx];nv=m.velocity_cells+1
  exact_rounded=solve(list(map(dec,matrix[1])),list(map(dec,matrix[0,1:])),list(map(dec,load[:,2])))
  diag=[D(0)]*nv;off=[D(0)]*(nv-1);rhs=[D(0)]*nv
  scale=dec(state.moments[...,overlap.number[index]].max());mass=dec(source.mass[index]);pm=dec(source.momentum_scale)
  start=np.searchsorted(state.edges,m.reference_faces[0]);coarse=m.reference_faces[m.velocity_indices]
  for fine in range(m.cells):
   c=int(m.owner[fine]);k=start+fine;lo,hi=m.reference_faces[fine:fine+2];h=coarse[c+1]-coarse[c]
   l0=(coarse[c+1]-lo)/h;l1=(coarse[c+1]-hi)/h;r0=(lo-coarse[c])/h;r1=(hi-coarse[c])/h
   nn=sum(dec(w)*dec(n)/scale for w,n in zip(state.weights[k],state.moments[k,:,overlap.number[index]]))
   diag[c]+=mass*nn*(l0*l0+l0*l1+l1*l1)/3
   diag[c+1]+=mass*nn*(r0*r0+r0*r1+r1*r1)/3
   off[c]+=mass*nn*(2*l0*r0+l0*r1+l1*r0+2*l1*r1)/6
   for q,x in enumerate(state.points[k]):
    p=dec(state.weights[k,q])*dec(state.moments[k,q,source.momentum[index,2]])/scale*pm
    rhs[c]+=p*(coarse[c+1]-x)/h;rhs[c+1]+=p*(x-coarse[c])/h
  exact=solve(diag,off,rhs);v=s.velocity[0,:,2];high=np.array([float(x) for x in exact]);rounded=np.array([float(x) for x in exact_rounded])
  # Independent constant-velocity limit using the assembled mass row sums.
  const=original(matrix,300.*(matrix[1]+np.r_[matrix[0,1:],0]+np.r_[0,matrix[2,:-1]]))
  report['arithmetic'].append(dict(species=index,production_range=[v.min(),v.max()],rounded_100digit_range=[rounded.min(),rounded.max()],full_100digit_range=[high.min(),high.max()],relative_rounded_solve_error=float(abs(v-rounded).max()/max(abs(rounded).max(),1)),absolute_assembly_velocity_error=float(abs(rounded-high).max()),constant_velocity_error=float(abs(const-300).max()),worst_nodes=np.argsort(abs(v-high))[-5:].tolist(),production=v.tolist(),independent=high.tolist()))

pr.solve_mass_banded=original
from scipy.constants import elementary_charge as EV
from plasma_bench.precise_geometry import fraction
import matplotlib.pyplot as plt
report['methods']={}
for method in ('galerkin','bounded_moments','constrained_velocity'):
 controls=dict(energy_projection='projected_frame') if method=='galerkin' else dict(projection_method=method)
 md,st,rec=pr.project_precise_reactions(overlap,state,velocity_cells=32,velocity_knots=knots,**controls)
 rows=[]
 for b,m,s in zip(md.species_indices,md.models,st.blocks):
  first=np.searchsorted(state.edges,m.reference_faces[0]);last=first+m.cells
  number=state.moments[first:last,:,overlap.number[b]]
  velocity=np.divide(state.moments[first:last,:,source.momentum[b]],number[...,None],out=np.zeros(number.shape+(3,)),where=number[...,None]>0)*(source.momentum_scale/source.mass[b])
  coarse=m.reference_faces[m.velocity_indices]
  a=fraction(state.points[first:last],coarse[m.owner,None],coarse[m.owner+1,None])
  actual=(1-a[...,None])*s.velocity[0,m.owner,None]+a[...,None]*s.velocity[0,m.owner+1,None]
  temperature=(m.gamma-1)*s.thermal[0]/m.relative_numbers/EV
  source_scale=state.moments[...,overlap.number[b]].max()
  expected_number=np.sum(state.weights[first:last]*(number/source_scale))
  p_source=np.sum(state.weights[first:last,:,None]*(state.moments[first:last,:,source.momentum[b]]/source_scale*source.momentum_scale),axis=(0,1))
  vl,vr=m.fine_velocity(s.velocity)
  p_target=np.sum((m.mass[0]*m.relative_numbers)[:,None]*(vl+vr)/2,axis=0)
  p_norm=np.sum(state.weights[first:last,:,None]*abs(state.moments[first:last,:,source.momentum[b]]/source_scale*source.momentum_scale))
  rows.append(dict(species=int(b),number_scale_m2=m.number_scale,velocity_cells=m.velocity_cells,density_cells=m.cells,
    source_velocity_range_m_s=[float(velocity[number>0].min()),float(velocity[number>0].max())],
    projected_velocity_range_m_s=[float(s.velocity.min()),float(s.velocity.max())],
    maximum_velocity_error_m_s=float(abs(actual-velocity)[number>0].max()),
    minimum_temperature_ev=float(temperature.min()),maximum_temperature_ev=float(temperature.max()),
    normalized_number_error=float(abs(m.relative_numbers.sum()-expected_number)/expected_number),
    normalized_momentum_error=float(abs(p_target-p_source).max()/p_norm) if p_norm else float(abs(p_target-p_source).max())))
 source_total=overlap.totals(state)
 target=md.totals(st)
 energy=target['energy_j_m2']+rec['binding_energy_j_m2']+rec['photon_energy_j_m2']-rec['input_heat_j_m2']
 report['methods'][method]=dict(blocks=rows,projection=rec,relative_total_energy_error=(energy-source_total['energy_j_m2'])/source_total['energy_j_m2'])
 save_states(out/(method+'.npz'),[md],[st])
 print(method, 'relative energy error',report['methods'][method]['relative_total_energy_error'],flush=True)
report['source_sha256']=identity['source_sha256']
report['terminal']=True
(out/'summary.json').write_text(json.dumps(report,indent=2,default=lambda x:x.item())+'\n')
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for method,row in report['methods'].items():
 ids=[r['species'] for r in row['blocks']]
 axes[0].semilogy(ids,[r['maximum_velocity_error_m_s'] for r in row['blocks']],'.-',label=method)
 axes[1].semilogy(ids,[r['maximum_temperature_ev'] for r in row['blocks']],'.-',label=method)
axes[0].set_ylabel('Maximum source-point velocity error [m/s]')
axes[1].set_ylabel('Maximum projected temperature [eV]')
for ax in axes:
 ax.set_xlabel('Species index (19 = electron)');ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.suptitle('Same rejected source and inherited knots — snapshot comparisons')
fig.savefig(out/'projection_comparison.png',dpi=160)
fig.savefig(out/'projection_comparison.pdf')
plt.close(fig)
if identity!=provenance(['examples.reaction_projection_diagnosis']):raise RuntimeError('Source changed during snapshot study')
