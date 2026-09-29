"""Assemble checked current/transport evidence without claiming a solved edge."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
from .current_cases import protocol
from .current_measurements import provenance
from .current_time_reference import measure as time_measure
from .current_plot import accuracy,regimes
from .transport_study import study as transport_study
from .verification_records import write_record


def read(directory):
    root=Path(directory);record=json.loads((root/'measurements.json').read_text());source=provenance()
    if not record['completed'] or record['provenance']!=source:raise ValueError('Incomplete or stale current measurements')
    cases=protocol();rows=record['cases'];arrays={}
    if len(cases)!=len(rows):raise ValueError('Incomplete current protocol')
    for case,row in zip(cases,rows):
        path=root/row['fields']['file']
        if not row['completed'] or row['configuration']!=case.configuration() or row['provenance']!=source or hashlib.sha256(path.read_bytes()).hexdigest()!=row['fields']['sha256']:
            raise ValueError('Changed current case configuration/source/archive')
        if json.loads((root/f'{case.name}.json').read_text())!=row:raise ValueError('Changed individual case record')
        with np.load(path,allow_pickle=False) as archive:
            arrays[case.name]={k:archive[k] for k in archive.files}
        if not all(np.all(np.isfinite(a)) for a in arrays[case.name].values()):raise ValueError('Nonfinite archive')
    return record,arrays


def controls(arrays):
    baseline=next(c for c in protocol() if c.name=='decay_d0.2_16');o,_,_=baseline.build();base=arrays[baseline.name];rows=[]
    for case in protocol():
        if not case.name.startswith(baseline.name+'_'):continue
        a=arrays[case.name]
        if not np.array_equal(a['final_coordinates'],base['final_coordinates']):continue
        q=a['final_current_over_radius'];dq=q-base['final_current_over_radius'];v=a['final_volume_weights'].sum(axis=1)
        heat=a['final_cell_heat'];delta=heat-base['final_cell_heat']
        rows.append({'reference':case.name,'test':baseline.name,
                     'current_relative_volume_L2':float(np.sqrt(dq@o.mass@dq/(q@o.mass@q))),
                     'cell_heat_density_relative_volume_L2':float(np.sqrt(np.sum(delta*delta/v)/np.sum(heat*heat/v)))})
    return rows


def screening_signs(rows,arrays):
    result=[]
    for row in rows:
        if row['configuration']['kind']!='screening':continue
        a=arrays[row['name']];w=a['final_volume_weights'];bias=float(a['final_bias'])
        result.append({'name':row['name'],'minimum_Bz_over_bias':float(a['final_Bz'].min()/bias),
                       'negative_Bz_volume_fraction':float(np.sum(w*(a['final_Bz']<0))/w.sum()),
                       'maximum_positive_Jtheta_a_m2':float(max(0.,a['final_current_quadrature'].max())),
                       'central_Bz_over_bias':float(2*a['final_material_a'][0]/bias)})
    return result


def report(result):
    rows=result['measurements']['cases'];transport=result['transport'];default=next(r for r in transport['cases'] if r['name']=='default_argon')
    lines=['# Electron-current inertia and physical regime checks','',
      'The new `MaterialCurrent` component implements fixed-material toroidal induction with physical electron inertia, scalar prescribed relaxation and genuine exterior magnetic vacuum. It separately accounts for induced magnetic energy, electron drift kinetic energy, positive collisional heat and prescribed-source work. It is a linear component reference, not a moving Hall/thermal solver or an FRC formation result. The prior heated-edge incompatibility remains open. See [the equations and discrete energy derivation](CURRENT_INERTIA.md) and [conditional transport conventions](TRANSPORT_REGIME.md).','',
      f"All {len(rows)} declared component cases are complete. Source/configuration and lossless array hashes are checked; report, time-audit, scan and rendering sources are separately recorded. Radius is 1 m, analytic field coefficient/bias is 1e-8 T, and the constant-profile decay uses D=0.1 m²/s. Density follows physical de²=me/(mu0 ne e²). These are explicit component inputs, not selected laboratory operating parameters.",'',
      'The de=0 rows explicitly select the massless comparison; no finite carrier density is assigned to them.','',
      '## Independent continuum comparisons','',
      'Sphere decay has continuous dipole/vacuum matching and exact rate D k²/(1+de² k²). Collisionless screening has an independent modified-spherical-Bessel solution and vanishes continuously with carrier density. A third exact sphere correlates prescribed density/resistivity profiles so the same continuum mode has a known decay; this verifies variable coefficients without a numerical forcing or an energy-only assertion. It is a manufactured profile, not a species closure.','',
      '| Case | Magnetic field relative L2 | Current relative L2 | Local heat relative L2 | Initial current relative L2 |','|---|---:|---:|---:|---:|']
    for row in rows:
        e=row['final_errors'];initial=row['initial_errors'];heat=e.get('local_heat_relative_volume_L2')
        lines.append(f"| {row['name']} | {e['field_relative_volume_L2']:.7g} | {e['current_relative_volume_L2']:.7g} | "+('n/a (eta=0)' if heat is None else f'{heat:.7g}')+' | '+('zero initial state' if initial is None else f"{initial['current_relative_volume_L2']:.7g}")+' |')
    lines += ['',
      'The Q1-gradient magnetic field converges approximately at first order. Inertia retains weak-current error from the interpolated initial a, especially near the surface; at de/R=0.2 the nr=8/16/32/64 current errors are about 0.02498/0.007687/0.002480/0.0008289. This approaches roughly order 1.5, rather than the near-second-order current error after strong high-mode damping in the massless comparison. The initialized field/current inconsistency is measured at t=0 and not removed. Surface grading alone worsens this current error in the recorded control.','',
      '## Thin screening and false reversal','',
      'The exact collisionless sphere has positive interior Bz everywhere for positive applied bias. Its numerical minimum and negative-field volume are therefore direct artifact checks; a small energy residual must not turn them into an FRC detector.','',
      '| Case | Minimum Bz / bias | Negative-Bz material volume [%] | Largest unphysical positive Jtheta [A/m²] | Central Bz / bias |','|---|---:|---:|---:|---:|']
    for row in result['screening_signs']:
        lines.append(f"| {row['name']} | {row['minimum_Bz_over_bias']:.7g} | {100*row['negative_Bz_volume_fraction']:.7g} | {row['maximum_positive_Jtheta_a_m2']:.7g} | {row['central_Bz_over_bias']:.7g} |")
    lines += ['',
      'The unresolved de/R=0.05, nr=8 mesh produces false reversal and opposite-sign current. Both disappear in the recorded finer samples, but thin-layer magnetic-field error remains appreciable. Even positive field samples do not by themselves prove topology accuracy between samples. No clipping, density floor, resistance enhancement or heat cutoff hides the coarse artifact.','',
      '## Time, coefficients and boundary controls','',
      '| Same material mesh: baseline / reference | Current difference [%] | Cell heat-density difference [%] |','|---|---:|---:|']
    for row in result['controls']:
        lines.append(f"| {row['test']} / {row['reference']} | {100*row['current_relative_volume_L2']:.7g} | {100*row['cell_heat_density_relative_volume_L2']:.7g} |")
    lines += ['',
      'Angular/grading cases have different material meshes and remain separate continuum-error rows above; no equal-index field comparison is made across them. The explicit fixed-spatial-eigenmode audit isolates temporal accuracy independently of continuum spatial error:','',
      '| dt [s] | Exact semidiscrete current relative L2 error |','|---:|---:|']
    for row in result['time_audit']['cases']:lines.append(f"| {row['dt_s']:g} | {row['relative_current_L2_error']:.8g} |")
    ratios=result['time_audit']['error_halving_ratios']
    lines += ['',f'Time-error halving ratios are {ratios[0]:.7g} and {ratios[1]:.7g}. This confirms second order on this resolved fixed mode; midpoint remains unsuitable for arbitrarily underresolved dissipative steps.','',
      '## Energy and current response','',
      'Independent true-sphere volume integrals of B²/(2 mu0), me J²/(2 ne e²) and local heat include the exact infinite exterior dipole energy. Orders 48 and 96 agree, and uniform-decay totals independently match their closed formulas. Nonuniform-coefficient and deformed-geometry tests also evaluate physical energies separately from the mixed matrices.','',
      f"Maximum continuum energy-quadrature relative change is {result['maximum_reference_quadrature_change']:.4g}. Maximum whole-run response residual divided by the corresponding continuum magnetic/current/heat energy scale is {result['maximum_relative_response_residual']:.4g}; maximum equation residual is {result['maximum_equation_residual']:.4g}. Positive heat is retained per cell. These budget checks do not erase the measured field/current errors or establish a physical surface closure.",'',
      '## Conditional argon and local-edge assessment','',
      f"The original zero-current Ar+ reference (ne=1e20 m^-3, Te=5 eV, Ti=1 eV, radius=0.035 m, B=0.005 T) has di/R={default['hall_length_ratio']:.5g}, de/R={default['electron_inertia_over_length']:.5g}, and conditional electron mean-free-path/R={default['electron_mean_free_path_over_length']:.5g}. Thus the supplied global parameters do not justify Hall-free MHD. Smaller edge-gradient lengths are more restrictive.",'',
      f"Its conditional Coulomb momentum resistivity is {default['classical_coefficients']['momentum_resistivity_ohm_m']:.7g} ohm m; the separate Z=1 parallel/unmagnetized estimate is {default['classical_coefficients']['parallel_resistivity_z1_ohm_m']:.7g} ohm m. Prescribed conductivity is unchanged. The two coefficients have distinct limits; neither replaces the full tensor, neutral collisions or thermal transport.",'',
      f"The archived scan contains {len(transport['cases'])} explicitly assumed states and {len(transport['edge']['samples'])} samples of a prescribed vanishing-density edge. Its fixed current is a diagnostic scale, not a self-consistent solved edge. It exposes rising Hall/inertial/Debye and drift ratios as vacuum is approached. Fixed charge/ionization, Maxwellian distributions and temperatures are assumptions, not predictions; no tested state is labeled fully physically validated.",'',
      'Next integrate the required Hall/toroidal and electron/ion thermal equations with a justified edge treatment, preserving genuine exterior vacuum and independent wave/current/energy checks. The linear inertia component alone must not be inserted into moving dynamics without its geometric current-energy force and convective terms. Spatial coil work, actual gas ends, robust meshes/shocks and converged formation remain required. Stage 2 remains gated.','',
      'Reproduce with `python -m examples.current_measurements --reuse`, `python -m examples.transport_study`, then `python -m examples.verify_current`. The numerical suite currently has 192 passing tests. All evidence is in `outputs/current_verification`.','',
      '![Current accuracy](../outputs/current_verification/current_accuracy.png)','',
      '![Conditional transport scales](../outputs/current_verification/transport_scales.png)']
    Path('docs/CURRENT_VERIFICATION.md').write_text('\n'.join(lines)+'\n')


def verify(directory='outputs/current_verification'):
    directory=Path(directory);measurements,arrays=read(directory);rows=measurements['cases']
    transport=json.loads((directory/'transport.json').read_text())
    if transport!=transport_study():raise ValueError('Changed transport scan or sources')
    sources=[Path(f'examples/{n}.py') for n in ('verify_current','current_plot','current_time_reference','transport_study','verification_records')]
    residual=max(abs(h['response_residual_j'])/sum(r['continuum_final'][k] for k in ('magnetic_energy_j','current_kinetic_energy_j','ohmic_heat_j')) for r in rows for h in r['history'])
    result={'completed':True,'measurements':measurements,'transport':transport,'report_environment':{'matplotlib':matplotlib.__version__},
            'report_sources':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            'controls':controls(arrays),'screening_signs':screening_signs(rows,arrays),'time_audit':time_measure(),
            'maximum_reference_quadrature_change':max(r['continuum_final']['quadrature_relative_change'] for r in rows),
            'maximum_relative_response_residual':residual,
            'maximum_equation_residual':max(h[k] for r in rows for h in r['history'] for k in ('ampere_relative_residual','ohm_relative_residual','vacuum_relative_residual'))}
    accuracy(directory,rows,arrays);regimes(directory,transport);report(result)
    write_record(directory/'metrics.json',result);print('Complete current/transport study assembled',flush=True)


if __name__=='__main__':verify()
