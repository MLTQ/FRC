"""Named experiments for examples/nucte_formation.py: `python -m examples.nucte_formation OUT preset=NAME [k=v ...]`.

Each preset is a dict of Nucte overrides (see examples/nucte_formation.py); explicit arguments override it.
Geometry and banks follow docs/NUCTE_BENCHMARK.md (NUCTE-III) and docs/DESIGN_SCAN.md (12-inch device).
"""
from examples.nucte_formation import AuxCoil, Bank

# NUCTE-III (Takahashi 2000). Case 1: driver defaults (two-radius coil, mirror 1.10).
NUCTE_CASE1 = dict()
# Case 6: graded coil matching the published mirror ratio and field index; bias and peak from NUCTE's own
# flux estimate and pressure balance; bank coupled to the plasma.
NUCTE_CASE6 = dict(coil_radius=.18, end_radius=.1267, taper_start=.1, filaments=8, bias_t=-.05, peak_t=.3,
                   conduction_limit=.1, drive="coupled", end_s=30e-6)

# The user's 12-inch device: 12" ID tube (r_w 15.24 cm) x 2.2 m, single-turn coil 1.2 m (24 x 5 cm straps on one
# collector), eight 6 uF modules (30 nH each, 8 nH crowbars, ~2 mOhm, unidirectional), NUCTE-like bias bank.
FRC12 = dict(wall_radius=.1524, half_length=1.1, coil_radius=.165, end_radius=.165, centre_straps=24, end_straps=0,
             strap_width=.05, filaments=4, nr=24, nz=144, exterior_cells=2, end_cells=2, pi_half_length=.6,
             bias_t=-.05, peak_t=.2, drive="coupled", main_capacitance=48e-6, main_stray=3.75e-9, main_modules=8,
             crowbar_inductance=1e-9, crowbar_decay_s=300e-6, quarter_period_s=3.1e-6, conduction_limit=.1, end_s=30e-6)
H2_FILL = dict(gas="H2", fill_pa=.4533)        # 3.4 mTorr
D2_FILL = dict(gas="D2", fill_pa=.4666)        # 3.5 mTorr
UNIFORM = dict()                               # coil hugs the tube, mirror ratio 1.00
GRADED = dict(coil_radius=.20, taper_start=.1)  # 20 cm at |z| <= 0.1 m tapering to 16.5 cm at the ends, mirror 1.27

INNER_TUBE = dict(inner_radius=.0381)           # 3" OD inner tube (snapped to the grid: 3.81 cm at nr=24)
# 40-turn solenoid, r 2.5 cm, 0.8 m long, inside the inner tube; 300 uF at 2 kV (0.6 kJ), quarter period ~40 us.
FLUX_CORE = AuxCoil(name="core", z_center=0., length=.8, radius=.025, turns=40,
                    bank=Bank(capacitance=300e-6, voltage=2e3, stray=1e-6, fire_time=4e-6, crowbar_decay_s=300e-6))
CORE_PUNCH = AuxCoil(name="punch", z_center=0., length=.3, radius=.025, turns=10,
                     bank=Bank(capacitance=1.6e-6, voltage=19e3, stray=.2e-6, fire_time=8e-6, crowbar_decay_s=100e-6))

PRESETS = {
    "nucte_case1": NUCTE_CASE1,
    "nucte_case6": NUCTE_CASE6,
    # Option A: uniform coil, H2, 0.2 T (11.6 kV, 3.2 kJ). Raise peak_t=0.3 for 16.2 kV / 6.3 kJ.
    "frc12_uniform_h2": {**FRC12, **UNIFORM, **H2_FILL},
    # Option B: mild graded coil, D2, 0.2 T (15.2 kV, 5.6 kJ); peak_t=0.3 gives 21.3 kV / 10.9 kJ.
    "frc12_graded_d2": {**FRC12, **GRADED, **D2_FILL},
    "frc12_graded_h2": {**FRC12, **GRADED, **H2_FILL},
    "frc12_uniform_d2": {**FRC12, **UNIFORM, **D2_FILL},
    # Playground: a second bank on the main coil fired after formation to restore the decaying field.
    "frc12_reboost": {**FRC12, **UNIFORM, **H2_FILL,
                      "extra_main_banks": (Bank(capacitance=48e-6, voltage=6e3, stray=3.75e-9, fire_time=15e-6, modules=8,
                                                crowbar_inductance=1e-9, crowbar_decay_s=300e-6),)},
    # Playground: end "plug" coils (2 turns each side, in series) just beyond the main coil, fired once the FRC
    # has formed, to raise the end mirror and hold the plasma axially.
    "frc12_end_plugs": {**FRC12, **UNIFORM, **H2_FILL,
                        "aux_coils": (AuxCoil(name="plugs", z_center=.635, length=.05, radius=.17, turns=1, pair=True,
                                              bank=Bank(capacitance=20e-6, voltage=15e3, stray=60e-9, fire_time=3e-6,
                                                        crowbar_decay_s=200e-6)),)},
    # Playground: a coaxial 3-inch OD glass tube on the axis (gas only between it and the 12-inch tube). Control
    # for the flux-core runs below.
    "frc12_inner_tube": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE},
    # Playground: a pulsed solenoid inside the inner tube ("flux core"), fired once the FRC has formed, in the
    # main-field direction. Its rising flux lifts the separatrix level at the tube while the O-point flux is
    # frozen into the plasma, so it adds closed flux against resistive decay (loop voltage ~ dPhi/dt).
    "frc12_flux_core": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "aux_coils": (FLUX_CORE,)},
    # Playground: a fast "punch" from the inner tube: a short 10-turn coil on a 1.6 uF / 19 kV bank (2 us quarter
    # period), fired at 8 us. Same polarity as the main field; voltage < 0 reverses it (removes closed flux).
    "frc12_core_punch": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "aux_coils": (CORE_PUNCH,)},
}
