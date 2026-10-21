"""Named experiments for examples/nucte_formation.py: `python -m examples.nucte_formation OUT preset=NAME [k=v ...]`.

Each preset is a dict of Nucte overrides (see examples/nucte_formation.py); explicit arguments override it.
Geometry and banks follow docs/NUCTE_BENCHMARK.md (NUCTE-III) and docs/DESIGN_SCAN.md (12-inch device).
"""
from dataclasses import replace
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

# Lifetime levers (docs/EXPERIMENT_PLAYGROUND.md, "Lifetime levers").
# Tuned flux core: the same solenoid on a 1 kV charge (150 J).
FLUX_CORE_1KV = replace(FLUX_CORE, bank=replace(FLUX_CORE.bank, voltage=1e3))
# Main-coil re-boost step for the 2 kV core: 8 x 6 uF at 3 kV (216 J) at 15 us, no crowbar of its own (start switches
# open at the current zero; the main crowbars carry on). A crowbarred re-boost mostly halves the crowbar resistance.
# Sized by replaying the coupled circuit with the recorded plasma image currents: the main current then stays within
# +2.9/-2.0 kA of the inner-tube control over 15-40 us (core alone: -6 kA by 40 us).
CORE_REBOOST = Bank(capacitance=48e-6, voltage=3e3, stray=3.75e-9, fire_time=15e-6, modules=8, crowbar=False)
# Slow mirror pair just beyond the main coil (z = +-0.635 m, r 0.17 m, 4 turns each in series) on 200 uF at 3 kV
# (0.9 kJ, ~90 us quarter period), fired with the bias lead, crowbarred at its current peak.
SLOW_MIRROR = AuxCoil(name="mirror", z_center=.635, length=.05, radius=.17, turns=4, pair=True,
                      bank=Bank(capacitance=200e-6, voltage=3e3, stray=.5e-6, fire_time=-100e-6, crowbar_decay_s=300e-6))

# Timed core program (docs/EXPERIMENT_PLAYGROUND.md, "Timed core program"). Same 40-turn core winding throughout.
# Double swing: one core bank charged negative and fired with the bias (-108 us) through a switch that conducts both
# ways, no crowbar until the positive peak (beyond 100 us). Its current rings to about -11 kA at formation (pre-bias,
# in the bias direction) and swings through zero near 60 us, so the positive flux ramp after formation (2.2 mWb by
# 100 us in vacuum) is about twice the peak core flux. 450 uF at -1.0 kV (225 J).
DOUBLE_SWING = replace(FLUX_CORE, bank=Bank(capacitance=450e-6, voltage=-1e3, stray=1e-6, fire_time=-108e-6,
                                            crowbar_decay_s=300e-6, unidirectional=False, crowbar_one_way=False,
                                            crowbar_time=140e-6))
# Staged core: four banks on the core winding, all at 1.05 kV, diode-isolated (each start switch opens at its current
# zero), fired at 4, 36.3, 59.8 and 77.2 us; crowbar on the last stage only. The capacitances hold the vacuum core loop
# voltage at 21-24 V from 10 to 100 us (2.1 mWb), where the single 1 kV bank falls from 22 to 14 V by 40 us and is spent
# at 66 us. 3.17 mF in all (1.75 kJ).
STAGED_CORE = replace(FLUX_CORE, bank=Bank(capacitance=723e-6, voltage=1.05e3, stray=1e-6, fire_time=4e-6, crowbar=False),
                      extra_banks=(Bank(capacitance=880e-6, voltage=1.05e3, stray=1e-6, fire_time=36.3e-6, crowbar=False),
                                   Bank(capacitance=1021e-6, voltage=1.05e3, stray=1e-6, fire_time=59.8e-6, crowbar=False),
                                   Bank(capacitance=550e-6, voltage=1.05e3, stray=1e-6, fire_time=77.2e-6,
                                        crowbar_decay_s=300e-6)))
# Shape program for the staged core: four main-coil boosts, 8 x 6 uF at 3.5 kV each (294 J), fired at 20 us and with
# core stages 2-4, each with its own one-way crowbar. In vacuum they hold the main current at 195-207 kA from 20 to
# 100 us (without them it sags to 133 kA).
MAIN_HOLD = tuple(Bank(capacitance=48e-6, voltage=3.5e3, stray=3.75e-9, fire_time=t, modules=8, crowbar_inductance=1e-9,
                       crowbar_decay_s=300e-6) for t in (20e-6, 36.3e-6, 59.8e-6, 77.2e-6))

# Toroidal field (docs/EXPERIMENT_PLAYGROUND.md, "Toroidal field"): an axial current on a 1 cm rod inside the inner
# tube, returning through a 20 cm coaxial cage outside the main coil, so B_theta = mu0 I/(2 pi r) between them and
# nothing outside. Line L = mu0 l/(2 pi) ln(b/a) = 1.32 uH over 2.2 m; bank 20 uF with 0.1 uH stray (quarter period
# 8.4 us), one-way crowbar at the current peak, 200 us L/R (6.6 mOhm in the crowbar loop). The charge sets the peak:
# 2.66 / 5.33 / 10.65 kV give 10 / 20 / 40 kA (71 J / 0.28 kJ / 1.13 kJ), B_theta(8 cm) = 0.025 / 0.05 / 0.1 T.
# "before": fired at -20 us, peaking at -11.6 us and crowbarred ahead of pre-ionization (fired -50 to -10 us), so the
# gas is ionized with the field in place (0.94 of the peak at t = 0). "after": fired at 6 us, once the FRC has
# formed, peaking at 14.4 us.
def rod_bank(peak_ka, fire_time):
    return Bank(capacitance=20e-6, voltage={10: 2.663e3, 20: 5.326e3, 40: 10.65e3}[peak_ka], stray=.1e-6,
                fire_time=fire_time, crowbar_decay_s=200e-6)


ROD_FIRE = {"before": -20e-6, "after": 6e-6}

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
    # Lifetime levers. Flux core at 1 kV: about the loop voltage that holds the trapped flux.
    "frc12_flux_core_1kv": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "aux_coils": (FLUX_CORE_1KV,)},
    # 2 kV flux core plus a main-coil re-boost that holds the external field near the inner-tube control.
    "frc12_core_reboost": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "aux_coils": (FLUX_CORE,),
                           "extra_main_banks": (CORE_REBOOST,)},
    # Option A with the slow mirror pair energized before the main bank. The parallel main straps spread the mirror's
    # flux (+0.043 T at the centre at t = 0), so the bias is raised to -0.0946 T on the main coil alone (-0.05 T in
    # total; bias bank -4.4 kV, 19 kJ) and peak_t lowered to 0.169 (0.2 T total peak).
    "frc12_mirror_bias": {**FRC12, **UNIFORM, **H2_FILL, "bias_t": -.0946, "peak_t": .169, "aux_coils": (SLOW_MIRROR,)},
    # Timed core program: double-swing core (pre-biased negative with the bias, swung through zero after formation).
    "frc12_core_double_swing": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "aux_coils": (DOUBLE_SWING,)},
    # Timed core program: four staged core banks holding the core loop voltage near the measured ~21-25 V loss.
    "frc12_core_staged": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "aux_coils": (STAGED_CORE,)},
    # Staged core plus main-coil boosts that hold the external field (shape program).
    "frc12_core_staged_shape": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "aux_coils": (STAGED_CORE,),
                                "extra_main_banks": MAIN_HOLD},
    # Toroidal field from the rod (2.5D, implies toroidal): frc12_inner_tube plus the rod bank, 10/20/40 kA peak,
    # fired before (-20 us) or after (6 us) formation. Without Hall, so outputs/fine_inner_tube is the exact control.
    **{f"frc12_rod{ka}_{when}": {**FRC12, **UNIFORM, **H2_FILL, **INNER_TUBE, "rod_bank": rod_bank(ka, fire)}
       for ka in (10, 20, 40) for when, fire in ROD_FIRE.items()},
}
