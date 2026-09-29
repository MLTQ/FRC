# Boundary verification report assembly

Requires the completed 12-case pulse study with hashed field archives, 15
analytic skin cases, four gas-end/matching controls and the nr=384 strip shock.
Earlier compatible shocks and strip controls are read with their original
records. No running process is started or restarted by report assembly.
Engine and driver hashes are checked, and incomplete cases are rejected.
End records must match the full source/probe/pulse/end provenance, not just the
end driver's own hash; their field archives are checked on load as well.
All reused numerical records must also match the frozen engine environment
(Python, platform and recorded package versions), including earlier shock cases.

Builds full-field volume-weighted comparisons, same-grid method/time/magnetic
matching comparisons, and coincident-window gas-end sensitivity. Writes a
lossless compact metrics record, the detailed Markdown report and scientific
figure. Analysis/report/plot source hashes are separate from solver hashes.

Interpretation distinguishes analytic error, inter-grid differences, numerical
dissipation, changed physical gas boundaries and full-space response budgets.
Kinematic skin tests do not stand in for gas dynamics; shock energy improvements
do not establish formation. The report retains unresolved physical and numerical
requirements and does not promote a production default or reopen Stage 2.
The interpretation explicitly computes the finest same-grid method difference
and both finest adjacent-grid velocity differences; favorable energy accounting
cannot hide a substantial change in the driven gas response.
Fine-grid density ranges, measured current widths in metres/cells and local
numerical/physical eta ratios accompany the integrated work and heat tables.
