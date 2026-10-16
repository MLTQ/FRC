# Argon radiation-trapping experimental reference

Primary report: J. E. Talmage Jr., G. S. Hurst, M. G. Payne and E. B. Wagner,
*Transport of Resonance Radiation in Argon*, ORNL-TM-4368 (October 1973).
Public source: https://www.osti.gov/servlets/purl/4429456 . The downloaded original
is retained as `ornl-tm-4368.pdf`; SHA256 is in `source.json`. No account or terms
acceptance was used. This reference is unrelated to the pending LXCat download.

`table1.csv` is transcribed from printed page 30 (PDF page 41), visually checked
against its rendered page. Pressures are Torr; rates are inverse microseconds.
The run number at 0.093 Torr is 194 (OCR was incorrect). A true impurity flag
preserves the table's superscript a; it is not a numerical correction. Printed
page 54 estimates errors below ±3%, except ±20% at 0.0015 Torr; these are stated
measurement error estimates, not specified standard deviations. Later a-marked
runs had changed pumping/impurity conditions, with rates roughly 0.03/us higher
below 0.3 Torr. Do not silently merge them with clean runs or subtract that value.

Relevant geometry (printed pages 12-13 and Figure 9, page 47): cylinder length
7.0 cm, radius 1.11 cm; axial proton-beam excitation with 3 mm collimator upstream;
side aperture diameter 0.714 cm halfway along the cell, LiF flush with the wall;
monochromator entrance slit width 50 micrometres and height 1 cm. Figure 9 shows
the narrow axial slit and vertical aperture limiting the observed volume.
The actual beam profile, angular acceptance and pulse response need treatment,
not an unqualified point source or total-wall flux. Pressure broadening,
quenching, spectral response and thermal conditions must remain explicit.

The observable is the EARLY-time semilog slope of detector counts, selected by
extending the fit window until the slope changes (pages 27, 53-54). Above about
0.5 Torr, a distinct slower late-time slope appears. It is neither the whole-
gas survival curve nor inverse mean photon escape time. Before comparing model
predictions, reproduce the viewed volume and fit-window dependence. The current
1 m centre-source component study has not been validated against this table.
