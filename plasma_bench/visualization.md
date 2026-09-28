# Scientific visualization
`snapshot` produces two r-z field maps with psi contours and O/X markers, probe
traces and energy histories. Fixed color limits across a run permit honest time
comparison; only the chamber is displayed, excluding singular source filaments.
`viewer` exports an HTML animation with play/pause/step/slider controls and
assumption warnings; optional GIF uses the same data. It subsamples snapshots
for display only. `comparison_plot` overlays three normalized experiments and
vacuum counterfactual center probes. `phase_diagram` shows four explicit evidence
classes without claiming FRC success. No field or topology is artistically
generated. Density, pressure and velocity plots await actual evolving models.
The HTML export replaces Matplotlib's CDN icon dependency with text labels so
all frames and playback controls work offline and remain readable.
Contours use a contrasting outline across the diverging colormap; extra levels
around detected O points make small resolved closed regions visible. Color
limits remain fixed, while contour levels adapt to each actual flux slice.
