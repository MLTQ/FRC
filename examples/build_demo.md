# Reproduce the delivered demonstration
`main` writes the default configuration, verified-stage vacuum example, three
normalized conductor runs, static figures, offline playback/GIF, 27-case sweep
and findings. `findings` renders measured summaries into the research report,
with assumptions and convergence limitations kept adjacent to conclusions.
Run only after the test suite; output locations are explicit under outputs/.
Markdown tables use unambiguous absolute-value labels for current magnitude.
The findings explicitly include the measured topology-duration sensitivity from
examples.verify; if the default fixtures change, regenerate/review both reports.
