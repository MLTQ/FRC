# Standalone trajectory viewer assembly

Validates the complete archived display-data SHA and embeds it with local
HTML/renderer/controller assets into one offline HTML document. Saves a separate
presentation manifest with hashes for every asset, the data and final HTML.
There is no simulation call, external runtime, CDN or data upload.

The optional localhost-only server serves the new output folder. Its only
write endpoint accepts a bounded PNG with valid signature at a fixed poster.png
path. The historical FRC poster and every scientific source remain untouched.
The output HTML can be reopened without a server; then PNG uses a browser
download. A server is only needed for direct saving beside the HTML.

CLI: `python -m examples.trajectory_display [--serve] [--port 8786]`.
UI assets may be refined before final visual inspection; the assembly manifest
records their exact final hashes separately from frozen data provenance.
