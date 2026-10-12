"""Preserve the archived figure calculation while clearing overlapping legends."""
from pathlib import Path
from . import curl_spectrum_report as report
from .implicit_plasma_records import provenance


def run(destination, curl, spectra):
    original_close, original_provenance = report.plt.close, report.provenance
    def close(figure=None):
        if hasattr(figure, 'axes') and len(figure.axes) == 4:
            figure.axes[2].set_ylim(0, 9)
            figure.axes[3].set_ylim(7e-13, 1e-10)
            figure.axes[3].legend(fontsize=8, loc='upper center')
            figure.savefig(Path(destination)/'verification.png', dpi=160)
        return original_close(figure)
    try:
        report.plt.close = close
        report.provenance = lambda seeds: provenance([*seeds, 'examples.curl_spectrum_report_layout'])
        return report.run(destination, curl, spectra)
    finally:
        report.plt.close, report.provenance = original_close, original_provenance
