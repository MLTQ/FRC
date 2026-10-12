"""Move the memory legend clear of the execution-budget line."""
from matplotlib.axes import Axes
from . import workspace_refinement_report as report
from .implicit_plasma_records import provenance


def run(destination):
    saved = Axes.legend, report.provenance
    def legend(axis, *args, **kwargs):
        if axis.get_ylabel() == 'Own-process peak RSS (GiB)':
            kwargs.update(loc='upper center', bbox_to_anchor=(.5, 1.19), ncol=3)
        return saved[0](axis, *args, **kwargs)
    try:
        Axes.legend = legend
        report.provenance = lambda seeds: provenance([*seeds, 'examples.workspace_refinement_report_layout',
                                                       'examples.workspace_refinement_report'])
        return report.run(destination)
    finally:
        Axes.legend, report.provenance = saved
