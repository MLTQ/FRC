"""Command-line entry points; no notebook or GUI session is required."""
import argparse
from pathlib import Path
from dataclasses import replace
import json
from .config import Experiment
from .results import Result
from .simulation import run
from .experiments import compare, sweep


def main():
    parser = argparse.ArgumentParser(description="Virtual plasma bench · SI units · stages 0 and 1")
    sub = parser.add_subparsers(dest="command", required=True)
    init = sub.add_parser("init", help="Write the complete default experiment configuration")
    init.add_argument("path", type=Path)
    for name in ("run", "compare", "sweep"):
        p = sub.add_parser(name)
        p.add_argument("--config", type=Path)
        p.add_argument("--out", type=Path, default=Path("outputs")/name)
        if name == "run":
            p.add_argument("--stage", type=int, choices=(0, 1, 2))
            p.add_argument("--viewer", action="store_true")
        if name == "compare":
            p.add_argument("--peak-field", type=float, default=.015, help="Peak absolute vacuum center perturbation [T]")
    p = sub.add_parser("view", help="Render a saved archive into a standalone HTML time viewer")
    p.add_argument("directory", type=Path)
    p.add_argument("--gif", action="store_true")
    args = parser.parse_args()
    if args.command == "init":
        args.path.parent.mkdir(parents=True, exist_ok=True)
        Experiment().save(args.path)
        print(args.path)
        return
    # Plotting backend is selected only for CLI output; library users retain control.
    import matplotlib
    matplotlib.use("Agg")
    from .visualization import snapshot, viewer, comparison_plot, phase_diagram
    if args.command == "view":
        result = Result.load(args.directory)
        viewer(result, args.directory/"viewer.html", gif=args.directory/"animation.gif" if args.gif else None)
        print(args.directory/"viewer.html")
        return
    config = Experiment.load(args.config) if args.config else Experiment()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.command == "run":
        if args.stage is not None:
            config = replace(config, stage=args.stage)
        result = run(config)
        result.save(args.out)
        snapshot(result, args.out/"snapshot.png")
        if args.viewer:
            viewer(result, args.out/"viewer.html")
        print(json.dumps(result.metadata["summary"], indent=2))
        for warning in result.metadata["warnings"]:
            print("MODEL NOTE:", warning)
    elif args.command == "compare":
        results = compare(args.out, base=config, peak_field=args.peak_field)
        comparison_plot(results, args.out/"comparison.png")
        for name, result in results.items():
            snapshot(result, args.out/name/"snapshot.png")
        print(json.dumps({name: result.metadata["summary"] for name, result in results.items()}, indent=2))
    elif args.command == "sweep":
        rows = sweep(args.out, base=config)
        phase_diagram(rows, args.out/"phase_diagram.png")
        print(f"Saved {len(rows)} cases to {args.out}")


if __name__ == "__main__":
    main()
