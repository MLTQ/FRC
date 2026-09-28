"""Scriptable two-coil vacuum/diffusion experiment; run from project root."""
from dataclasses import replace
from pathlib import Path
from plasma_bench.config import Experiment
from plasma_bench.simulation import run
from plasma_bench.visualization import snapshot


def main():
    for stage in (0, 1):
        config = replace(Experiment(), stage=stage)
        output = Path("outputs")/f"example_stage{stage}"
        result = run(config)
        result.save(output)
        snapshot(result, output/"snapshot.png")
        print(stage, result.metadata["summary"])


if __name__ == "__main__":
    main()
