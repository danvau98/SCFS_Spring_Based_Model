from pathlib import Path
import json
import subprocess
import sys

import numpy as np

BASE_DIR = Path(__file__).resolve().parent
EXPERIMENTAL_DIR = BASE_DIR.parent / "single_simulation_OFF"

if BASE_DIR.name == "single_simulation_OFF":
    FORCE_FILE = "SimulationDataForce_OFF.txt"
    EXPERIMENTAL_FILE = "avg_OFF_all_locs.txt"
elif BASE_DIR.name == "single_simulation_ON":
    FORCE_FILE = "SimulationDataForce_ON.txt"
    EXPERIMENTAL_FILE = "avg_ON_all_locs.txt"
else:
    raise ValueError(
        "Place this script in single_simulation_OFF/ "
        "or single_simulation_ON/."
    )


def simulation(parameters):
    args = [
        str(BASE_DIR / "Subprocess"),
        *["%.16g" % value for value in parameters],
        "1.0",  # div
    ]

    completed = subprocess.run(
        args,
        cwd=BASE_DIR,
        capture_output=True,
        text=True,
    )
    
    if completed.returncode != 0:
        raise RuntimeError(
            f"C++ executable failed with code {completed.returncode}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )

    force = np.loadtxt(
        BASE_DIR / FORCE_FILE, dtype=float
    ).reshape(-1)

    experimental = np.loadtxt(
        EXPERIMENTAL_DIR / EXPERIMENTAL_FILE,
        dtype=float,
        ndmin=2,
    )[:, 1]  # Second column: force

    if force.size == 0 or force.shape != experimental.shape:
        raise ValueError(
            f"Simulated and experimental arrays differ: "
            f"{force.shape} vs {experimental.shape}"
        )

    residuals = force - experimental

    if not np.all(np.isfinite(residuals)):
        raise ValueError("Non-finite simulated or experimental values.")

    return {
        "mse": float(np.mean(residuals**2)),
        "mae": float(np.mean(np.abs(residuals))),
        "n": int(force.size),
    }


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(
            "Usage: python single_sim.py pu beta hat_eta hat_theta"
        )

    parameters = np.asarray(sys.argv[1:], dtype=float)
    print(json.dumps(simulation(parameters)))