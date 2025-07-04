import yaml
import numpy as np
from dataclasses import dataclass
from typing import List

@dataclass
class SimulationConfig:
    pw: bool
    control_step: int
    repeat_run: int
    nsteps: int
    temperatures: List[float]
    diffusion_barriers: List[float]
    jump_frequencies: np.ndarray
    distances: List[float]

def load_config(filename="input.yaml") -> SimulationConfig:
    with open(filename, "r") as f:
        data = yaml.safe_load(f)
    p = data["simulation_parameters"]
    return SimulationConfig(
        pw=p["PWKMC"],
        control_step=int(p["control_step"]),
        repeat_run=int(p["repeat_run"]),
        nsteps=int(p["nsteps"]),
        temperatures=p["temperatures"],
        diffusion_barriers=p["diffusion_barriers"],
        jump_frequencies=np.array(p["jump_frequencies"], dtype=np.float64) * 1e12,
        distances=p["distances"]
    )