from config import load_config
from events import gen_events
from kmc_core import run_kmc_at_temperature
from kmc_core import kmc_iteration
import numpy as np

def main():
    cfg = load_config()
    repeating, events = gen_events(cfg.distances, mode='read')
    events_num = sum(repeating)
    _ = kmc_iteration(1, 1, np.array([[0, 1, 0],[1, 0, 0]]), 1, 0, 1e12, np.array([1e12]))
    for T in cfg.temperatures:
        run_kmc_at_temperature(T, cfg, repeating, events, events_num)

if __name__ == "__main__":
    main()