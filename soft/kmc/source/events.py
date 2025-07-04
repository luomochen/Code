import numpy as np
import pandas as pd
from ase.io import read

def list_extend(lst, repeating_num):
    return [v for i, v in enumerate(lst) for _ in range(repeating_num[i])]

def gen_events(select_path, mode=True):
    if mode is True:
        sites = read("./sites.vasp")
        n_atoms = len(sites)
        events = []
        for i in range(n_atoms):
            event = [i]
            repeating = []
            for path in select_path:
                targets = [j for j in range(n_atoms)
                           if abs(sites.get_distance(i, j, mic=True) - path) < 5E-4]
                event.extend(targets)
                repeating.append(len(targets))
            event.extend([0] * sum(repeating))
            events.append(event)
        pd.DataFrame(events).to_csv("events.csv", index=False)
        pd.DataFrame(repeating).to_csv("repeating_num.csv", index=False)
    else:
        events = pd.read_csv("events.csv").values.tolist()
        repeating = pd.read_csv("repeating_num.csv").values.flatten().astype(int).tolist()
    return repeating, np.array(events)