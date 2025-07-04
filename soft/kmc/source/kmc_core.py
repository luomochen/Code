import numpy as np
import pandas as pd
from numba import njit
from scipy.constants import Boltzmann
from ase.io import read
import time, copy, multiprocessing as mp
from events import list_extend

@njit
def kmc(events, events_num, ini, total_rate, rates):
    event = events[ini]
    rho1 = np.random.random()
    acc = 0.0
    for i in range(events_num):
        acc += rates[i]
        if rho1 * total_rate < acc:
            fin = event[1 + i]
            event[1 + events_num + i] += 1
            break
    events[ini] = event
    dt = -np.log(np.random.random()) / total_rate
    return dt, fin, events

@njit
def kmc_iteration(t_limit, nsteps, events, events_num, ini, total_rate, rates):
    t, steps = 0.0, 0
    while steps < nsteps:
        dt, ini, events = kmc(events, events_num, ini, total_rate, rates)
        t += dt
        steps += 1
        if t > t_limit:
            break
    return t, steps, events

def weighted_sampling(rates, total_rate, control_step):
    probs = rates / total_rate
    weights = np.ones_like(rates)
    sorted_unique = sorted(set(probs), reverse=True)
    cutoff = sorted_unique[control_step]
    sum_small = sum(p for p in probs if p < cutoff)
    for idx, p in enumerate(probs):
        if p >= cutoff:
            weights[idx] *= round(p / sum_small)
    return weights

def compute_rates(T, barriers, freqs, pw, control_step):
    beta = 1 / (Boltzmann / 1.60218e-19 * T)
    rates = freqs * np.exp(-np.array(barriers) * beta)
    total = np.sum(rates)
    if pw:
        weights = weighted_sampling(rates, total, control_step)
        rates /= weights
        total = np.sum(rates)
    else:
        weights = np.ones_like(rates)
    df = pd.DataFrame({"rate": rates, "probability": rates / total, "weights": weights})
    df.to_csv(f"{int(T)}_rate.csv", index=False)
    return total, rates

def cal_displacement(events, events_num):
    sites = read("./sites.vasp")
    disp = np.zeros(3)
    for e in events:
        for i in range(events_num):
            vec = sites.get_distance(e[0], e[1 + i], mic=True, vector=True)
            disp += vec * e[1 + events_num + i]
    return disp

def kmc_loop(queue, *args):
    t, steps, events = kmc_iteration(*args)
    queue.put((t, steps, events))

def run_kmc_at_temperature(temperature, cfg, repeating_num, events_origin, events_num):
    print(f"Running {temperature} K...")
    start = time.time()
    barriers = list_extend(cfg.diffusion_barriers, repeating_num)
    freqs = list_extend(cfg.jump_frequencies, repeating_num)
    total_rate, rates = compute_rates(temperature, barriers, freqs, cfg.pw, cfg.control_step)

    queue = mp.Manager().Queue()
    pool = mp.Pool()
    for _ in range(cfg.repeat_run):
        events = copy.deepcopy(events_origin)
        ini = np.random.randint(0, len(events))
        pool.apply_async(kmc_loop, (queue, 1e30, cfg.nsteps, events, events_num, ini, total_rate, rates))
    pool.close()
    pool.join()
    queue.put(None)

    records = []
    while True:
        item = queue.get()
        if item is None:
            break
        t, steps, ev = item
        dx, dy, dz = cal_displacement(ev, events_num)
        records.append([t, steps, dx, dy, dz])
    df = pd.DataFrame(records, columns=["t", "steps", "d_x", "d_y", "d_z"])
    df.to_csv(f"{int(temperature)}.csv", index=False)
    print(f"{temperature} K done in {round(time.time()-start, 2)} s")