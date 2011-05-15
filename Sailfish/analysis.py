#!/usr/bin/python

import math
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import bisect, ridder, newton

fname = sys.argv[1]

def load(fname, num):
    try:
        return np.load('%s_%07d.npz' % (fname, num))
    except:
        return np.load('%s_%08d.npz' % (fname, num))

prev_pos = None

kappa = 0.04
const_A = 0.04
gamma = math.sqrt(8.0 * kappa * const_A / 9.0)
tau_liq = 2.5
mu_liq = (tau_liq - 0.5) / 3.0

def measure(data, time_steps):
    phi = data['phi'][1:-1,1:-1,:]

    y = phi[-1,-1,:]
    vx_y = data['v'][0,-2,-2,:]
    max_x = phi.shape[-1] - 1
    x = np.arange(0, phi.shape[-1])
    f = interp1d(x, y)
    f_vel = interp1d(x, vx_y)
    idxs = np.where(y < 0)[0]

    if idxs[0] > 0:
        bub_x0 = ridder(f, max(idxs[0] - 20, 0), min(idxs[0] + 20, max_x))
        bub_x1 = ridder(f, max(idxs[-1] - 20, 0), min(idxs[-1] + 20, max_x))
        bub_len = bub_x1 - bub_x0
        bub_mid_idx = idxs[len(idxs) / 2]
    else:
        end_idx = np.where(idxs - np.roll(idxs, -1) != -1)[0][0]
        bub_x1 = ridder(f, max(idxs[end_idx] - 20, 0),
                        min(idxs[end_idx] + 20, max_x))

        try:
            idx0 = idxs[end_idx+1]
            bub_x0 = ridder(f, max(idx0 - 20, 0), min(idx0 + 20, max_x))
            part2 = max_x+1 - bub_x0
        except IndexError:
            if f(0) * f(20) < 0:
                bub_x0 = ridder(f, max(idx0 - 20, 0), min(idx0 + 20, max_x))
            else:
                bub_x0 = 0.0
            part2 = 0.0

        bub_len = bub_x1 + part2

        if part2 > bub_x1:
            bub_mid_idx = int(bub_x0 + bub_len / 2)
        else:
            bub_mid_idx = int(bub_x1 - bub_len / 2)

    cutplane = phi[:,:,bub_mid_idx]

    y = cutplane[:,-1]
    x = np.arange(0, phi.shape[0])
    f1 = interp1d(x, y)

    r_axis = phi.shape[0] - ridder(f1, 0, phi.shape[0]-1)
    hx, hy = np.mgrid[0:phi.shape[0], 0:phi.shape[1]]

    y = cutplane[hx == hy]
    x = np.arange(0, len(y))
    f2 = interp1d(x, y)

    r_diag = len(y) - ridder(f2, 0, len(y)-1)
    bub_vel = 0.0
    global prev_pos
    if prev_pos is not None:
        if bub_x0 > prev_pos:
            bub_vel = (bub_x0 - prev_pos) / time_steps
        else:
            bub_vel = (bub_x0 + max_x + 1 - prev_pos) / time_steps

    global gamma, mu_liq
    Re = bub_vel * (phi.shape[0]-1) / mu_liq
    Ca = bub_vel * mu_liq / gamma

    prev_pos = bub_x0
    return bub_x0, bub_vel, bub_len, r_axis, r_diag, Re, Ca, f_vel(bub_x0), f_vel(bub_x1)

for i in range(0, 500):
    num = i * 10000
    data = load(fname, num)
    feat = [num]
    feat.extend(measure(data, 10000))
    print ' '.join(str(x) for x in feat)
