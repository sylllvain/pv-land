#!/usr/bin/env python3
"""Compare two pv_land_output CSV files column by column (common timesteps).

Usage: compare_outputs.py A.csv B.csv [--meas]
Prints, for each common column, the mean and max |A-B|; with --meas also
the RMSE of both runs against the measured H, LE, Rn, G.
"""
import csv
import math
import sys


def load(p):
    r = list(csv.reader(open(p)))
    h = r[0]
    rows = {}
    for l in r[1:]:
        if l:
            v = [float(x) for x in l]
            rows[int(v[0])] = v
    return h, rows


ha, a = load(sys.argv[1])
hb, b = load(sys.argv[2])
steps = sorted(set(a) & set(b))
print(f'{len(steps)} common timesteps')
print(f"{'column':20s} {'mean|A-B|':>11s} {'max|A-B|':>11s} {'mean A':>10s} {'mean B':>10s}")
for n in ha[1:]:
    if n not in hb:
        continue
    ia, ib = ha.index(n), hb.index(n)
    d = [abs(a[s][ia] - b[s][ib]) for s in steps]
    print(f"{n:20s} {sum(d)/len(d):11.4g} {max(d):11.4g} "
          f"{sum(a[s][ia] for s in steps)/len(steps):10.4g} {sum(b[s][ib] for s in steps)/len(steps):10.4g}")
if '--meas' in sys.argv:
    print()
    for q in ['H', 'LE', 'Rn', 'G']:
        im, iA, iB = ha.index(q + '_meas'), ha.index(q + '_model'), hb.index(q + '_model')
        ok = [s for s in steps if abs(a[s][im]) < 9000]
        ra = math.sqrt(sum((a[s][iA] - a[s][im])**2 for s in ok) / len(ok))
        rb = math.sqrt(sum((b[s][iB] - a[s][im])**2 for s in ok) / len(ok))
        print(f'RMSE {q:3s} vs measurements: A = {ra:7.2f}   B = {rb:7.2f} W/m2')
