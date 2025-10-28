#!/usr/bin/python3
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    prog = 'gen_plot_single',
    description = 'Plot a single decay rate vs acceleration',
)
parser.add_argument('filename')
parser.add_argument('aprop', type=float)
parser.add_argument('--amin', type=float, default=0.03)
parser.add_argument('--amax', type=float, default=3.0)
parser.add_argument('--rmin', type=float, default=0.5)
parser.add_argument('--rmax', type=float, default=5)
parser.add_argument('--output', '-o', default='plot.png')

args = parser.parse_args()

v = np.array(eval(open(args.filename,'r').read())).transpose()

matplotlib.rcParams['font.size'] = 12

fig, ax = plt.subplots(figsize=(8,4.4), layout='constrained')
ax.plot(v[0]*args.aprop, v[1])
ax.set_xscale('log')
ax.set_xlabel("$a (m/s^2)$")
ax.set_ylabel("$\\Gamma(a)/\\Gamma(0)$")
ax.set_xlim(args.amin*args.aprop, args.amax*args.aprop)
ax.set_ylim(args.rmin, args.rmax)

fig.savefig(args.output, dpi=200)
plt.show()
