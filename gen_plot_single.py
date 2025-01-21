import matplotlib
import matplotlib.pyplot as plt
import numpy as np

v = np.array(eval(open('plotdata.txt','r').read())).transpose()

matplotlib.rcParams['font.size'] = 12

fig, ax = plt.subplots(figsize=(8,4.4), layout='constrained')
ax.plot(v[0], v[1])
ax.set_xscale('log')
ax.set_xlabel("$a/m$")
ax.set_ylabel("$\\Gamma(a)/\\Gamma(0)$")
ax.set_xlim(0.03,30)
ax.set_ylim(0.5,5)

fig.savefig("plot.png", dpi=200)
plt.show()
