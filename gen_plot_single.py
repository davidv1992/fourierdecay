import matplotlib.pyplot as plt
import numpy as np

v = np.array(eval(open('plotdata.txt','r').read())).transpose()

fig, ax = plt.subplots(figsize=(8,4.4), layout='constrained')
ax.plot(v[0], v[1])
ax.set_xscale('log')
ax.set_xlabel("$a/m$")
ax.set_ylabel("$\\Gamma(a)/\\Gamma(0)$")

fig.savefig("plot.png", dpi=200)
plt.show()
