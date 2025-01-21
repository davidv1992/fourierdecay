import matplotlib
import matplotlib.pyplot as plt
import numpy as np

v1 = np.array(eval(open('plotdata1.txt','r').read())).transpose()
v2 = np.array(eval(open('plotdata2.txt','r').read())).transpose()

matplotlib.rcParams['font.size'] = 12

fig, ax = plt.subplots(figsize=(8,4.4), layout='constrained')
ax.plot(v2[0], v2[1], label='Analytic expression', color='C1')
ax.plot(v1[0], v1[1], label='Numerical integration', linestyle=':', color='k')
ax.set_xscale('log')
ax.set_xlabel("$a/m$")
ax.set_ylabel("$\\Gamma(a)/\\Gamma(0)$")
ax.set_xlim(0.03,3)
ax.set_ylim(0,2)
ax.legend()

fig.savefig("plot.png", dpi=200)
plt.show()
