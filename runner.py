import matplotlib.pyplot as plt
import numpy as np
from gotran import load_cell, ODESolver
import time

start = time.time()
cell = load_cell('PhePol.ode')
solver = ODESolver(cell, methods='cvode')

t, y = solver.solve(np.arange(0, 2000, 1))

print(" -------------- %s seconds --------------" % (time.time() - start))

#fig, ax = plt.subplots(7, 4)
#axs = ax.flatten()
#for i, s in enumerate(cell.state_symbols):
#    axs[i].plot(t, y[:, solver.module.state_indices(s)])
#    axs[i].set_title(s)

for i, s in enumerate(cell.state_symbols):
    if s=="Cai":
        plt.subplots(1,1)
        plt.plot(t,y[:, solver.module.state_indices(s)])
        plt.title(s)

plt.show()