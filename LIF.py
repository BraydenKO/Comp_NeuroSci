import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin


# Membrane Equation:
#τm ddtV(t)=EL−V(t)+RI(t)

t_max = 150e-3   # second
dt = 1e-3        # second
#membrane time constant
tau = 20e-3      # second
#leak potential
el = -60e-3      # milivolt
#reset voltage
vr = -70e-3      # milivolt
#firing threshold
vth = -50e-3     # milivolt
#membrane resistence
r = 100e6        # ohm
#mean current imput
i_mean = 25e-11  # ampere


def I():
  i = i_mean * (1 + 0.1 * (t_max / dt)**(0.5) * np.random.random([n,step_end]) * 2 -1)
  #i[:,:20] = i[:,:20]*10
  return i


def V(i, v):
  return v + (dt/tau)*(el - v + r*i)

np.random.seed(1948)

step_end = int(t_max/dt) -1
n = 50
t_range = np.linspace(0, t_max, num=step_end, endpoint=False)
v_n = el * np.ones([n,step_end])
i = I()

for step in range(1,step_end):

  v_n[:, step] = V(i[:,step], v_n[:,step-1])

v_mean = np.mean(v_n, axis = 0)
v_std = np.std(v_n, axis = 0)

fig, axs = plt.subplots(2)

axs[0].plot(t_range, v_n.T, '-', alpha = 0.3, color = "lightgray")
axs[0].plot(t_range, v_n[-1], '-', alpha=0.3, color = "lightgray", label = "V(t)")
axs[0].plot(t_range, v_mean, 'C0', alpha=0.8, label='mean')
axs[0].plot(t_range, v_mean+v_std, 'C7', alpha=0.8)
axs[0].plot(t_range, v_mean-v_std, 'C7', alpha=0.8, label='mean $\pm$ std')
axs[1].plot(t_range, i[:1].T, 'g-', alpha = 0.5)

axs[0].set_title('$V_m$ with random I(t)')
axs[1].set_title('Example Current for simulation 0')
plt.xlabel('time (s)')
axs[0].set_ylabel('$V_m$ (V)')
axs[1].set_ylabel('I (a)')
plt.legend()
plt.show()
