import numpy as np
import matplotlib.pyplot as plt


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

# Set random number generator
np.random.seed(1948)

def I():
  i = i_mean * (1 + 0.1 * (t_max / dt)**(0.5) * np.random.random([n,step_end]) * 2 -1)
  return i
def V(i, v):
  return v + (dt/tau)*(el - v + r*i)

# Initialize t_range, step_end, n, v_n, i and nbins
t_range = np.arange(0, t_max, dt)
step_end = len(t_range)
n = 500
v_n = el * np.ones([n, step_end])
nbins = 32
i = I()



# Loop over time steps
for step, t in enumerate(t_range):

  # Skip first iteration
  if step==0:
    continue

  # Compute v_n
  v_n[:, step] = V(i[:, step], v_n[:, step - 1])

# Initialize the figure
fig, axs = plt.subplots(1,2)

plt.ylabel('$V_m$ (V)')
axs[0].set_xlabel('time (s)')
axs[1].set_xlabel('Frequency')

# Plot a histogram at t_max/10 (add labels and parameters histtype='stepfilled' and linewidth=0)
axs[1].hist(v_n[:,int(step_end / 10)], nbins, histtype='stepfilled', linewidth=0,label = f't= {t_max/10}s', orientation= 'horizontal')

# Plot a histogram at t_max (add labels and parameters histtype='stepfilled' and linewidth=0)
axs[1].hist(v_n[:,-1], nbins, histtype='stepfilled', linewidth=0,label = f't= {t_max}s', orientation= 'horizontal')

axs[0].plot(t_range, v_n.T, color = 'lightgray')

# Add legend
plt.legend()
plt.show()