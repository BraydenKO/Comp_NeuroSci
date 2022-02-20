import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.xmargin'] = 0

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

class LIFNeurons:
  """
  Keeps track of membrane potential for multiple realizations of LIF neuron,
  and performs single step discrete time integration.
  """
  def __init__(self, n, t_ref_mu=0.01, t_ref_sigma=0.002,
               tau=20e-3, el=-60e-3, vr=-70e-3, vth=-50e-3, r=100e6):

    # Neuron count
    self.n = n

    # Neuron parameters
    self.tau = tau        # second
    self.el = el          # milivolt
    self.vr = vr          # milivolt
    self.vth = vth        # milivolt
    self.r = r            # ohm

    # Initializes refractory period distribution
    self.t_ref_mu = t_ref_mu
    self.t_ref_sigma = t_ref_sigma
    self.t_ref = self.t_ref_mu + self.t_ref_sigma * np.random.normal(size=self.n)
    self.t_ref[self.t_ref<0] = 0

    # State variables
    self.v = self.el * np.ones(self.n)
    self.spiked = self.v >= self.vth
    self.last_spike = -self.t_ref * np.ones([self.n])
    self.t = 0.
    self.steps = 0

  def V(self):
    return self.v + dt / self.tau * (self.el - self.v + self.r * i)

  def ode_step(self, dt, i):

    # Update running time and steps
    self.t += dt
    self.steps += 1

    # One step of discrete time integration of dt
    self.v = V()

    # Spike and clamp
    self.spiked = (self.v >= self.vth)
    self.v[self.spiked] = self.vr
    self.last_spike[self.spiked] = self.t
    clamped = (self.t_ref > self.t-self.last_spike)
    self.v[clamped] = self.vr

    self.last_spike[self.spiked] = self.t

# Initialize step_end, t_range, n, v_n and i
t_range = np.arange(0, t_max, dt)
step_end = len(t_range)
n = 500
v_n = el * np.ones([n, step_end])
i = I()

# Initialize binary numpy array for raster plot
raster = np.zeros([n,step_end])

# Initialize neurons
neurons = LIFNeurons(n)

# Loop over time steps
for step, t in enumerate(t_range):

  # Call ode_step method
  neurons.ode_step(dt, i[:,step])

  # Log v_n and spike history
  v_n[:,step] = neurons.v
  raster[neurons.spiked, step] = 1.

v_mean = np.mean(v_n, axis=0)
spikes_mean =  raster.sum(axis=0)/ n

# Plot simulations and sample mean
fig, axs = plt.subplots(3)

for v in v_n:  
  axs[0].scatter(t_range, v, color = 'k', marker = '.', alpha = 0.01)
axs[0].plot(t_range, v_mean, 'C1', alpha=0.8, linewidth = 3)
axs[0].set_ylabel('$V_m$ (V)')

# for each neuron j: collect spike times and plot them at height j
axs[1].imshow(raster, cmap='Greys', origin='lower', aspect='auto', extent=[0,t_max, 0, n])
axs[1].set_ylabel('neuron')

# Plot firing rate
axs[2].plot(t_range, spikes_mean)
axs[2].set_xlabel('time (s)')
axs[2].set_ylabel('rate (Hz)')

plt.tight_layout()
plt.show()