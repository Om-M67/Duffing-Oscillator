"""
Duffing Oscillator Simulation & Fourier Analysis
==================================================
Numerical solution using 4th-order Runge-Kutta (RK4).
Discrete Fourier Transform (DFT) applied to analyze
frequency components of the oscillator's displacement.
"""

import numpy as np
import matplotlib.pyplot as plt

# ── PARAMETERS ────────────────────────────────────────────────────────────────
# Duffing equation: x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t)
delta = 0.3       # damping coefficient
alpha = -1.0      # linear stiffness (negative → double-well potential)
beta  = 1.0       # nonlinear stiffness
gamma = 0.5       # forcing amplitude
omega = 1.2       # forcing frequency

# Time parameters
t_start = 0.0
t_end   = 200.0
dt      = 0.01
t       = np.arange(t_start, t_end, dt)

# Initial conditions [x0, v0]
x0, v0 = 1.0, 0.0

# ── RK4 SOLVER ────────────────────────────────────────────────────────────────
def duffing(state, t):
    """
    Returns derivatives [dx/dt, dv/dt] for the Duffing equation.
    state = [x, v]  where v = dx/dt
    """
    x, v = state
    dxdt = v
    dvdt = gamma * np.cos(omega * t) - delta * v - alpha * x - beta * x**3
    return np.array([dxdt, dvdt])

def rk4_step(state, t, dt):
    """Single RK4 step."""
    k1 = duffing(state,            t)
    k2 = duffing(state + dt/2 * k1, t + dt/2)
    k3 = duffing(state + dt/2 * k2, t + dt/2)
    k4 = duffing(state + dt    * k3, t + dt)
    return state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

# Integrate
n      = len(t)
states = np.zeros((n, 2))
states[0] = [x0, v0]

for i in range(1, n):
    states[i] = rk4_step(states[i-1], t[i-1], dt)

x = states[:, 0]   # displacement
v = states[:, 1]   # velocity

# ── FOURIER ANALYSIS ──────────────────────────────────────────────────────────
# Use second half of signal (skip transient)
half = n // 2
x_steady = x[half:]
t_steady = t[half:]

# DFT
N       = len(x_steady)
X_fft   = np.fft.fft(x_steady)
freqs   = np.fft.fftfreq(N, d=dt)
power   = np.abs(X_fft) ** 2

# Keep positive frequencies only
pos_mask = freqs > 0
freqs_pos = freqs[pos_mask]
power_pos = power[pos_mask]

# ── PLOTS ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({"figure.dpi": 150, "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False})

fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle("Duffing Oscillator Simulation & Fourier Analysis",
             fontsize=14, fontweight="bold", y=1.01)

# 1 — Time series (full)
ax = axes[0, 0]
ax.plot(t, x, color="#4C72B0", linewidth=0.6)
ax.set_xlabel("Time")
ax.set_ylabel("Displacement x(t)")
ax.set_title("Time Series (full)")

# 2 — Steady-state time series
ax = axes[0, 1]
ax.plot(t_steady, x_steady, color="#55A868", linewidth=0.8)
ax.set_xlabel("Time")
ax.set_ylabel("Displacement x(t)")
ax.set_title("Steady-State Response")

# 3 — Phase portrait
ax = axes[1, 0]
ax.plot(x[half:], v[half:], color="#E05C5C", linewidth=0.5, alpha=0.8)
ax.set_xlabel("Displacement x")
ax.set_ylabel("Velocity v")
ax.set_title("Phase Portrait (steady-state)")

# 4 — Power spectrum
ax = axes[1, 1]
ax.semilogy(freqs_pos, power_pos, color="#8172B2", linewidth=0.8)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Power (log scale)")
ax.set_title("Power Spectrum (DFT)")
ax.set_xlim(0, 2.0)

plt.tight_layout()
plt.savefig("duffing_analysis.png", bbox_inches="tight")
plt.show()
print("Plot saved → duffing_analysis.png")

# ── PRINT KEY RESULTS ─────────────────────────────────────────────────────────
dominant_freq = freqs_pos[np.argmax(power_pos)]
print("\n" + "=" * 50)
print("RESULTS")
print("=" * 50)
print(f"  Forcing frequency (omega/2pi) : {omega / (2*np.pi):.4f} Hz")
print(f"  Dominant frequency (DFT)      : {dominant_freq:.4f} Hz")
print(f"  Max displacement              : {x[half:].max():.4f}")
print(f"  Min displacement              : {x[half:].min():.4f}")
print(f"  Time steps simulated          : {n:,}")
print("=" * 50)
