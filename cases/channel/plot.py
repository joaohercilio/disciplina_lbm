import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file = "output/out_100000.tsv"
df = pd.read_csv(file, sep=r"\s+")
vx_numerical = df["vx"]

N = len(vx_numerical)

# Physical units
nu_phy = 0.001
L = 1.0
gx = 1.0
U_phy = L*L*gx / (8.0 * nu_phy)
ReMax = L*U_phy / nu_phy

# Lattice units
tau = 0.8
nu_latt = 1.0 / 3.0 * (tau - 0.5)
h = L / N
delta = nu_latt / nu_phy * h**2
U_latt = U_phy * delta/h
    
y = np.linspace(0, L, 30)
vx_analytical = 4*U_phy/L**2 * y*(L-y)
plt.plot(y, vx_analytical)

y = np.linspace(0, L, N)
vx_numerical = vx_numerical * U_phy/U_latt
plt.plot(y, vx_numerical)

plt.xlabel("y")
plt.ylabel("vx")
plt.grid(True)
plt.show()