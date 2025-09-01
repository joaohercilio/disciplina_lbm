import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

tau = 0.8
L_phy = 1.0
nu_phy = 1.0e-3
Umax_phy = 0.005

N_list = [8, 16, 32, 64]
errors = []

plt.figure(figsize=(10, 6))

for N in N_list:
    file = f"outputTSV/{N}.tsv"
    df = pd.read_csv(file, sep=r"\s+")
    
    u_latt = df["vx"][1:-1].values
    
    nu_latt = (1/3)*(tau - 0.5)
    h = L_phy / N 
    delta = nu_latt / nu_phy * h**2
    u_phy = u_latt * h / delta
    
    y_nodes = (np.arange(N) + 0.5) * h
    
    u_ana = 4*Umax_phy * y_nodes * (L_phy - y_nodes) / (L_phy**2)
    
    plt.plot(y_nodes, u_phy, 'o', label=f"N={N}")
    
    L2 = np.sqrt(np.sum((u_ana - u_phy)**2) / N)
    errors.append(L2)
    print(f"N = {N}, L2 Error = {L2:.6e}")

y_continuous = np.linspace(0, L_phy, 200)
u_ana_continuous = 4*Umax_phy * y_continuous * (L_phy - y_continuous) / (L_phy**2)
plt.plot(y_continuous, u_ana_continuous, 'k--', label="Analytical Solution", linewidth=2)

plt.xlabel('y')
plt.ylabel('u')
plt.title('Poiseuille Flow')
plt.legend()
plt.grid(True)
plt.show()

log_N = np.log(N_list)
log_errors = np.log(errors)

coefficients = np.polyfit(log_N, log_errors, 1)
slope = coefficients[0]
intercept = coefficients[1]

plt.figure(figsize=(8, 6))
plt.loglog(N_list, errors, 'o', markersize=8, label='Erro L2')
plt.loglog(N_list, np.exp(intercept) * np.array(N_list)**slope, 'r--', 
           label=f'Slope = {slope:.2f}')
plt.xlabel('N')
plt.ylabel('L2 error')
plt.title(f'Convergence analysis')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.show()