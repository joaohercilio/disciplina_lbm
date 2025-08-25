import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

tau = 0.8
L_phy = 1.0
nu_phy = 1.0e-3
Umax_phy = 0.01

N_list = [4, 8, 16, 32, 128]
errors = []

for N in N_list:
    file = f"halfway/{N}.tsv"
    df = pd.read_csv(file, sep=r"\s+")
    u_latt = df["vx"][1:-1]

    nu_latt = (1/3)*(tau - 0.5)
    h = L_phy / N            
    delta = nu_latt / nu_phy * h**2
    u_phy = u_latt * h / delta        

    y_nodes = (np.arange(N) + 0.5) * h
    u_ana = 4*Umax_phy/L_phy**2 * y_nodes*(L_phy - y_nodes)
    
    plt.plot(y_nodes, u_phy, 'o', label = f"LBM {N} nodes")

    L2 = 1/np.sqrt(N) * np.sum( np.sqrt( (u_ana - u_phy)**2 ) )
    errors.append(L2)

h = L_phy / 200            
y_nodes = (np.arange(200) + 0.5) * h
u_ana = 4*Umax_phy/L_phy**2 * y_nodes*(L_phy - y_nodes)
plt.plot(y_nodes, u_ana, label = "Analytic")

plt.legend()
plt.show()

log_Nx = np.log10(N_list)
log_L2 = np.log10(errors)

coef_angular, intercept = np.polyfit(log_Nx, log_L2, 1)
trendline = 10**(intercept) * (N_list)**coef_angular

plt.figure(figsize=(8, 6))
plt.loglog(N_list, errors, 'o', markersize=8,)
plt.loglog(N_list, trendline, 'r--', label=f'Coef angular: {coef_angular:.2f}')
plt.xlabel('Nx', fontsize=12)
plt.ylabel('L2 Error', fontsize=12)
plt.title(f'Coeficiente angular = {coef_angular:.2f}', fontsize=14)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.show()



