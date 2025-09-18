import numpy as np
import matplotlib.pyplot as plt

N = 8
L = 1.0
U0 = 0.05
k0 = 2*np.pi / L
h = L / N
Snu = 1 / 0.506

xi = (np.arange(0, N) + 0.5) * h
yi = (np.arange(0, N) + 0.5) * h
Xi, Yi = np.meshgrid(xi, yi, indexing='xy')  

ux = - U0 * np.cos(k0 * Xi) * np.sin(k0 * Yi)
uy =   U0 * np.sin(k0 * Xi) * np.cos(k0 * Yi)

# --- Analytic ---
dux_dx_analytic =  U0 * k0 * np.sin(k0 * Xi) * np.sin(k0 * Yi)   
duy_dy_analytic = -U0 * k0 * np.sin(k0 * Xi) * np.sin(k0 * Yi)  
pxx_analytic = -2.0 / (3.0 * Snu) * (dux_dx_analytic - duy_dy_analytic)


# --- (cell-centered finite differences) ---
dux_dx_fd = (np.roll(ux, -1, axis=1) - np.roll(ux, 1, axis=1)) / (2*h)  # axis=1 -> x
duy_dy_fd = (np.roll(uy, -1, axis=0) - np.roll(uy, 1, axis=0)) / (2*h)  # axis=0 -> y
pxx_fd = -2.0 / (3.0 * Snu) * (dux_dx_fd - duy_dy_fd)

# ---- Spectral derivative ----
kx = 2*np.pi * np.fft.fftfreq(N, d=h)  
ky = 2*np.pi * np.fft.fftfreq(N, d=h)
kx_grid = np.reshape(kx, (1, N))        
ky_grid = np.reshape(ky, (N, 1))       
ux_hat = np.fft.fft2(ux)   
uy_hat = np.fft.fft2(uy)
dux_dx_spec = np.fft.ifft2(1j * kx_grid * ux_hat).real
duy_dy_spec = np.fft.ifft2(1j * ky_grid * uy_hat).real

pxx_spec = -2.0 / (3.0 * Snu) * (dux_dx_spec - duy_dy_spec)

def report(name, a, b):
    diff = a - b
    RMS = np.sqrt(np.mean(diff**2))
    print(f"{name}: max_abs={np.max(np.abs(diff)):.6e}, RMS = {RMS:.6e}")
    
report("FD   vs Spec    ", pxx_spec, pxx_fd)
report("FD   vs Analytic", pxx_fd, pxx_analytic)
report("Spec vs Analytic", pxx_spec, pxx_analytic)

'''
plt.figure(figsize=(6,5))
plt.streamplot(Xi, Yi, ux, uy, color=umag, cmap="plasma", density=1.2)
plt.colorbar(label="|u|")
plt.show()

plt.figure(figsize=(6,5))
plt.pcolormesh(Xi, Yi, p, shading="auto", cmap="RdBu_r") 
plt.colorbar(label="p")
plt.show()
'''
