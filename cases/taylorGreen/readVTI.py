import os
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt

L = 1
kL = 2 * np.pi / L
U0 = 0.05
nu = 0.002

def processVTI(directory, P2_list, K_list):
    vti_files = sorted([f for f in os.listdir(directory) if f.endswith(".vti")])
    for vti_file in vti_files:
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(os.path.join(directory, vti_file))
        reader.Update()
        data = reader.GetOutput()

        vel = vtk_to_numpy(data.GetPointData().GetArray("Velocity")) # shape = (Npoints, 3)
        density = vtk_to_numpy(data.GetPointData().GetArray('Density'))
        
        N = int(np.sqrt(len(density)))
        
        pressure = density / 3.0

        integrand_sum_pressure = 0.0
        integrand_sum_energy = 0.0

        for i in range(N):
            for j in range(N):
            
                x = L / N * (i + 0.5)
                y = L / N * (j + 0.5)

                k = i * N + j

                cos1 = np.cos(kL * (x + y))
                cos2 = np.cos(kL * (x - y))

                integrand_sum_pressure += pressure[k] * cos1 * cos2

                integrand_sum_energy += vel[k][0]**2 + vel[k][1]**2

        P2_list.append(-16.0 / (U0**2 * N**2) * integrand_sum_pressure)
        K_list.append(2.0 / (U0**2 * N**2) * integrand_sum_energy)

P2_iterated = []
P2_non_iterated = []
K_iterated = []
K_non_iterated = []

processVTI("iterated/", P2_iterated, K_iterated)
processVTI("constant/", P2_non_iterated, K_non_iterated)

time = np.linspace(0,200, len(P2_iterated))

plt.plot(time, P2_iterated/P2_iterated[0], label = "iterated")
plt.plot(time, P2_non_iterated/P2_iterated[0], label = "constant")
plt.xlim(0,len(P2_iterated))
plt.ylim(0,2)
plt.legend()
plt.show()

plt.plot(time, K_iterated/K_iterated[0], label = "iterated")
plt.plot(time, K_non_iterated/K_non_iterated[0], label = "constant")
plt.xlim(0,200)
plt.ylim(0,1)
plt.legend()
plt.show()