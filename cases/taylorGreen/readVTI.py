import os
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt

directory = "outputVTI/"

vti_files = sorted([f for f in os.listdir(directory) if f.endswith(".vti")])

K = []

for vti_file in vti_files:
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(os.path.join(directory, vti_file))
    reader.Update()
    data = reader.GetOutput()

    vel_vtk = data.GetPointData().GetArray("Velocity")
    vel = vtk_to_numpy(vel_vtk)   # shape = (Npoints, 3)

    u2 = np.sum(vel**2, axis=1)
    sum = np.sum(u2)

    K.append(sum)

np.savetxt("K.txt", K)
plt.plot(np.linspace(0,100,len(K)), K/K[0])
plt.grid()
plt.xlim(0, 100)
plt.ylim(0, 1)
plt.show()