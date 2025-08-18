import pandas as pd
import matplotlib.pyplot as plt

arquivo = "output/out_1001.tsv"

df = pd.read_csv(arquivo, sep=r"\s+")

y = df["y"]
vx = df["vx"]

plt.plot(y, vx)
plt.xlabel("y")
plt.ylabel("vx")
plt.grid(True)
plt.show()