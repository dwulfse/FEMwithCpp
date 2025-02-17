from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("FEMwithClasses//OOPFEM1D//main//hp_error.csv")

degrees = df["p"].unique()

for p in degrees:
    subset = df[df["p"] == p]
    plt.loglog(subset["h"], subset["error"], label=f"p = {p}")

plt.xlabel("h")
plt.ylabel("error")
plt.title("Error vs. h for different polynomial degrees")
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()

avg_orders = []
for p in degrees:
    subset = df[df["p"] == p]
    errors = subset["error"].to_numpy()
    hs = subset["h"].to_numpy()
    orders = np.log(errors[:-1] / errors[1:]) / np.log(hs[:-1] / hs[1:])
    avg_orders.append(orders.mean())

print("Average orders of convergence:")
for p, avg_order in zip(degrees, avg_orders):
    print(f"p = {p}: {avg_order:.2f}")
