import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from helper_functions import Helper
helper = Helper()

M_ranges = np.arange(1, 0, -0.1)
def cycles(moment_range: float, a=8.392, b=5.392): return 10**(a-b*moment_range)


df = pd.DataFrame(columns=["M_range", "n", "n_cum"])
n = list()
n_cum = [0]
for M_range in M_ranges:
	n.append(cycles(M_range))
	n_cum.append(n_cum[-1]+n[-1])
n_cum = n_cum[1:]
pd.DataFrame({"M_range": M_ranges, "n": n, "n_cum": n_cum}).to_csv("../data/results/fatigue.dat", index=False)


def load(load_range, frequency, times): return load_range/2*np.sin(2*np.pi*frequency*times)


# fig, ax = plt.subplots()
# time_duration = 10
# time_resolution = 1000
# time = np.linspace(0, time_duration, time_resolution)
# combined_load = np.zeros(time_resolution)
# for M_range, ns in zip(M_ranges, n):
# 	combined_load += load(M_range, ns/6000, time)
# ax.plot(time, combined_load)
# plt.show()

fig, ax = plt.subplots()
ax.plot([1]+n_cum+[n_cum[-1]], [1]+M_ranges.tolist()+[0], label="no wind change")

n_cum = [1e3] + [n+1e3 for n in n_cum]
M_ranges = [2] + M_ranges.tolist()
pd.DataFrame({"M_range": M_ranges, "n": [1e3]+n, "n_cum": n_cum}).to_csv("../data/results/fatigue_wind.dat",
																		 index=False)
n_cum = [1] + n_cum + [n_cum[-1]]
M_ranges = [2] + M_ranges + [M_ranges[-1]]
ax.plot(n_cum, M_ranges, label="with wind change")
helper.handle_axis(ax, x_scale="log", grid=True, legend=True)
helper.handle_figure(fig, show=True, size=(6, 4))
