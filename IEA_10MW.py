import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FAST_to_pandas = False
if FAST_to_pandas:
    data_dir = "data/IEA_10MW/airfoils/FAST_format"
    file_blank = "IEA-10.0-198-RWT_AeroDyn15_Polar_"
    for i in range(30):
        idx = f"0{i}" if i < 10 else f"{i}"
        file = data_dir+"/"+file_blank+idx+".dat"
        df = pd.read_csv(file, skiprows=54, names=["alpha", "c_l", "c_d", "c_m"], delim_whitespace=True)
        df.to_csv("data/IEA_10MW/airfoils/pandas_format/"+file_blank+f"{i}"+".dat", index=None)


data_dir = "data/IEA_10MW"
blade_file = "blade_data.txt"
file_blank = "IEA-10.0-198-RWT_AeroDyn15_Polar_"
tsr = 7.55
R = 85
B = 3

df_blade = pd.read_csv(data_dir+"/"+blade_file, index_col=None)
df_blade["r/R"] = df_blade["BlSpn"]/df_blade["BlSpn"].max()
positions, twist, chord, l2ds = list(), list(), list(), list()
for i, row in df_blade.iterrows():
    if i < 1: continue
    df_airfoil = pd.read_csv(data_dir+"/airfoils/pandas_format/"+file_blank+f"{i}.dat", index_col=None)
    lift2drag = df_airfoil["c_l"]/df_airfoil["c_d"]
    c_l = df_airfoil["c_l"].iloc[np.argmax(lift2drag)]
    alpha = df_airfoil["alpha"].iloc[np.argmax(lift2drag)]

    positions.append(row["r/R"]*R)
    l2ds.append(lift2drag.iloc[np.argmax(lift2drag)])
    fac = tsr*row["r/R"]*(1+2/(9*tsr**2*row["r/R"]**2))
    chord.append(16*np.pi*R/(9*tsr*c_l*B*np.sqrt(4/9+fac**2)))
    inflow_angle = np.arctan(2/(3*fac))*180/np.pi
    twist.append(inflow_angle-alpha)

twist = [section_twist-twist[-1] for section_twist in twist]
fig, ax = plt.subplots(3,1)
ax[0].plot(positions, chord)
ax[1].plot(positions, twist)
ax[2].plot(positions, l2ds)
print(np.trapz(positions, l2ds))
plt.show()
