import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_dir = "../../data/NREL_5MW"
blade_file = "NREL_5MW_blade_data.txt"
tsr = 7.55
R = 85
B = 3


df_blade = pd.read_csv(data_dir+"/"+blade_file, index_col=None)
df_blade["r/R"] = df_blade["RNodes"]/df_blade["RNodes"].max()
positions, twist, chord, l2ds = list(), list(), list(), list()
for i, row in df_blade.iterrows():
    if i < 3: continue
    airfoil_file = row["Airfoil"]
    df_airfoil = pd.read_csv(data_dir+"/"+airfoil_file, index_col=None)
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


