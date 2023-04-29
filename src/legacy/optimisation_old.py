df_bem_results_new = pd.read_csv("../data/results/BEM_results_modified.dat")
df_blade_new = pd.read_csv("../data/FAST_integration/blade_aero_dyn_modified.dat")

df_bem_results_R_scaled = pd.read_csv("../data/results/BEM_results_R_scaled.dat")
df_blade_original = pd.read_csv("../data/FAST_integration/blade_aero_dyn_R_scaled.dat")

new_torque = integrate.simpson(df_bem_results_new["f_t"] * df_bem_results_new["r_centre"],
                               df_bem_results_new["r_centre"])
old_torque = integrate.simpson(df_bem_results_R_scaled["f_t"] * df_bem_results_R_scaled["r_centre"],
                               df_bem_results_R_scaled["r_centre"])
print("New minus old torque: ", new_torque - old_torque, ". Relative change in torque:",
      (new_torque - old_torque) / old_torque)
fig, ax = plt.subplots(3, 2)
# all about angles
# ax[0].plot(df["r_centre"], df["alpha"], label="alpha")
# ax[0].plot(df["r_centre"], df["alpha_max"], label="alpha max")
ax[0, 0].plot(df_bem_results_new["r_centre"], df_bem_results_new["alpha_max"] - df_bem_results_new["alpha"],
              label="stall margin")
ax[0, 0].plot(df_bem_results_new["r_centre"], df_bem_results_new["alpha_best"] - df_bem_results_new["alpha"],
              label="alpha best - alpha")
# ax[0].plot(df_bem_results_new["r_centre"], df_bem_results_new["alpha_best"], label="alpha best")

# all about forces
ax[1, 0].plot(df_bem_results_new["r_centre"], df_bem_results_new["f_t"], label="f_t (new)")
ax[1, 0].plot(df_bem_results_R_scaled["r_centre"], df_bem_results_R_scaled["f_t"], label="f_t (old)")
ax[2, 0].plot(df_bem_results_new["r_centre"], df_bem_results_new["f_n"], label="f_n (new)")
ax[2, 0].plot(df_bem_results_R_scaled["r_centre"], df_bem_results_R_scaled["f_n"], label="f_n (old)")

# change in twist distribution
ax[0, 1].plot(df_blade_original["BlSpn"], df_blade_original["BlTwist"], label="twist original")
ax[0, 1].plot(df_blade_new["BlSpn"], df_blade_new["BlTwist"], label="twist new")
# change in chord distribution
ax[0, 1].plot(df_blade_original["BlSpn"], df_blade_original["BlChord"], label="chord original")
ax[0, 1].plot(df_blade_new["BlSpn"], df_blade_new["BlChord"], label="chord new")

ax[1, 1].plot(df_bem_results_R_scaled["r_centre"], df_bem_results_R_scaled["a"], label="original")
ax[1, 1].plot(df_bem_results_new["r_centre"], df_bem_results_new["a"], label="new")

ax[2, 1].plot(df_bem_results_R_scaled["r_centre"], df_bem_results_R_scaled["a_prime"], label="original")
ax[2, 1].plot(df_bem_results_new["r_centre"], df_bem_results_new["a_prime"], label="new")

helper.handle_axis(ax, x_label="radial position (m)", grid=True, legend=True,
                   y_label=["angle in degree", "twist (°), chord (m)", "load in (N/m)", "axial induction",
                            "load in (N/m)", "tangential induction"])
helper.handle_figure(fig, size=(7, 7), close=False)

fig, ax = plt.subplots(2)
ax[0].plot(df_bem_results_R_scaled["r_centre"], df_bem_results_R_scaled["inflow_speed"], label="original")
ax[0].plot(df_bem_results_new["r_centre"], df_bem_results_new["inflow_speed"], label="new")

ax[1].plot(df_bem_results_R_scaled["r_centre"], df_bem_results_R_scaled["phi"], label="original")
ax[1].plot(df_bem_results_new["r_centre"], df_bem_results_new["phi"], label="new")

helper.handle_axis(ax, x_label="radial position", y_label=["inflow velocity in m/s", "inflow angle in degree"],
                   legend=True, grid=True)
helper.handle_figure(fig, size=(3, 5), close=False)
plt.show()