import data_handling
from blade_design import BladeApproximation
import matplotlib.pyplot as plt
import pandas as pd
from BEM import BEM
from helper_functions import Helper
from scipy import integrate
import numpy as np
helper = Helper()

rotor_radius = 90
number_of_blades = 3
tip_speed_ratio = 10.58

do = {
    "FAST_to_pandas": False,
    "openFAST_to_FAST": False,
    "scale_rotor": False,
    "modify_rotor": False,
    "modify_thickness_factor": True,
    "optimum_NREL": False,
    "optimum_DTU": False,
    "optimum_IEA": False,
    "plot_optimum_file": False,
    "compare_optimal_actual": False,
    "BEM": False,
    "plot_BEM_results": False,
    "t_distribution": False
}

if do["FAST_to_pandas"]:
    data_handling.FAST_to_pandas(dir_FAST_data="../data/IEA_10MW/airfoils/FAST",
                                 dir_save="../data/IEA_10MW/airfoils/pandas")

if do["openFAST_to_FAST"]:
    data_handling.prepare_openFAST_to_FAST(dir_openFAST_data="../data/openFAST",
                                           aero_dyn_blade_file="IEA-10.0-198-RWT_AeroDyn15_blade.dat",
                                           elasto_dyn_blade_file= "IEA-10.0-198-RWT_ElastoDyn_blade.dat",
                                           dir_FAST="../data/FAST_integration")

if do["scale_rotor"]:
    data_handling.scale_blade_by_R(dir_FAST="../data/FAST_integration", file_type="dat", old_radius=99.155,
                                   new_radius=rotor_radius)

if do["modify_rotor"]:
    data_handling.incorporate_modifications(dir_FAST="../data/FAST_integration", file="modifications_blade_J.dat")

if do["modify_thickness_factor"]:
    add_t = 0.00
    # add_t = "../data/FAST_integration/additional_t.dat"
    data_handling.add_thickness_factor(add_t)

if do["optimum_NREL"]:
    NREL = BladeApproximation(root_dir="../data",
                              blade_dir="NREL_5MW",
                              blade_filename="NREL_5MW_blade_data.txt",
                              save_dir="results",
                              blade_name="NREL_5MW")
    NREL.set_rotor(tip_speed_ratio= tip_speed_ratio,
                   rotor_radius=rotor_radius,
                   number_of_blades=number_of_blades)
    NREL.set_blade_columns(column_positions="RNodes",
                           column_airfoil_path="Airfoil")
    NREL.set_airfoil_columns()
    NREL.chord_and_twist(actual_radius=63, skip_first_percentage=15)

if do["optimum_DTU"]:
    DTU = BladeApproximation(root_dir="../data",
                             blade_dir="DTU_10MW",
                             blade_filename="blade_data_new.txt",
                             save_dir="results",
                             blade_name="DTU_10MW")
    DTU.set_rotor(tip_speed_ratio= tip_speed_ratio,
                  rotor_radius=rotor_radius,
                  number_of_blades=number_of_blades)
    DTU.set_blade_columns(column_positions="radius",
                          column_airfoil_path="airfoil",
                          interpolation_required=True)
    DTU.set_airfoil_columns()
    DTU.chord_and_twist(actual_radius=89.15, skip_first_percentage=15)

if do["optimum_IEA"]:
    IEA = BladeApproximation(root_dir="../data",
                             blade_dir="IEA_10MW",
                             blade_filename="blade_data.txt",
                             save_dir="results",
                             blade_name="IEA_10MW")
    IEA.set_rotor(tip_speed_ratio= tip_speed_ratio,
                  rotor_radius=rotor_radius,
                  number_of_blades=number_of_blades)
    IEA.set_blade_columns(column_positions="BlSpn",
                          column_airfoil_path="BlAFID")
    IEA.set_airfoil_columns()
    IEA.chord_and_twist(actual_radius=99.155, skip_first_percentage=15)

if do["plot_optimum_file"]:
    data_handling.plot_results(file_path="../data/results/optimum_results.dat",
                               plot_dir="../data/results")

if do["compare_optimal_actual"]:
    data_handling.compare_optimal_actual(dir_save="../data/results",
                                         file_optimum="../data/results/optimum_results.dat",
                                         file_actual="../data/IEA_10MW/blade_data.txt", actual_radius=99,
                                         aero_max_radius=90, skip_first_percentage=10)

if do["BEM"]:
    bem = BEM("../data/results")
    bem.set_constants(rotor_radius=90, root_radius=2.4*90/99.155, n_blades=3, air_density=1.225)
    bem.solve_TUD("../data/FAST_integration/blade_aero_dyn_modified.dat", wind_speed=8, tip_speed_ratio=10.58, pitch=0,
                  start_radius=4)
    data_handling.calculate_root_moments("../data/results/BEM_results.dat",
                                         json_file="../data/results/shaft_moments.json", turbine_name="modified")

if do["plot_BEM_results"]:
    root_radius = 2.4*90/99.155
    df_bem_results_new = pd.read_csv("../data/results/BEM_results_modified.dat")
    df_blade_new = pd.read_csv("../data/FAST_integration/blade_aero_dyn_modified.dat")

    df_bem_results_R_scaled = pd.read_csv("../data/results/BEM_results_R_scaled.dat")
    df_blade_original = pd.read_csv("../data/FAST_integration/blade_aero_dyn_R_scaled.dat")
    mu_blade = (root_radius+df_blade_original["BlSpn"])/90

    r = np.array([0, root_radius]+df_bem_results_new["r_centre"].tolist()+[90])
    mu= r/90
    f_t_mod = np.array([0, 0]+df_bem_results_new["f_t"].tolist()+[0])
    f_t_R = np.array([0, 0]+df_bem_results_R_scaled["f_t"].tolist()+[0])
    new_torque = integrate.simpson(f_t_mod*r, r)
    old_torque = integrate.simpson(f_t_R*r, r)
    print("New minus old torque: ", new_torque-old_torque, ". Relative change in torque:", (new_torque-old_torque)/old_torque)


    fig, ax = plt.subplots()
    # all about angles
    # ax[0].plot(df["r_centre"], df["alpha"], label="alpha")
    # ax[0].plot(df["r_centre"], df["alpha_max"], label="alpha max")
    ax.plot(mu[2:-1], df_bem_results_R_scaled["alpha_max"]-df_bem_results_R_scaled["alpha"],label=r"$\Delta\alpha_s$")
    ax.plot(mu[2:-1], df_bem_results_R_scaled["alpha_best"]-df_bem_results_R_scaled["alpha"], label=r"$\Delta\alpha_{max}$")
    helper.handle_axis(ax, x_label=r"$\mu_r$ (-)", y_label=r"angle ($\degree$)", grid=True, legend=True,
                       title=r"Stall margin and distance to $\alpha_{max}$", line_width=5, font_size=35)
    helper.handle_figure(fig, save_to="../data/results/optimisation/dalpha_R.png")

    fig, ax = plt.subplots()
    print(np.min((df_bem_results_new["alpha_max"]-df_bem_results_new["alpha"])[10:]))
    ax.plot(mu[2:-1], df_bem_results_new["alpha_max"]-df_bem_results_new["alpha"],label="stall margin")
    ax.plot(mu[2:-1], df_bem_results_new["alpha_best"]-df_bem_results_new["alpha"], label=r"$\Delta\alpha_{max}$")
    helper.handle_axis(ax, x_label=r"$\mu_r$ (-)", y_label=r"angle ($\degree$)", grid=True, legend=True,
                       title=r"Stall margin and distance to $\alpha_{max}$", line_width=5, font_size=35)
    helper.handle_figure(fig, save_to="../data/results/optimisation/dalpha_mod.png")

    # all about forces
    fig, ax = plt.subplots()
    ax.plot(mu, f_t_R, label="R scaled")
    ax.plot(mu, f_t_mod, label="modified")
    helper.handle_axis(ax, title="Tangential loads", x_label=r"$\mu_r$ (-)", y_label=r"$f_t$ (N/m)", line_width=5,
                       font_size=35, grid=True, legend=True)
    helper.handle_figure(fig, save_to="../data/results/optimisation/f_t.png")


    fig, ax = plt.subplots()
    f_a_mod = np.array([0, 0]+df_bem_results_new["f_n"].tolist()+[0])
    f_a_R = np.array([0, 0]+df_bem_results_R_scaled["f_n"].tolist()+[0])
    ax.plot(mu, f_a_R, label="R scaled")
    ax.plot(mu, f_a_mod, label="modified")
    helper.handle_axis(ax, title="Axial loads", x_label=r"$\mu_r$ (-)", y_label=r"$f_a$ (N/m)", line_width=5,
                       font_size=35, grid=True, legend=True)
    helper.handle_figure(fig, save_to="../data/results/optimisation/f_a.png")

    fig, ax = plt.subplots()
    # change in twist distribution
    ax.plot(mu_blade, df_blade_original["BlTwist"], label="R scaled")
    ax.plot(mu_blade, df_blade_new["BlTwist"], label="modified")
    helper.handle_axis(ax, title=r"Twist $\phi$", x_label=r"$\mu_r$ (-)", y_label="twist (-)", line_width=5,
                       font_size=35, grid=True, legend=True)
    helper.handle_figure(fig, save_to="../data/results/optimisation/twist.png")

    fig, ax = plt.subplots()
    ax.plot(mu[2:-1], df_bem_results_R_scaled["a"], label="R scaled")
    ax.plot(mu[2:-1], df_bem_results_new["a"], label="modified")
    helper.handle_axis(ax, title="Induction factors", x_label=r"$\mu_r$", y_label="a (-)", line_width=5, font_size=35,
                       grid=True, legend=True)
    helper.handle_figure(fig, save_to="../data/results/optimisation/a.png")

    # ax.plot(r, df_bem_results_R_scaled["a_prime"], label="original")
    # ax.plot(r, df_bem_results_new["a_prime"], label="new")
    #
    # helper.handle_axis(ax, x_label="radial position (m)", grid=True, legend=True,
    #                    y_label=["angle in degree", "twist (°), chord (m)", "load in (N/m)", "axial induction",
    #                             "load in (N/m)", "tangential induction"])
    # helper.handle_figure(fig, size=(7,7), close=False)
    #
    # fig, ax = plt.subplots(2)
    # ax[0].plot(r, df_bem_results_R_scaled["inflow_speed"], label="original")
    # ax[0].plot(r, df_bem_results_new["inflow_speed"], label="new")
    #
    # ax[1].plot(r, df_bem_results_R_scaled["phi"], label="original")
    # ax[1].plot(r, df_bem_results_new["phi"], label="new")
    #
    # helper.handle_axis(ax, x_label="radial position", y_label=["inflow velocity in m/s", "inflow angle in degree"],
    #                    legend=True, grid=True)
    # helper.handle_figure(fig, size=(3,5), close=False)
    plt.show()


if do["t_distribution"]:
    BEM_IEA10 = pd.read_csv("../data/results/BEM_results_original.dat")
    BEM_7_mod = pd.read_csv("../data/results/BEM_results_modified.dat")
    BEM_7_R = pd.read_csv("../data/results/BEM_results_R_scaled.dat")

    # f_t_IEA10 = np.append(0, BEM_IEA10["f_t"], 0)
    f_t_7 = np.array([0,0]+BEM_7_mod["f_t"].tolist()+[0]) # the second 0 is needed because the first blade element
    # after the root is skipped due to convergence problems
    f_t_7R = np.array([0,0]+BEM_7_R["f_t"].tolist()+[0])

    root_radius = 2.4*90/99.155
    r_IEA10 = np.append(BEM_IEA10["r_centre"], 99.155)
    r_7_mod = np.array([root_radius, root_radius+3.028133912677172]+BEM_7_mod["r_centre"].tolist()+[90])
    r_7_R = np.array([root_radius, root_radius+3.028133912677172]+BEM_7_R["r_centre"].tolist()+[90])

    R_IEA10 = 99
    R_7 = 90
    print(r_7_R.shape, f_t_7R.shape)
    moment_IEA10 = []
    moment_7_m = []
    moment_7_R = []
    for i in range(r_7_mod.shape[0]-1):
        # moment_IEA10.append(integrate.trapz(f_t_IEA10[::-1][:i+1],-r_IEA10[::-1][:i+1]))
        moment_7_m.append(integrate.trapz(f_t_7[::-1][:i+2],-r_7_R[::-1][:i+2]))
        moment_7_R.append(integrate.trapz(f_t_7R[::-1][:i+2],-r_7_R[::-1][:i+2]))
    t = np.append(1, np.array(moment_7_m)/np.array(moment_7_R))[::-1]

    mu = r_7_mod/90
    # pd.DataFrame({"mu_r": mu, "t":t}).to_csv("../data/results/t.dat", index=False)
    # fig, ax = plt.subplots()
    # ax.plot(mu, t)
    # helper.handle_axis(ax, title="Thickness factor t", x_label=r"$\mu_r$ (-)", y_label="t (-)", grid=True,
    #                    line_width=5, font_size=35)
    # helper.handle_figure(fig, save_to="../data/results/optimisation/t.png")