import data_handling
from blade_design import BladeApproximation
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from BEM import BEM
from helper_functions import Helper
from scipy import integrate
helper = Helper()

rotor_radius = 90
number_of_blades = 3
tip_speed_ratio = 10.58

do = {
    "FAST_to_pandas": False,
    "openFAST_to_FAST": False,
    "NREL": False,
    "DTU": False,
    "IEA": False,
    "plot_results_file": False,
    "BEM": True,
    "plot_BEM_results": True,
}

if do["FAST_to_pandas"]:
    data_handling.FAST_to_pandas(dir_FAST_data="../data/IEA_10MW/airfoils/FAST",
                                 dir_save="../data/IEA_10MW/airfoils/pandas")

if do["openFAST_to_FAST"]:
    data_handling.prepare_openFAST_to_FAST(dir_openFAST_data="../data/openFAST",
                                           aero_dyn_blade_file="IEA-10.0-198-RWT_AeroDyn15_blade.dat",
                                           elasto_dyn_blade_file= "IEA-10.0-198-RWT_ElastoDyn_blade.dat",
                                           dir_FAST="../data/FAST_integration",
                                           incorporate_external={"new_blade_data.txt": ["BlTwist"]})

if do["NREL"]:
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
    NREL.chord_and_twist(skip_first_percentage=15)

if do["DTU"]:
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
    DTU.chord_and_twist(skip_first_percentage=15)

if do["IEA"]:
    IEA = BladeApproximation(root_dir="../data",
                             blade_dir="IEA_10MW",
                             blade_filename="blade_data.txt",
                             save_dir="results",
                             blade_name="IEA_10MW")
    IEA.set_rotor(tip_speed_ratio= 10.77,
                  rotor_radius=rotor_radius,
                  number_of_blades=number_of_blades)
    IEA.set_blade_columns(column_positions="BlSpn",
                          column_airfoil_path="BlAFID")
    IEA.set_airfoil_columns()
    IEA.chord_and_twist(skip_first_percentage=15)

if do["plot_results_file"]:
    data_handling.plot_results(file_path="../data/results/results.dat",
                               plot_dir="../data/results")

if do["BEM"]:
    bem = BEM("../data/results")
    bem.set_constants(rotor_radius=99, root_radius=0, n_blades=3, air_density=1.225)
    bem.solve_TUD("../data/IEA_10MW/blade_data_J.txt", wind_speed=8, tip_speed_ratio=10.58, pitch=0, start_radius=10)

if do["plot_BEM_results"]:
    df_bem_results_new = pd.read_csv("../data/results/BEM_results.dat")
    df_blade_new = pd.read_csv("../data/IEA_10MW/blade_data_J.txt")

    df_bem_results_original = pd.read_csv("../data/results/BEM_results_original.dat")
    df_blade_original = pd.read_csv("../data/IEA_10MW/blade_data.txt")

    new_torque = integrate.simpson(df_bem_results_new["f_t"]*df_bem_results_new["r_centre"], df_bem_results_new["r_centre"])
    old_torque = integrate.simpson(df_bem_results_original["f_t"]*df_bem_results_original["r_centre"], df_bem_results_original["r_centre"])
    print("New minus old torque: ", new_torque-old_torque, ". Relative change in torque:", (new_torque-old_torque)/old_torque)
    fig, ax = plt.subplots(3,2)
    # all about angles
    # ax[0].plot(df["r_centre"], df["alpha"], label="alpha")
    # ax[0].plot(df["r_centre"], df["alpha_max"], label="alpha max")
    ax[0,0].plot(df_bem_results_new["r_centre"], df_bem_results_new["alpha_max"]-df_bem_results_new["alpha"],
                label="stall margin")
    ax[0,0].plot(df_bem_results_new["r_centre"], df_bem_results_new["alpha_best"]-df_bem_results_new["alpha"],
                label="alpha best - alpha")
    # ax[0].plot(df_bem_results_new["r_centre"], df_bem_results_new["alpha_best"], label="alpha best")

    # all about forces
    ax[1,0].plot(df_bem_results_new["r_centre"], df_bem_results_new["f_t"], label="f_t (new)")
    ax[1,0].plot(df_bem_results_original["r_centre"], df_bem_results_original["f_t"], label="f_t (old)")
    ax[2,0].plot(df_bem_results_new["r_centre"], df_bem_results_new["f_n"], label="f_n (new)")
    ax[2,0].plot(df_bem_results_original["r_centre"], df_bem_results_original["f_n"], label="f_n (old)")

    # change in twist distribution
    ax[0,1].plot(df_blade_original["BlSpn"], df_blade_original["BlTwist"], label="original")
    ax[0,1].plot(df_blade_new["BlSpn"], df_blade_new["BlTwist"], label="new")

    ax[1,1].plot(df_bem_results_original["r_centre"], df_bem_results_original["a"], label="original")
    ax[1,1].plot(df_bem_results_new["r_centre"], df_bem_results_new["a"], label="new")

    ax[2,1].plot(df_bem_results_original["r_centre"], df_bem_results_original["a_prime"], label="original")
    ax[2,1].plot(df_bem_results_new["r_centre"], df_bem_results_new["a_prime"], label="new")

    helper.handle_axis(ax, x_label="radial position (m)", grid=True, legend=True,
                       y_label=["angle in degree", "twist in degree", "load in (N/m)", "axial induction","load in (N/m)",
                                "tangential induction"])
    helper.handle_figure(fig, size=(7,7), close=False)

    fig, ax = plt.subplots(2)
    ax[0].plot(df_bem_results_original["r_centre"], df_bem_results_original["inflow_velocity"], label="original")
    ax[0].plot(df_bem_results_new["r_centre"], df_bem_results_new["inflow_velocity"], label="new")

    ax[1].plot(df_bem_results_original["r_centre"], df_bem_results_original["inflow_angle"], label="original")
    ax[1].plot(df_bem_results_new["r_centre"], df_bem_results_new["inflow_angle"], label="new")

    helper.handle_axis(ax, x_label="radial position", y_label=["inflow velocity in m/s", "inflow angle in degree"],
                       legend=True, grid=True)
    helper.handle_figure(fig, size=(3,5), close=False)
    plt.show()
