import data_handling
from blade_design import BladeApproximation
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from BEM import BEM
from helper_functions import Helper
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
    "BEM": False,
    "plot_BEM_results": True,
}

if do["FAST_to_pandas"]:
    data_handling.FAST_to_pandas(dir_FAST_data="../data/IEA_10MW/airfoils/FAST",
                                 dir_save="../data/IEA_10MW/airfoils/pandas")

if do["openFAST_to_FAST"]:
    data_handling.prepare_openFAST_to_FAST(dir_openFAST_data="../data/openFAST",
                                           aero_dyn_blade_file="IEA-10.0-198-RWT_AeroDyn15_blade.dat",
                                           elasto_dyn_blade_file= "IEA-10.0-198-RWT_ElastoDyn_blade.dat",
                                           dir_FAST="../data/FAST_integration")

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
    bem.solve_TUD("../data/IEA_10MW/blade_data_S.txt", wind_speed=8, tip_speed_ratio=10.58, pitch=0, start_radius=10)

if do["plot_BEM_results"]:
    df = pd.read_csv("../data/results/BEM_results.dat")
    df_original = pd.read_csv("../data/IEA_10MW/blade_data.txt")
    df_new = pd.read_csv("../data/IEA_10MW/blade_data_S.txt")
    fig, ax = plt.subplots(3,1)
    # all about angles
    # ax[0].plot(df["r_centre"], df["alpha"], label="alpha")
    # ax[0].plot(df["r_centre"], df["alpha_max"], label="alpha max")
    ax[0].plot(df["r_centre"], df["alpha_max"]-df["alpha"], label="stall margin")
    ax[0].plot(df["r_centre"], df["alpha_best"]-df["alpha"], label="alpha best - alpha")
    # ax[0].plot(df["r_centre"], df["alpha_best"], label="alpha best")

    # all about forces
    ax[1].plot(df["r_centre"], df["f_t"], label="f_n")
    ax[1].plot(df["r_centre"], df["f_n"], label="f_t")

    # change in twist distribution
    ax[2].plot(df_original["BlSpn"], df_original["BlTwist"], label="original")
    ax[2].plot(df_new["BlSpn"], df_new["BlTwist"], label="new")
    helper.handle_axis(ax, x_label="radial position (m)", grid=True, legend=True,
                       y_label=["angle in degree", "load in (N/m)","twist in degree"])
    helper.handle_figure(fig, size=(5,7), show=True)
