from src import data_handling
from src.blade_design import BladeApproximation

rotor_radius = 85
number_of_blades = 3
tip_speed_ratio = 7.55

do = {
    "FAST_to_pandas": False,
    "openFAST_to_FAST": True,
    "NREL": False,
    "DTU": False,
    "IEA": False,
    "plot_results_file": False
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
    IEA.set_rotor(tip_speed_ratio= tip_speed_ratio,
                  rotor_radius=rotor_radius,
                  number_of_blades=number_of_blades)
    IEA.set_blade_columns(column_positions="BlSpn",
                          column_airfoil_path="BlAFID")
    IEA.set_airfoil_columns()
    IEA.chord_and_twist(skip_first_percentage=15)

if do["plot_results_file"]:
    data_handling.plot_results(file_path="../data/results/results.dat",
                               plot_dir="../data/results")