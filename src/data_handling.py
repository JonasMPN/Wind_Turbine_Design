import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from scipy.integrate import simpson
import numpy as np
import pandas as pd
import os, shutil
from src.helper_functions import Helper
from copy import copy
import json
helper = Helper()


class AirfoilInterpolator:
    def __init__(self,
                 dir_data: str,
                 data_profiles: list or str,
                 **pandas_read_csv):
        if type(data_profiles) == str:
            data_profiles = [data_profiles]
        self.df = pd.DataFrame()
        for file in data_profiles:
            self.df = pd.concat([self.df, pd.read_csv(dir_data+"/"+file, **pandas_read_csv)], ignore_index=True)
        self.interpolator = dict()

    def get_df_profiles(self):
        return self.df

    def prepare_interpolation(self,
                              to_interpolate: dict,
                              constraints: dict = None) -> None:
        """
        Multidimensional interpolation. The data is taken from the files that are specified during the initialisation of
        the class' instance. The parameter 'to_interpolate' states which parameter(s) is(are) to be interpolated and on
        which arguments this interpolation is based. The keys of 'to_interpolate' specify the parameters that will be
        interpolated. The values of each key specify on which arguments the interpolation is based upon.
        :param to_interpolate: keys: values to interpolate, values: arguments for interpolation
        :param constraints: dict of dicts. First key is parameter to be interpolated, value of that key is again a
        dict. That dict has the arguments for the interpolation as keys and the ranges that should be used as value (
        as tuple)
        :return: None
        """
        constraints = constraints if type(constraints) == dict else {}
        for to_inter, arguments in to_interpolate.items():
            if type(arguments) == str:
                arguments = [arguments]
            points = {var: list() for var in arguments}
            values = list()
            df_data_points = self.df

            if to_inter in constraints.keys():
                df_data_points = pd.DataFrame(columns=df_data_points.columns.tolist())
                for arg_parameter, constraint_values in constraints[to_inter].items():
                    for constraint_value in constraint_values:
                        df_data_points = pd.concat([df_data_points,
                                                   self.df[self.df[arg_parameter]==constraint_value]],
                                                   ignore_index=True)
            for _, row in df_data_points.iterrows():
                for arg in arguments:
                    points[arg].append(row[arg])
                values.append(row[to_inter])
            points = np.asarray([all_points for all_points in points.values()]).T
            if points.shape[1] == 1:
                self.interpolator[to_inter] = interpolate.interp1d(points.flatten(), values)
            else:
                self.interpolator[to_inter] = interpolate.LinearNDInterpolator(points, values)

    def __getitem__(self, item):
        return self.interpolator[item]


def plot_results(file_path: str, plot_dir: str) -> None:
    df = pd.read_csv(file_path)
    names = df["name"].unique()
    fig, axes = plt.subplots(3,1)
    for name in names:
        df_tmp = df[df["name"]==name]
        axes[0].plot(df_tmp["radius"], df_tmp["chord"], label=name)
        axes[1].plot(df_tmp["radius"], df_tmp["twist"])
        axes[2].plot(df_tmp["radius"], df_tmp["l2d"])
    helper.handle_axis(axes,
                       title=["NREL 5MW, DTU 10MW, and IEA 10MW in comparison"], x_label="Radius (m)",
                       y_label=["chord (m)", "twist (°)", r"$c_l/c_d$ (-)"], legend=True, font_size=20,
                       line_width=4, grid=True)
    helper.handle_figure(fig, save_to=plot_dir+"/"+f"{names}"+".png", size=(18.5, 14))


def FAST_to_pandas(dir_FAST_data: str,
                   dir_save: str):
    file_blank = "IEA-10.0-198-RWT_AeroDyn15_Polar_"
    for i in range(30):
        idx = f"{i}" if i>9 else f"0{i}"
        file = dir_FAST_data+"/"+file_blank+idx+".dat"
        df = pd.read_csv(file, names=["alpha", "c_l", "c_d", "c_m"], skiprows=54, delim_whitespace=True)
        df.to_csv(dir_save+"/"+file_blank+f"{i}"+".dat", index=False)


def prepare_openFAST_to_FAST(dir_openFAST_data: str,
                             aero_dyn_blade_file: str,
                             elasto_dyn_blade_file: str,
                             dir_FAST: str,
                             file_type: str="dat") -> None:
    """
    This function copies the coordinates and polars directory while changing the filetype of the directories' files
    to "file_type".
    It also writes out the stall angles, c_n slope, and critical c_n and saves it to a file called
    "airfoil_additional_information" (file_type = dat)
    :param dir_openFAST_data:
    :param dir_coordinates:
    :param dir_polars:
    :param dir_FAST:
    :return:
    """
    helper.create_dir(dir_FAST)
    helper.create_dir(dir_FAST+"/coordinates", overwrite=True)
    helper.create_dir(dir_FAST+"/polars", overwrite=True)
    lines_to_use = {
        "StallAngle1": 18,
        "StallAngle2": 19,
        "CnSlope": 21,
        "CritCn1": 36,
        "CritCn2": 37}
    airfoil_additional_data = {col: list() for col in ["polar_file_name", "coord_file_name", *lines_to_use]}
    for file in os.listdir(dir_openFAST_data+"/Airfoils"):
        new_file = file[::-1][file[::-1].find("."):][::-1]+file_type
        if "Polar" in file:
            df = pd.read_csv(dir_openFAST_data+"/Airfoils/"+file, names=["alpha", "c_l", "c_d", "c_m"], skiprows=54,
                             delim_whitespace=True)
            df.to_csv(dir_FAST+"/polars/"+new_file, index=False)
            with open(dir_openFAST_data+"/Airfoils/"+file, "r") as f:
                lines = f.readlines()
                for column, line_number in lines_to_use.items():
                    value = str()
                    for char in lines[line_number]:
                        if char != " ":
                            value += char
                        else:
                            break
                    airfoil_additional_data[column].append(float(value))
        elif "Coords" in file:
            shutil.copyfile(dir_openFAST_data+"/Airfoils/"+file, dir_FAST+"/coordinates/"+new_file)

    airfoil_additional_data["polar_file_name"] = os.listdir(dir_FAST + "/polars")
    airfoil_additional_data["coord_file_name"] = os.listdir(dir_FAST + "/coordinates")
    pd.DataFrame(airfoil_additional_data).to_csv(dir_FAST + f"/airfoil_additional_information.{file_type}", index=False)

    df_blade_aero = pd.read_csv(dir_openFAST_data+"/"+aero_dyn_blade_file, skiprows=6, delim_whitespace=True,
                                names=["BlSpn", "BlCrvAC", "BlSwpAC", "BlCrvAng", "BlTwist", "BlChord", "BlAFID"])
    df_blade_aero.to_csv(dir_FAST+f"/blade_aero_dyn_original.{file_type}", index=False)

    use_lines = {"skiprows": None, "nrows": None}
    with open(dir_openFAST_data+"/"+elasto_dyn_blade_file, "r") as f:
        lines = f.readlines()
        for row_number, line in enumerate(lines):
            if "DISTRIBUTED BLADE PROPERTIES" in line:
                use_lines["skiprows"] = row_number+3
            if use_lines["skiprows"] is not None:
                if "----" in line:
                    use_lines["nrows"] = row_number-use_lines["skiprows"]
    df_structure = pd.read_csv(dir_openFAST_data+"/"+elasto_dyn_blade_file, delim_whitespace=True,
                               names=["BlFract","PitchAxis","StrcTwst", "BMassDen", "FlpStff", "EdgStff"], **use_lines)
    df_structure.to_csv(dir_FAST+f"/blade_elasto_dyn_original.{file_type}", index=False)


def scale_blade_by_R(dir_FAST: str, file_type: str, old_radius: float, new_radius: float) -> None:
    df_aero_original = pd.read_csv(dir_FAST+"/blade_aero_dyn_original."+file_type)
    df_elasto_original = pd.read_csv(dir_FAST+"/blade_elasto_dyn_original."+file_type)
    df_aero_R_scaled = copy(df_aero_original)
    df_elasto_R_scaled = copy(df_elasto_original)
    scaling_fac = new_radius/old_radius

    # scale aero file
    for parameter in ["BlSpn", "BlChord"]:
        df_aero_R_scaled[parameter] *= scaling_fac

    # scale elasto file
    for parameter, order in {"BMassDen": 2, "FlpStff": 4, "EdgStff": 4}.items():
        df_elasto_R_scaled[parameter] *= scaling_fac**order

    # save scaled files
    df_aero_R_scaled.to_csv(dir_FAST+"/blade_aero_dyn_R_scaled."+file_type, index=False)
    df_elasto_R_scaled.to_csv(dir_FAST+"/blade_elasto_dyn_R_scaled."+file_type, index=False)
    return

def incorporate_modifications(dir_FAST: str, file: str) -> None:
    df_aero_R_scaled = pd.read_csv(dir_FAST+"/blade_aero_dyn_R_scaled.dat")
    df_elasto_R_scaled = pd.read_csv(dir_FAST+"/blade_elasto_dyn_R_scaled.dat")
    df_aero_modified = copy(df_aero_R_scaled)
    df_elasto_modified = copy(df_elasto_R_scaled)
    df_modifications = pd.read_csv(dir_FAST+"/"+file)
    if not df_modifications["BlSpn"].equals(df_aero_R_scaled["BlSpn"]):
        raise ValueError("The radial positions of the R-scaled blade and the modifications must be the same, "
                         "but they are not.")
    # change blade properties
    df_aero_modified["BlChord"] = df_modifications["BlChord"]
    df_aero_modified["BlTwist"] = df_modifications["BlTwist"]
    df_elasto_modified["StrcTwst"] = df_modifications["BlTwist"]

    # scale blade properties from chord change
    chord_fac = df_modifications["BlChord"]/df_aero_R_scaled["BlChord"]
    for parameter, order in {"BMassDen": 2, "FlpStff": 4, "EdgStff": 4}.items():
        df_elasto_modified[parameter] *= chord_fac**order

    # save modified/scaled files
    df_aero_modified.to_csv(dir_FAST+"/blade_aero_dyn_modified.dat", index=False)
    df_elasto_modified.to_csv(dir_FAST+"/blade_elasto_dyn_modified.dat", index=False)


def compare_optimal_actual(dir_save: str,
                           file_optimum: str,
                           file_actual: str,
                           actual_radius: float,
                           aero_max_radius: float,
                           skip_rows: int=0,
                           skip_first_percentage: float=None):
    df_optimum= pd.read_csv(file_optimum)
    df_optimum = df_optimum[df_optimum["name"] == "IEA_10MW"]
    df_optimum.reset_index()
    df_actual = pd.read_csv(file_actual)
    fig, axes = plt.subplots(2)

    axes[0].plot(df_actual["BlSpn"]/actual_radius, df_actual["BlChord"]/actual_radius, label="real")
    axes[0].plot(df_optimum["radius"]/aero_max_radius, df_optimum["chord"]/aero_max_radius, label="aero max")

    axes[1].plot(df_actual["BlSpn"]/actual_radius, df_actual["BlTwist"], label="real twist")
    axes[1].plot(df_optimum["radius"]/aero_max_radius, df_optimum["twist"], label="aero max twist")

    positions, stall_margin = list(), list()
    df_optimum["r/R"] = df_optimum["radius"]/aero_max_radius
    for i, row in df_optimum.iterrows():
        if i < skip_rows and skip_first_percentage is None:
            continue
        if skip_first_percentage is not None:
            if row["r/R"] < skip_first_percentage/100:
                continue
        df = pd.read_csv(row["airfoil"])
        alpha_cl_max = df["alpha"][np.argmax(df["c_l"])]
        alpha_l2d_max = df["alpha"][np.argmax(df["c_l"]/df["c_d"])]
        stall_margin.append(alpha_cl_max-alpha_l2d_max)
        positions.append(row["r/R"])

    axes[1].plot(positions , stall_margin, label="aero max stall margin")
    helper.handle_axis(axes, title="Comparison of the scaled 10MW IEA to its aerodynamic optimised pendant",
                       x_label="r/R (-)", y_label=["chord/R (-)", "angle (°)"], grid=True, legend=True, line_width=4,
                       font_size=20)
    helper.handle_figure(fig, save_to=dir_save+"/actual_optimum_comparison.png", size=(18.5, 14))


def exists_already(df, **kwargs) -> bool:
    if df.query(order_to_query(kwargs)).empty:
        return False
    else:
        return True

def slice_results(df, **kwargs):
    for parameter, value in kwargs.items():
        df = df[df[parameter]==value]
    return df


def order_to_query(order: dict, negate_order: bool=False) -> str:
    compare_by = "==" if not negate_order else "!="
    query = str()
    for param, value in order.items():
        if type(value) in [str, list]:
            query += f"{param}{compare_by}'{value}' and "
        else:
            query += f"{param}{compare_by}{value} and "
    return query[:-5]


def calculate_root_moments(BEM_results_file: str, json_file: str, turbine_name: str) -> None:
    df_BEM = pd.read_csv(BEM_results_file)
    axial_loading = np.asarray([0]+df_BEM["f_n"]+[0])
    azimuthal_loading = np.asarray([0]+df_BEM["f_t"]+[0])
    radial_positions = np.asarray([df_BEM["root_radius"][0]]+df_BEM["r_centre"]+[df_BEM["rotor_radius"][0]])
    if os.path.isfile(json_file):
        with open(json_file) as f:
            moments = json.load(f)
            moments[turbine_name] = {"axial": simpson(axial_loading*radial_positions, radial_positions),
                                     "azimuthal": simpson(azimuthal_loading*radial_positions, radial_positions)}
    else:
        moments = {turbine_name: {"axial": simpson(axial_loading*radial_positions, radial_positions),
                                     "azimuthal": simpson(azimuthal_loading*radial_positions, radial_positions)}}
    with open(json_file, "w") as f:
        json.dump(moments, f)
    return