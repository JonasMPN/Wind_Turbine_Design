import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import numpy as np
import pandas as pd
import os, shutil
from src.helper_functions import Helper
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
    fig, axes = plt.subplots(4,1)
    for name in names:
        df_tmp = df[df["name"]==name]
        axes[0].plot(df_tmp["radius"], df_tmp["chord"], label=name)
        axes[1].plot(df_tmp["radius"], df_tmp["twist"])
        axes[2].plot(df_tmp["radius"], df_tmp["l2d"])
        axes[3].plot(df_tmp["radius"], df_tmp["l2d"]*df_tmp["chord"])
        power_indicator = np.trapz(df_tmp["l2d"]*df_tmp["chord"], df_tmp["radius"])
    helper.handle_axis([ax for ax in axes],
                       title=f"NREL 5MW, DTU 10MW, and IEA 10MW in comparison",
                       x_label="Radius in m",
                       y_label=["chord in m", "twist in °", r"$c_l/c_d$", r"load distribution"],
                       legend=True,
                       font_size=20,
                       line_width=4)
    helper.handle_figure(fig, file_figure=plot_dir+"/"+f"{names}"+".png", size=(18.5, 14))


def FAST_to_pandas(dir_FAST_data: str,
                   dir_save: str):
    file_blank = "IEA-10.0-198-RWT_AeroDyn15_Polar_"
    for i in range(30):
        idx = f"{i}" if i>9 else f"0{i}"
        file = dir_FAST_data+"/"+file_blank+idx+".dat"
        df = pd.read_csv(file, names=["alpha", "c_l", "c_d", "c_m"], skiprows=54, delim_whitespace=True)
        df.to_csv(dir_save+"/"+file_blank+f"{i}"+".dat", index=False)


def prepare_openFAST_to_FAST(dir_openFAST_data: str,
                             dir_coordinates: str,
                             dir_polars: str,
                             file_airfoil_names: str,
                             file_aero: str,
                             file_structure: str,
                             dir_FAST: str,
                             file_type: str="dat") -> None:
    """
    THIS FUNCTION OVERWRITES THE DIRECTORY "dir_FAST"!!!
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
    dir_coordinates = dir_openFAST_data+"/"+dir_coordinates
    dir_polars = dir_openFAST_data+"/"+dir_polars
    helper.create_dir(dir_FAST)
    helper.create_dir(dir_FAST+"/coordinates", overwrite=True)
    for file in os.listdir(dir_coordinates):
        new_file = file[::-1][file[::-1].find("."):][::-1]+file_type
        shutil.copyfile(dir_coordinates+"/"+file, dir_FAST+"/coordinates/"+new_file)

    helper.create_dir(dir_FAST+"/polars", overwrite=True)
    lines_to_use = {
        "StallAngle1": 18,
        "StallAngle2": 19,
        "CnSlope": 21,
        "CritCn1": 36,
        "CritCn2": 37}
    data = {col: list() for col in ["airfoil_name", "polar_file_name", "coord_file_name", *lines_to_use]}
    for file in os.listdir(dir_polars):
        new_file = file[::-1][file[::-1].find("."):][::-1]+file_type
        df = pd.read_csv(dir_polars+"/"+file, names=["alpha", "c_l", "c_d", "c_m"], skiprows=54, delim_whitespace=True)
        df.to_csv(dir_FAST+"/polars/"+new_file, index=False)
        with open(dir_polars+"/"+file, "r") as f:
            lines = f.readlines()
            for column, line_number in lines_to_use.items():
                value = str()
                for char in lines[line_number]:
                    if char != " ":
                        value += char
                    else:
                        break
                data[column].append(float(value))
    data["airfoil_name"] = pd.read_csv(dir_openFAST_data+"/"+file_airfoil_names)["name"]
    data["polar_file_name"] = os.listdir(dir_FAST+"/polars")
    data["coord_file_name"] = os.listdir(dir_FAST+"/coordinates")
    pd.DataFrame(data).to_csv(dir_FAST+f"/airfoil_additional_information.{file_type}", index=False)

    shutil.copyfile(dir_openFAST_data+f"/{file_aero}",dir_FAST+f"/blade_aero.{file_type}")
    shutil.copyfile(dir_openFAST_data+f"/{file_structure}",dir_FAST+f"/blade_structure.{file_type}")