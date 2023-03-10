import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import numpy as np
import pandas as pd
from helper_functions import Helper
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
                   dir_pandas_save: str):
    file_blank = "IEA-10.0-198-RWT_AeroDyn15_Polar_"
    for i in range(30):
        idx = f"{i}" if i>9 else f"0{i}"
        file = dir_FAST_data+"/"+file_blank+idx+".dat"
        df = pd.read_csv(file, names=["alpha", "c_l", "c_d", "c_m"], skiprows=54, delim_whitespace=True)
        df.to_csv(dir_pandas_save+"/"+file_blank+f"{i}"+".dat", index=False)
