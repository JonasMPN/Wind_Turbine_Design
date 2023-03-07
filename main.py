import numpy as np
import pandas as pd
import scipy.interpolate as interpolate


class AirfoilInterpolator:
    def __init__(self,
                 dir_data: str,
                 data_profiles: list or str,
                 **pandas_read_csv):
        if type(data_profiles) == str:
            data_profiles = [data_profiles]
        self.dir_data = dir_data
        self.df_profiles = list()
        for file in data_profiles:
            self.df_profiles.append(pd.read_csv(dir_data + "/" + file, **pandas_read_csv))

    def interpolate(self,
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
        for to_inter, arguments in to_interpolate.items():
            if type(arguments) == str:
                arguments = [arguments]
            points = {var: list() for var in arguments}
            values = list()
            df_data_points = self.df_profiles

            if to_inter in constraints.keys():
                df_data_points = pd.DataFrame(columns=df_data_points.columns)
                for arg_parameter, interval in constraints[to_inter].items():
                    if arg_parameter not in arguments:
                        raise ValueError(f"Interval {arg_parameter}={interval} cannot be set for the "
                                         f"interpolation of {to_inter}, because {to_inter} is only based on "
                                         f"{arguments}.")
                    df_data_points = pd.concat(df_data_points,
                                               df_data_points[df_data_points[arg_parameter==interval[0]] &
                                               df_data_points[arg_parameter==interval[1]]], ignore_index=True)

            for df in df_data_points:
                for _, row in df.iterrows():
                    for arg in arguments:
                        points[arg].append(row[arg])
                    values.append(row[to_inter])
            points = np.asarray([all_points for all_points in points.values()]).T
            self.__dict__[to_inter] = interpolate.LinearNDInterpolator(points, values)
