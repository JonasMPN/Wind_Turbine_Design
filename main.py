import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

class AirfoilInterpolator:
    def __init__(self,
                 dir_data: str,
                 data_profiles: list or str,
                 **pandas_read_csv):
        if type(data_profiles) == str:
            data_profiles = [data_profiles]
        self.dir_data = dir_data
        self.df = pd.DataFrame()
        for file in data_profiles:
            self.df = pd.concat([self.df, pd.read_csv(dir_data+"/"+file, **pandas_read_csv)],
                                         ignore_index=True)

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
                self.__dict__[to_inter] = interpolate.interp1d(points.flatten(), values)
            else:
                self.__dict__[to_inter] = interpolate.LinearNDInterpolator(points, values)


dir_data = "data"
file = "combined_data_new_v2.txt"
tsr = 8
interp = AirfoilInterpolator(dir_data=dir_data, data_profiles=file, **{"sep":","})
df_blade_data = pd.read_csv("data/blade_data_new.txt")
df_blade_data["r/R"] = df_blade_data["radius"]/df_blade_data["radius"].max()

positions, chord, twist = list(), list(), list()
# line below this assumes airfoil data is listed from the smallest rel_thickness to largest
df_blades_data =pd.read_csv("data/combined_data_new_v2.txt")
airfoil_rel_thicknesses = df_blades_data["rel_thickness"].unique()[::-1]
alphas = df_blades_data[df_blades_data["rel_thickness"] == 241]["alpha"].unique()
for _, row in df_blade_data.iterrows():
    positions.append(row["r/R"])
    rel_thickness = row["rel_thickness"]
    interp_args, constraints = None, None
    for inner, outer in zip(airfoil_rel_thicknesses[:-1], airfoil_rel_thicknesses[1:]):
        if rel_thickness == inner or rel_thickness == outer:
            interp_args, constraints = ["alpha"], {"c_l": {"rel_thickness": [rel_thickness]},
                                                   "c_d": {"rel_thickness": [rel_thickness]}}
        elif rel_thickness < inner and rel_thickness > outer:
            interp_args, constraints = ["alpha", "rel_thickness"], {"c_l": {"rel_thickness": [inner, outer]},
                                                                    "c_d": {"rel_thickness": [inner, outer]}}
        if interp_args is not None:
            interp.interpolate({"c_l": interp_args, "c_d": interp_args}, constraints=constraints)
            c_ls, lift2drag = list(), list()
            for alpha in alphas:
                if len(interp_args) == 1:
                    c_ls.append(interp.c_l(alpha))
                    lift2drag.append(interp.c_l(alpha)/interp.c_d(alpha))
                else:
                    c_ls.append(interp.c_l(alpha, rel_thickness))
                    lift2drag.append(interp.c_l(alpha, rel_thickness)/interp.c_d(alpha, rel_thickness))

            c_l = c_ls[np.argmax(lift2drag)]
            fac = tsr*row["r/R"]*(1+2/(9*(tsr*row["r/R"])**2))
            chordtoradius = 16*np.pi/(9*tsr*c_l*3)
            twist.append(np.arctan(2/(3*tsr*row["r/R"]*(1+2/(9*(tsr*row["r/R"])**2))))*180/np.pi)
            break


plt.plot(positions, twist)
plt.show()
# interp.interpolate(to_interpolate={"c_l": ["alpha", "rel_thickness"]}, constraints={"c_l": {"rel_thickness": [0,2]}})

