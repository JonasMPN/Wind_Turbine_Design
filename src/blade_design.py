import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from helper_functions import Helper
from data_handling import AirfoilInterpolator
helper = Helper()


class BladeApproximation:
    def __init__(self,
                 root_dir: str,
                 blade_dir: str,
                 blade_filename: str,
                 save_dir: str,
                 blade_name: str):
        self.root = root_dir
        self.dir_blade = root_dir+"/"+blade_dir
        self.file_blade = self.dir_blade+"/"+blade_filename
        self.dir_save = root_dir+"/"+save_dir
        self.blade_name = blade_name
        self.file_save = self.dir_save+"/optimum_results.dat"
        helper.create_dir(root_dir+"/"+save_dir)
        try:
            self.df_results = pd.read_csv(self.file_save)
        except FileNotFoundError:
            self.df_results = pd.DataFrame(columns=["name", "radius", "chord", "twist", "l2d", "airfoil"])

        self.df_blade = pd.read_csv(self.file_blade)
        self.tsr = None
        self.R = None
        self.B = None
        self.interp = None

        self.column_positions = None
        self.column_airfoil_path = None
        self.column_rel_thickness = None

        self.column_cl = None
        self.column_cd = None
        self.column_alpha = None

    def set_rotor(self,
                  tip_speed_ratio: float,
                  rotor_radius: float,
                  number_of_blades: int) -> None:
        self.tsr = tip_speed_ratio
        self.R = rotor_radius
        self.B = number_of_blades

    def set_blade_columns(self,
                          column_positions: str,
                          column_airfoil_path: str,
                          column_rel_thickness: str="rel_thickness",
                          interpolation_required: bool=False,
                          **pandas_read_csv):
        self.column_positions = column_positions
        self.column_airfoil_path = column_airfoil_path
        self.column_rel_thickness = column_rel_thickness
        if interpolation_required:
            # Interpolation only works on relative thickness and angle of attack (alpha). Both the blade file and the
            # airfoil files need to have the same column name for the relative thickness.
            self.interp = AirfoilInterpolator(dir_data=self.dir_blade,
                                              data_profiles=self.df_blade[column_airfoil_path].unique(),
                                              **pandas_read_csv)

    def set_airfoil_columns(self,
                            column_cl: str="c_l",
                            column_cd: str="c_d",
                            column_alpha: str="alpha"):
        self.column_cl = column_cl
        self.column_cd = column_cd
        self.column_alpha = column_alpha

    def chord_and_twist(self,
                        actual_radius: float,
                        skip_rows: int=0,
                        skip_first_percentage: float=None) -> None:
        """
        save_to_file overrides any data that has the same blade name
        :param skip_rows:
        :param skip_first_percentage:
        :param save_to_file:
        :param plot:
        :return:
        """
        self.df_blade["r/R"] = self.df_blade[self.column_positions]/actual_radius
        if self.interp is not None:
            df_profiles = self.interp.get_df_profiles()
            alphas = df_profiles[self.column_alpha].unique()
            rel_thicknesses = df_profiles[self.column_rel_thickness].unique()
            rel_thicknesses = np.sort(np.asarray(rel_thicknesses))[::-1]
        positions, twist, chord, l2ds, all_c_l, all_c_d = list(), list(), list(), list(), list(), list()
        airfoils_used = list()
        inflow_angles = list()
        for i, row in self.df_blade.iterrows():
            if i < skip_rows and skip_first_percentage is None:
                continue
            if skip_first_percentage is not None:
                if row["r/R"] < skip_first_percentage/100:
                    continue

            if self.interp is None:
                airfoil_path = row[self.column_airfoil_path]
                airfoil_path = airfoil_path if type(airfoil_path) is str else str(int(airfoil_path))
                airfoils_used.append(self.dir_blade+"/"+airfoil_path)
                df_airfoil = pd.read_csv(self.dir_blade+"/"+airfoil_path)
                lift2drag = df_airfoil[self.column_cl]/df_airfoil[self.column_cd]
                c_l = df_airfoil[self.column_cl].iloc[np.argmax(lift2drag)]
                c_d = df_airfoil[self.column_cd].iloc[np.argmax(lift2drag)]
                alpha = df_airfoil[self.column_alpha].iloc[np.argmax(lift2drag)]
            else:
                rel_thickness = row[self.column_rel_thickness]
                airfoils_used.append("interpolated")
                for inner, outer in zip(rel_thicknesses[:-1], rel_thicknesses[1:]):
                    if rel_thickness == inner or rel_thickness == outer:
                        c_ls = df_profiles[df_profiles[self.column_rel_thickness] == rel_thickness][self.column_cl]
                        lift2drag = c_ls/df_profiles[df_profiles[self.column_rel_thickness] == rel_thickness][self.column_cd]
                        break
                    if inner > rel_thickness and rel_thickness > outer:
                        interp_args = [self.column_alpha, self.column_rel_thickness]
                        constraints = {self.column_cl: {self.column_rel_thickness: [inner, outer]},
                                       self.column_cd: {self.column_rel_thickness: [inner, outer]}}
                        self.interp.prepare_interpolation({self.column_cl: interp_args,
                                                           self.column_cd: interp_args}, constraints=constraints)
                        c_ls, c_ds, lift2drag = list(), list(), list()
                        for alpha in alphas:
                            c_l = self.interp[self.column_cl](alpha, rel_thickness)
                            c_d = self.interp[self.column_cd](alpha, rel_thickness)
                            c_ls.append(c_l)
                            c_ds.append(c_d)
                            lift2drag.append(c_l/c_d)
                        break
                c_l = c_ls[np.argmax(lift2drag)]
                c_d = c_ds[np.argmax(lift2drag)]
                alpha = alphas[np.argmax(lift2drag)]

            positions.append(row["r/R"]*self.R)
            l2ds.append(max(lift2drag))
            fac = self.tsr*row["r/R"]*(1+2/(9*self.tsr**2*row["r/R"]**2))
            chord.append(16*np.pi*self.R/(9*self.tsr*c_l*self.B*np.sqrt(4/9+fac**2)))
            inflow_angle = np.arctan(2/(3*fac))*180/np.pi
            inflow_angles.append(inflow_angle)
            twist.append(inflow_angle-alpha)
            all_c_l.append(c_l)
            all_c_d.append(c_d)

        if self.blade_name in self.df_results["name"].unique():
            self.df_results = self.df_results[self.df_results["name"] != self.blade_name]
        self.df_results = pd.concat([self.df_results,
                                     pd.DataFrame({"name": self.blade_name,
                                                   "radius": positions,
                                                   "chord": chord,
                                                   "twist": twist,
                                                   "l2d": l2ds,
                                                   "airfoil": airfoils_used})])
        self.df_results.to_csv(self.file_save, index=False)
        return


