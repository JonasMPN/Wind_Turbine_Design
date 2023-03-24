import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as interpolate
from helper_functions import Helper
import scipy
helper=Helper()


class BEM:
    def __init__(self,
                 data_root: str):
        self.rotor_radius = None
        self.root_radius = None
        self.n_blades = None
        self.air_density = None
        self.root = data_root
        self.a_prime = 0
        try:
            self.df_results = pd.read_csv(data_root+"/BEM_results.dat")
        except FileNotFoundError:
            self.df_results = pd.DataFrame()

        self.interp = {"c_l": None, "c_d": None}
        self.implemented_glauert_correction = ["none", "dtu", "tud"]
        self.implemented_tip_correction = ["none", "dtu", "tud"]
        self.implemented_root_correction = ["none", "tud"]
        self.blade_end_correction = 1

    def set_constants(self,
                      rotor_radius: float,
                      root_radius: float,
                      n_blades: int,
                      air_density: float) -> None:
        self._set(**{param: value for param, value in locals().items() if param != "self"})
        return None

    def solve_TUD(self,
                  blade_file: str,
                  wind_speed: float,
                  tip_speed_ratio: float,
                  pitch: float or np.ndarray,
                  start_radius: float = None,
                  max_convergence_error: float=1e-6,
                  max_iterations: int=500,
                  tip_loss_correction: bool=True,
                  root_loss_correction: bool=True) -> None:
        """
        Function to run the BEM loop
        Glauert_correction: either 'tud' (TU Delft) or 'dtu' (Denmark's TU). Same for blade_end_correction
        All angles must be in rad.
        :param wind_speed:
        :param tip_speed_ratio:
        :param pitch: IN DEGREE
        :param resolution:
        :return:
        """
        start_radius = start_radius if start_radius is not None else self.root_radius
        # Initialise the result containers
        results = {
            "r_centre": list(),     # radius used for the calculations
            "a": list(),            # Axial Induction factor
            "a_prime": list(),      # Tangential induction factor
            "f_n": list(),          # Forces normal to the rotor plane in N/m
            "f_t": list(),          # Forces tangential in the rotor plane in N/m
            "bec": list(),          # blade end correction (depending on 'tip' and 'root')
            "C_T": list(),          # thrust coefficient
            "alpha": list(),        # angle of attack
            "alpha_max": list(),    # maximum angle of attack
            "alpha_best": list(),   # angle of attack for maximum L/D
            "circulation": list(),  # magnitude of the circulation using Kutta-Joukowski
            "v0": list(),           # flow velocity normal to rotor plane
            "tsr": list(),          # tip speed ratio
            "pitch": list(),         # pitch in degree
        }
        # delete data with same wind speed, tip speed ratio and pitch angle.
        try:
            self.df_results = self.df_results.loc[~((self.df_results["tsr"]==tip_speed_ratio) &
                                                    (self.df_results["v0"]==wind_speed) &
                                                    (self.df_results["pitch"]==pitch))]
        except KeyError:
            pass
        pitch = np.deg2rad(pitch)
        # Calculate the rotational speed
        omega = tip_speed_ratio*wind_speed/self.rotor_radius
        # Loop along the span of the blade (blade sections)
        df = pd.read_csv(blade_file)
        print(f"Doing BEM for v0={wind_speed}, tsr={tip_speed_ratio}, pitch={pitch}")
        for idx, row in df.iterrows():      # Take the left and right radius of every element
            if row["BlSpn"] < start_radius:
                continue
            # insert airfoil polars into the interpolator
            df_tmp = pd.read_csv("../data/IEA_10MW/"+row["BlAFID"])
            self.interp = {"c_l": interpolate.interp1d(df_tmp["alpha"], df_tmp["c_l"]),
                           "c_d": interpolate.interp1d(df_tmp["alpha"], df_tmp["c_d"])}
            a, a_new, a_prime, converged = 0.8, 0.1, 0, False
            for i in range(max_iterations):
                # get inflow angle and speed for the airfoil
                phi, inflow_speed = self._flow(a=a, a_prime=a_prime, wind_speed=wind_speed, rotational_speed=omega,
                                               radius=row["BlSpn"])
                # get combined lift and drag coefficient projected into the normal and tangential direction
                _, _, _, c_n, c_t = self._phi_to_aero_values(phi=phi, twist=np.deg2rad(row["BlTwist"]), pitch=pitch,
                                                             tip_seed_ratio=tip_speed_ratio, university="tud")
                # get thrust force (in N/m) of the whole turbine at the current radius
                thrust = self._aero_force(inflow_speed, row["BlChord"], c_n)*self.n_blades
                # calculate thrust coefficient that results from the blade element
                C_T = thrust/(self.air_density*wind_speed**2*np.pi*row["BlSpn"])
                # get Glauert corrected axial induction factor
                a_new = self._a(C_T=C_T)
                # get the combined (tip and root) correction factor
                blade_end_correction = self._blade_end_correction(which="tud", tip=tip_loss_correction,
                                                                  root=root_loss_correction, radius=row["BlSpn"],
                                                                  tip_seed_ratio=tip_speed_ratio, a=a_new)
                # correct the Glauert corrected axial induction factor with the blade end losses
                a_new /= blade_end_correction
                # update the axial induction factor for the next iteration
                a = 0.75*a+0.25*a_new
                # get the tangential force (in N/m) of the whole turbine at the current radius
                f_tangential = self._aero_force(inflow_speed, row["BlChord"], c_t)*self.n_blades
                # get the tangential induction factor that corresponds to this force AND correct it for tip losses
                a_prime = self._a_prime(f_tangential, row["BlSpn"], wind_speed, a, tip_speed_ratio)/blade_end_correction
                # check if the axial induction factor has converged. If it has, the tangential induction factor has too
                if np.abs(a-a_new) < max_convergence_error:
                    converged = True
                    break
            # notify user if loop did not converge, but was stopped by the maximum number of iterations
            if not converged:
                print(f"BEM did not converge for the blade element {row['BlSpn']}m. Current change after"
                      f" {max_iterations}: {np.abs(a-a_new)}.")
            # Now that we have the converged axial induction factor, we can get the rest of the values
            phi, inflow_speed = self._flow(a=a, a_prime=a_prime, wind_speed=wind_speed, rotational_speed=omega,
                                           radius=row["BlSpn"])
            alpha, c_l, _, c_n, c_t = self._phi_to_aero_values(phi=phi, twist=np.deg2rad(row["BlTwist"]),
                                                               pitch=pitch,  radius=row["BlSpn"],
                                                               tip_seed_ratio=tip_speed_ratio, university="tud")
            # Assemble the result output structure
            results["r_centre"].append(row["BlSpn"])
            results["a"].append(a)
            results["a_prime"].append(a_prime)
            results["f_n"].append(self._aero_force(inflow_velocity=inflow_speed, chord=row["BlChord"], force_coefficient=c_n))
            results["f_t"].append(self._aero_force(inflow_velocity=inflow_speed, chord=row["BlChord"], force_coefficient=c_t))
            results["bec"].append(blade_end_correction)
            results["C_T"].append(self._C_T(a))
            results["alpha"].append(alpha)
            results["alpha_max"].append(df_tmp["alpha"].loc[df_tmp["c_l"].idxmax()])
            results["alpha_best"].append(df_tmp["alpha"].loc[(df_tmp["c_l"]/df_tmp["c_d"]).idxmax()])
            results["circulation"].append(1/2*inflow_speed*c_l*row["BlChord"])
            results["v0"].append(wind_speed)
            results["tsr"].append(tip_speed_ratio)
            results["pitch"].append(pitch)
        self.df_results = pd.concat([self.df_results, pd.DataFrame(results)])
        self.df_results.to_csv(self.root+"/BEM_results.dat", index=False)
        return None

    def _calculate_thrust(self, f_n, radial_positions):
        """
            Calculate thrust from the normal forces. Account for f_t = 0 at the tip.
        f_n: normal forces
        radial_positions: radial position along the blade matching the positions of f_n
        n_blades:   number of blades
        radius:     max radius
        """
        thrust = self.n_blades*scipy.integrate.simpson([*f_n, 0], [*radial_positions, self.rotor_radius])
        return thrust

    def _calculate_power(self, f_t, radial_positions, omega):
        """
            Calculate power from the normal forces. Account for f_n = 0 at the tip.
        f_t: tangential forces
        radial_positions: radial position along the blade matching the positions of f_n
        n_blades:   number of blades
        radius:     max radius
        omega:      [rad/s] rotational speed
        """
        power = omega*self.n_blades*scipy.integrate.simpson(np.multiply([*radial_positions, self.rotor_radius], [*f_t, 0]),
                                                           [*radial_positions, self.rotor_radius])
        return power

    def _calc_ct(self, thrust, velocity):
        """
            Calculate the thrust coefficient ct
        """
        return thrust/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**2))

    def _calc_cp(self, power, velocity):
        """
            Calculate the power coefficient ct
        """
        return power/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**3))

    def _calc_ct_distribution(self, f_n, velocity):
        """
        Calculate the distribution of the thrust coefficient along the blade
        f_n: normal forces along the blade
        radius: maximum radius of the Blade
        velocity: fluid velocity od V0
        """
        f_n = np.array(f_n)
        return f_n/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**3))

    def _calc_cp_distribution(self, f_t, velocity):
        """
        Calc the distribution of the power coeff. along the blade
        f_t: tangential forces along the blade
        radius: maximum radius of the Blade
        velocity: fluid velocity od V0
        """
        f_t = np.array(f_t)
        return f_t/(0.5*np.pi*(self.rotor_radius**2)*self.air_density*(velocity**3))

    def _phi_to_aero_values(self, phi: float, twist:float or np.ndarray, pitch: float,tip_seed_ratio: float,
                            university: str,
                            radius: float=None, a: float=None, blade_end_correction_type: str=None, tip: bool=None,
                            root: bool=None) -> tuple:
        """
        phi, twist and pitch need to be in rad. Alpha is returned as degree
        :param phi:
        :param twist:
        :param pitch:
        :param tip_seed_ratio:
        :param university:
        :param radius:
        :param a:
        :param blade_end_correction_type:
        :param tip:
        :param root:
        :return:
        """
        alpha = np.rad2deg(phi-(twist+pitch))
        c_l = self.interp["c_l"](alpha)
        c_d = self.interp["c_d"](alpha)
        c_n = self._c_normal(phi, c_l, c_d)
        c_t = self._c_tangent(phi, c_l, c_d)
        if university == "tud":
            return alpha, c_l, c_d, c_n, c_t
        elif university == "dtu":
            return alpha, c_l, c_d, c_n, c_t, self._blade_end_correction(which=blade_end_correction_type, tip=tip,
                                                                         root=root, phi=phi, radius=radius, a=a,
                                                                         tip_seed_ratio=tip_seed_ratio)

    def _set(self, **kwargs) -> None:
        """
        Sets parameters of the instance. Raises an error if a parameter is trying to be set that doesn't exist.
        :param kwargs:
        :return:
        """
        existing_parameters = [*self.__dict__]
        for parameter, value in kwargs.items():
            if parameter not in existing_parameters:
                raise ValueError(f"Parameter {parameter} cannot be set. Settable parameters are {existing_parameters}.")
            self.__dict__[parameter] = value
        return None

    def _assert_values(self):
        not_set = list()
        for variable, value in vars(self).items():
            if value is None:
                not_set.append(variable)
        if len(not_set) != 0:
            raise ValueError(f"Variable(s) {not_set} not set. Set all variables before use.")

    def _update_a_prime(self, local_solidity: float, c_tangential: float, blade_end_correction: float,
                        inflow_angle: float) -> None:
        self.a_prime = local_solidity*c_tangential*(1+self.a_prime)/(4*blade_end_correction*np.sin(inflow_angle)*
                                                                     np.cos(inflow_angle))

    def _aero_force(self, inflow_velocity: float, chord: float, force_coefficient: float):
        """
        Calculates the tangential force per unit span.
        :param inflow_velocity:
        :param chord:
        :param c_normal:
        :return:
        """
        return 1/2*self.air_density*inflow_velocity**2*chord*force_coefficient

    def _equate_blade_element_and_momentum(self, glauert_correction: str, a: float, blade_end_correction: float,
                                           phi: float, local_solidity: float, c_normal: float):
        """
        Function to calculate the difference of the CT values (Blade element vs momentum theory) from a given axial induction
        """
        if glauert_correction == "none":
            return 1/((4*blade_end_correction*np.sin(phi)**2)/(local_solidity*c_normal)+1)-a
        elif glauert_correction == "dtu":
            if a < 1/3:
                return 1/((4*blade_end_correction*np.sin(phi)**2)/(local_solidity*c_normal)+1)-a
            else:
                return local_solidity*((1-a)/np.sin(phi))**2*c_normal-4*a*blade_end_correction*(1-a/4*(5-3*a))
        elif glauert_correction == "tud":
            CT_1 = 1.816
            if a < 1-np.sqrt(CT_1)/2:
                return 4*a*(1-a)
            else:
                return local_solidity*((1-a)/np.sin(phi))**2*c_normal-(CT_1-4*(np.sqrt(CT_1)-1)*(1-a))

    def _a_prime(self, F_tangential: float, radius: float, wind_speed: float, a: float, tip_speed_ratio: float):
        return F_tangential/(4*self.air_density*np.pi*radius**2/self.rotor_radius*wind_speed**2*(1-a)*tip_speed_ratio)

    def _blade_end_correction(self, which: str, tip: bool, root: bool, radius: float,
                              tip_seed_ratio: float, a: float, phi: float=None) -> float:
        """
        Different Prandtl correction methods.
        :param which:
        :param tip:
        :param root:
        :param phi:
        :param radius:
        :param tip_seed_ratio:
        :param a:
        :return:
        """
        F = 1
        if tip:
            if which == "dtu":
                if np.sin(np.abs(phi)) < 0.001:
                    pass
                else:
                    F = 2/np.pi*np.arccos(np.exp(-(self.n_blades*(self.rotor_radius-radius))/(2*radius*np.sin(np.abs(phi)))))
            elif which == "tud":
                d = 2*np.pi/self.n_blades*(1-a)/(np.sqrt(tip_seed_ratio**2+(1-a)**2))
                F = 2/np.pi*np.arccos(np.exp(-np.pi*((self.rotor_radius-radius)/self.rotor_radius)/d))
            else:
                raise ValueError(f"Parameter 'blade_end_correction' must be one of "
                                 f"{self.implemented_tip_correction}, but was {which}.")
        if root:
            if which == "tud":
                d = 2*np.pi/self.n_blades*(1-a)/(np.sqrt(tip_seed_ratio**2+(1-a)**2))
                F *= 2/np.pi*np.arccos(np.exp(-np.pi*((radius-self.root_radius)/self.rotor_radius)/d))
            else:
                raise ValueError(f"Parameter 'blade_end_correction' must be one of "
                                 f"{self.implemented_root_correction}, but was {which}.")
        return F

    @staticmethod
    def _a(C_T: float) -> float:
        C_T1 = 1.816
        CT_2 = 2*np.sqrt(C_T1)-C_T1
        if C_T < CT_2:
            return 1/2-np.sqrt(1-C_T)/2
        else:
            return 1+(C_T-C_T1)/(4*np.sqrt(C_T1)-4)

    @staticmethod
    def _C_T(a: float):
        C_T1 = 1.816
        if a < 1-np.sqrt(C_T1)/2:
            return 4*a*(1-a)
        else:
            return C_T1-4*(np.sqrt(C_T1)-1)*(1-a)


    @staticmethod
    def _c_normal(phi: float, c_lift: float, c_drag: float) -> float:
        """
        Calculates an aerodynamic "lift" coefficient according to a coordinate transformation with phi
        :param phi: angle between flow and rotational direction in rad
        :param c_lift: lift coefficient old coordinate system
        :param c_drag: lift coefficient old coordinate system
        :return: Normal force in Newton
        """
        return c_lift*np.cos(phi)+c_drag*np.sin(phi)

    @staticmethod
    def _c_tangent(phi: float, c_lift: float, c_drag: float) -> float:
        """
        Calculates an aerodynamic "drag" coefficient according to a coordinate transformation with phi
        :param phi: angle between flow and rotational direction in rad
        :param c_lift: lift coefficient old coordinate system
        :param c_drag: lift coefficient old coordinate system
        :return: Normal force in Newton
        """
        return c_lift*np.sin(phi)-c_drag*np.cos(phi)

    @staticmethod
    def _local_solidity(chord: float, radius: float, n_blades: int) -> float:
        """
        Calculates the local solidity
        :param chord: in m
        :param radius: distance from rotor axis to blade element in m
        :param n_blades: number of blades
        :return: local solidity
        """
        return n_blades*chord/(2*np.pi*radius)

    @staticmethod
    def _flow(a: float, a_prime: float, wind_speed: float, rotational_speed: float, radius: float) -> [float, float]:
        """
        Function to calculate the inflow angle and inflow velocity based on the two induction factors, the inflow
        velocity, radius and angular_velocity.
        :param a:
        :param a_prime:
        :return:
        """
        v_axial = wind_speed*(1-a)
        v_azimuthal = rotational_speed*radius*(1+a_prime)
        return np.tan(v_axial/v_azimuthal), np.sqrt(v_axial**2+v_azimuthal**2)

    @staticmethod
    def _get_twist(r, r_max):
        """
            function to get the twist along the blade in radians
            r: radial position along the blade in [m]
            r_max: maximum radius of the blade in [m]
            out: twist in radians
        """
        return np.radians(14*(1-r/r_max))

    @staticmethod
    def _get_chord(r, r_max):
        """
            function to calculate chord along the span in m
            r: radial position along the blade in [m]
            r_max: maximum radius of the blade in [m]
        """
        return 3*(1-r/r_max)+1

    @staticmethod
    def _tangential_induction_factor(phi: float, local_solidity: float, c_tangent: float, tip_loss_correction: float)\
            -> float:
        return 1/((4*tip_loss_correction*np.sin(phi)*np.cos(phi))/(local_solidity*c_tangent)-1)
