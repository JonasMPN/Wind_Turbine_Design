clear
addpath('openFAST_to_FAST_functions\');

%% user input
data_root = "../data/FAST_integration";
openFAST_data_file_type = "dat";
base_file = append(data_root, "/controlled7MW4.mat");
new_file = append(data_root, "/controlled7MW4_v2.mat");
new_radius = 90;
ref_radius = 99;
new_rated_power = 7e6;
max_tip_speed = 90;
initial_C_P = 0.49;
density = 1.225;
tip_speed_ratio = 10.58;
scale_fac = new_radius/ref_radius;

% Nacelle
shaft_tilt = 6; % degrees
hub_radius = 2.4*90/99.155; % metres

m_ref_yam_bearing = 93457;
m_ref_nacelle_turret_and_nose = 109450;
m_ref_shaft = 78894;
m_ref_hub = 81707;

m_generator = 164918.8144;
m_converter = 5250;

% Drivetrain
gen_efficiency = 0.978; 
gearbox_ratio = 1;
gen_inertia = 686056.1061;

% Blade
rotor_precone = 4; % degrees

% Controls
cut_in_wind_speed = 4;
control_transition_margin = 0.05;
limit_to_demanded_torque = 4.7403/4.3094;
ref_demanded_torque = 11.7e6;
min_omega = 6;

% Tower
eff_density = 8500;
hub_height = 108;

% Losses
cable_efficiency = 0.9991;
converter_efficiency = 0.98;

% Certification
IEA_wind_class = 1;
IEA_turbulence_class = 2; % A=1, B=2, C=3

%% the following lines assume a certain directory structure
FAST_object = load(base_file);
aero_data = readtable(append(data_root, "/blade_aero_dyn_modified.", openFAST_data_file_type));
structure_data = readtable(append(data_root, "/blade_elasto_dyn_modified.", openFAST_data_file_type));
tower_data = readtable(append(data_root, "/tower.", openFAST_data_file_type));
file_additional_information = append(data_root, "/airfoil_additional_information.",openFAST_data_file_type);
scalar_info = readtable(file_additional_information);
dir_polars = append(data_root, "/polars");
dir_coordinates = append(data_root, "/coordinates");


%% creating data tables
% drivetrain
elec_efficiency = gen_efficiency*converter_efficiency*converter_efficiency;
drivetrain_data = table(elec_efficiency, gearbox_ratio, gen_inertia);

% controls
opt_mode_gain = pi*density*new_radius^5*initial_C_P/(2*tip_speed_ratio^3);
omega_C = max_tip_speed/new_radius*60/(2*pi);
omega_B2 = omega_C*(1-control_transition_margin);
omega_A = min_omega;
omega_B = omega_A*(1+control_transition_margin);
torque_demanded = new_rated_power/(omega_C*2*pi/60*elec_efficiency);
torque_max = limit_to_demanded_torque*torque_demanded;
torque_min = 2.5e4;
control_data = table(torque_demanded, omega_A, omega_B, omega_B2, omega_C, opt_mode_gain, ...
    torque_max, torque_min);

% nacelle
m_yam_bearing = m_ref_yam_bearing*scale_fac^3;
m_nacelle_turret_and_nose = m_ref_nacelle_turret_and_nose*scale_fac^3;
m_shaft = m_ref_shaft*scale_fac*(torque_demanded/ref_demanded_torque)^(2/3);
nacelle_mass = m_yam_bearing+m_nacelle_turret_and_nose+m_generator+m_converter+m_shaft;
hub_mass = m_ref_hub*scale_fac^3;
nacelle_data = table(shaft_tilt, hub_radius, hub_mass, nacelle_mass);

% certification
certification_data = table(new_radius, IEA_wind_class, IEA_turbulence_class);

%% change hub and shaft parameters. This needs to happen before the blade changes!
disp("Changing hub and shaft")
FAST_object = nacelle(FAST_object, nacelle_data);

%% change blade properties (this needs to be run before the airfoil changes!)
disp("Changing blade properties")
FAST_object = blade(FAST_object, aero_data, structure_data, rotor_precone);

%% change airfoil properties
disp("Changing airfoil properties")
FAST_object = airfoil(FAST_object, scalar_info, dir_coordinates, dir_polars);

%% change drivetrain
disp("Changing drivetrain")
FAST_object = drivetrain(FAST_object, drivetrain_data);

%% change controls
disp("Changing controls")
FAST_object = controls(FAST_object, control_data);

%% change tower 
disp("Changing tower")
FAST_object = tower(FAST_object, tower_data, hub_height, eff_density);

%% change certification
disp("Changing certification")
FAST_object = certification(FAST_object, certification_data);

%% save changes
save(new_file, "-struct", "FAST_object")
disp("File created")


