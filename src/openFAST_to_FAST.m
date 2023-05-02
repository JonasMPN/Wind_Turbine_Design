clear
addpath('openFAST_to_FAST_functions\');

%% user input
data_root = "../data/FAST_integration";
openFAST_data_file_type = "dat";
base_file = append(data_root, "/IEA_7MW.mat");
new_file = append(data_root, "/IEA_7MW.mat");
new_radius = 90;
old_radius = 99;
scale_fac = new_radius/old_radius;

% Nacelle
shaft_tilt = 6; % degrees
hub_radius = 2.4*90/99.155; % metres

m_yam_bearing = 93457*scale_fac^3;
m_nacelle_turret_and_nose = 109450*scale_fac^3;
m_generator = 151651;
m_converter = 5250;
m_shaft = 78894*scale_fac*(3910381.2837784197/5125004.865491613)^(2/3);

nacelle_mass = m_yam_bearing+m_nacelle_turret_and_nose+m_generator+m_converter+m_shaft;
hub_mass = 81707*scale_fac^3;


% Drivetrain
gen_efficiency = 0.9236922238; 
gearbox_ratio = 1;

% Blade
rotor_precone = 4; % degrees

% Controls
torque_max = 8000000;
torque_rated = 7578283;
optimal_mode_gain = 4.7012e+6;
omega_A = 3.3677;
omega_B = 3.5;
omega_B2 = 9.3;
omega_C = 9.5493;

% Tower
eff_density = 8500;
hub_height = 108;


%% the following lines assume a certain directory structure
FAST_object = load(base_file);
aero_data = readtable(append(data_root, "/blade_aero_dyn_modified.", openFAST_data_file_type));
structure_data = readtable(append(data_root, "/blade_elasto_dyn_modified.", openFAST_data_file_type));
tower_data = readtable(append(data_root, "/tower.", openFAST_data_file_type));
file_additional_information = append(data_root, "/airfoil_additional_information.",openFAST_data_file_type);
scalar_info = readtable(file_additional_information);
dir_polars = append(data_root, "/polars");
dir_coordinates = append(data_root, "/coordinates");

%% change hub and shaft parameters. This needs to happen before the blade changes!
disp("Changing hub and shaft")
FAST_object = nacelle(FAST_object, shaft_tilt, hub_radius, hub_mass, nacelle_mass);

%% change blade properties (this needs to be run before the airfoil changes!)
disp("Changing blade properties")
FAST_object = blade(FAST_object, aero_data, structure_data, rotor_precone);

%% change airfoil properties
disp("Changing airfoil properties")
FAST_object = airfoil(FAST_object, scalar_info, dir_coordinates, dir_polars);

%% change drivetrain
disp("Changing drivetrain")
FAST_object = drivetrain(FAST_object, gen_efficiency, gearbox_ratio);

%% change controls
disp("Changing controls")
FAST_object = controls(FAST_object, torque_max, torque_rated, optimal_mode_gain, omega_A, omega_B, omega_B2, omega_C);

%% change tower 
disp("Changing tower")
FAST_object = tower(FAST_object, tower_data, hub_height, eff_density);

%% save changes
save(new_file, "-struct", "FAST_object")
disp("File created")


