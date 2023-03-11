clear
radius = 85; % this is the radius to which all data from the the directory FAST_integration will be scaled
openFAST_data_file_type = "dat";
data_root = "../data/FAST_integration";
file_to_change = append(data_root, "/IEA_10MW.mat");
base = load(file_to_change);
aero_data = readtable(append(data_root, "/blade_aero_dyn.",openFAST_data_file_type));
structure_data = readtable(append(data_root, "/blade_elasto_dyn.",openFAST_data_file_type));
n_positions = size(aero_data.BlSpn);
scaling_fac = radius/max(aero_data.BlSpn);

%% change Radius
base.Blade.Radius = aero_data.BlSpn.*scaling_fac;

%% change Chord
base.Blade.Chord = aero_data.BlChord.*scaling_fac;

%% change Twist
base.Blade.Twist = aero_data.BlTwist;

%% change NFoil
base.Blade.NFoil = 1:n_positions;

%% change IFoil
base.Blade.IFoil = 1:n_positions;

%% change Mass density
base.Blade.Mass = structure_data.BMassDen.*scaling_fac^2;

%% change flap stiffness
base.Blade.EIflap = structure_data.FlpStff.*scaling_fac;

%% change edge stiffness
base.Blade.EIedge = structure_data.EdgStff.*scaling_fac;

%% change pitch axis
base.Blade.PitchAxis = structure_data.PitchAxis;

%% change thickness
base.Blade.Thickness = ones(30,1);

%% SAVE CHANGES
save(file_to_change, "-struct", "base")