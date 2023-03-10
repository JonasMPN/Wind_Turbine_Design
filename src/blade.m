clear
data_root = "../data/FAST_integration";
file_to_change = append(data_root, "/IEA_10MW.mat");
base = load(file_to_change);
aero_data = readtable(append(data_root, "/blade_aero.txt"));
structure_data = readtable(append(data_root, "/blade_structure.txt"));
n_positions = size(aero_data.BlSpn);

%% change Radius
base.Blade.Radius = aero_data.BlSpn;

%% change Chord
base.Blade.Chord = aero_data.BlChord;

%% change Twist
base.Blade.Twist = aero_data.BlTwist;

%% change NFoil
file_additional_information = append(data_root, "/airfoil_additional_information.dat");
airfoil_names = readtable(file_additional_information).airfoil_name;
id_airfoils = dictionary(airfoil_names.', 1:length(airfoil_names));
base.Blade.NFoil = double(n_positions).';
for i = 1:n_positions(1)
    base.Blade.NFoil(i) = id_airfoils(aero_data.BlAFID(i));
end

%% change IFoil
base.Blade.IFoil = double(n_positions).';
for i = 1:n_positions(1)
    base.Blade.IFoil(i) = i;
end
%% change Mass density
base.Blade.Mass = structure_data.BMassDen;

%% change flap stiffness
base.Blade.EIflap = structure_data.FlpStff;

%% change edge stiffness
base.Blade.EIedge = structure_data.EdgStff;

%% change pitch axis
base.Blade.PitchAxis = structure_data.PitchAxis;

%% change thickness
base.Blade.Thickness = ones(30,1);

%% SAVE CHANGES
save(file_to_change, "-struct", "base")