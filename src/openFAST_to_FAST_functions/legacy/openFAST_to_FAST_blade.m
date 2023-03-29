% This file needs to be run before openFAST_to_FAST_airfoil
clear

data_root = "../data/FAST_integration";
openFAST_data_file_type = "dat";
file_to_change = append(data_root, "/IEA_7MW.mat");
base = load(file_to_change);
aero_data = readtable(append(data_root, "/blade_aero_dyn_modified.", openFAST_data_file_type));
structure_data = readtable(append(data_root, "/blade_elasto_dyn_modified.", openFAST_data_file_type));

%% check whether an element is at r=0m
if aero_data.BlSpn(1) == 0
    aero_data(1,:) = [];
    structure_data(1,:) = [];
end
n_positions = size(aero_data.BlSpn);

%% change Radius
base.Blade.Radius = aero_data.BlSpn;

%% change Chord
base.Blade.Chord = aero_data.BlChord;

%% change Twist and necessary pitch control values
fine_pitch = aero_data.BlTwist(end);
base.Blade.Twist = aero_data.BlTwist-fine_pitch;
base.Control.Pitch.Fine = fine_pitch;
base.Control.Pitch.Min = -2+fine_pitch;
print("The twist distribution has been reset to be zero at the tip." + ...
    "The fine pitch was adjusted accordingly.")

%% change NFoil
base.Blade.NFoil = double(1:n_positions)';

%% change IFoil
base.Blade.IFoil = double(1:n_positions)';

%% change Mass density
base.Blade.Mass = structure_data.BMassDen;

%% change flap stiffness
base.Blade.EIflap = structure_data.FlpStff;

%% change edge stiffness
base.Blade.EIedge = structure_data.EdgStff;

%% change pitch axis
base.Blade.PitchAxis = structure_data.PitchAxis;

%% change thickness
base.Blade.Thickness = ones(n_positions(1),1);

%% SAVE CHANGES
save(file_to_change, "-struct", "base")