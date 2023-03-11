clear
openFAST_data_file_type = "dat";
data_root = "../data/FAST_integration";
file_to_change = append(data_root, "/IEA_10MW.mat");
base = load(file_to_change);
file_additional_information = append(data_root, "/airfoil_additional_information.",openFAST_data_file_type);
scalar_info = readtable(file_additional_information);
dir_polars = append(data_root, "/polars");
dir_coordinates = append(data_root, "/coordinates");
n_airfoils = size(scalar_info.polar_file_name);


%% change Name
names = string(n_airfoils);
for i = 1:n_airfoils
    names(i) = string(i);
end
base.Airfoil.Name = cellstr(names).';

%% change Geometry
base.Airfoil.Geometry = cell(flip(n_airfoils));
for i = 1:n_airfoils
    file_name = scalar_info.coord_file_name{i};
    coords = readtable(append(dir_coordinates,"/",file_name), 'NumHeaderLines', 8);
    base.Airfoil.Geometry{i} = table2array(coords).';
end

%% change Alpha
base.Airfoil.Alpha = cell(flip(n_airfoils));
for i = 1:n_airfoils
    file_name = scalar_info.polar_file_name{i};
    data = readtable(append(dir_polars,"/",file_name));
    base.Airfoil.Alpha{i} = data.alpha;
end

%% change Cl
base.Airfoil.Cl = cell(flip(n_airfoils));
for i = 1:n_airfoils
    file_name = scalar_info.polar_file_name{i};
    data = readtable(append(dir_polars,"/",file_name));
    base.Airfoil.Cl{i} = data.c_l;
end

%% change Cd
base.Airfoil.Cd = cell(flip(n_airfoils));
for i = 1:n_airfoils
    file_name = scalar_info.polar_file_name{i};
    data = readtable(append(dir_polars,"/",file_name));
    base.Airfoil.Cd{i} = data.c_d;
end

%% change Cm
base.Airfoil.Cm = cell(flip(n_airfoils));
for i = 1:n_airfoils
    file_name = scalar_info.polar_file_name{i};
    data = readtable(append(dir_polars,"/",file_name));
    base.Airfoil.Cm{i} = data.c_m;
end

%% change SnSlope
base.Airfoil.CnSlope = scalar_info.CnSlope.';

%% change StallAngle1
base.Airfoil.StallAngle1 = scalar_info.StallAngle1.';

%% change StallAngle2
base.Airfoil.StallAngle2 = scalar_info.StallAngle2.';

%% change CritCn1
base.Airfoil.CritCn1= scalar_info.CritCn1.';

%% change CritCn2
base.Airfoil.CritCn2= scalar_info.CritCn2.';

%% SAVE CHANGES
save(file_to_change, "-struct", "base")















