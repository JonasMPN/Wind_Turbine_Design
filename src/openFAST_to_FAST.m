clear

data_root = "../data/FAST_integration";
openFAST_data_file_type = "dat";
file_to_change = append(data_root, "/IEA_7MW.mat");
base = load(file_to_change);

%% change gearbox ratio
base.Drivetrain.Gearbox.Ratio = 1;
