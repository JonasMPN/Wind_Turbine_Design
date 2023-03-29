function FAST_object = airfoil(FAST_object, scalar_info, dir_coordinates, dir_polars)
    %% check whether an element is at r=0m
    if length(FAST_object.Blade.Radius) ~= length(scalar_info.polar_file_name)
        scalar_info(1,:) = [];
    end
    n_airfoils = size(scalar_info.polar_file_name);
    
    %% change Name
    names = string(n_airfoils);
    for i = 1:n_airfoils
        names(i) = string(i);
    end
    FAST_object.Airfoil.Name = cellstr(names).';
    
    %% change Geometry
    FAST_object.Airfoil.Geometry = cell(flip(n_airfoils));
    for i = 1:n_airfoils
        file_name = scalar_info.coord_file_name{i};
        coords = readtable(append(dir_coordinates,"/",file_name), 'NumHeaderLines', 8);
        FAST_object.Airfoil.Geometry{i} = table2array(coords).';
    end
    
    %% change Alpha
    FAST_object.Airfoil.Alpha = cell(flip(n_airfoils));
    for i = 1:n_airfoils
        file_name = scalar_info.polar_file_name{i};
        data = readtable(append(dir_polars,"/",file_name));
        FAST_object.Airfoil.Alpha{i} = data.alpha;
    end
    
    %% change Cl
    FAST_object.Airfoil.Cl = cell(flip(n_airfoils));
    for i = 1:n_airfoils
        file_name = scalar_info.polar_file_name{i};
        data = readtable(append(dir_polars,"/",file_name));
        FAST_object.Airfoil.Cl{i} = data.c_l;
    end
    
    %% change Cd
    FAST_object.Airfoil.Cd = cell(flip(n_airfoils));
    for i = 1:n_airfoils
        file_name = scalar_info.polar_file_name{i};
        data = readtable(append(dir_polars,"/",file_name));
        FAST_object.Airfoil.Cd{i} = data.c_d;
    end
    
    %% change Cm
    FAST_object.Airfoil.Cm = cell(flip(n_airfoils));
    for i = 1:n_airfoils
        file_name = scalar_info.polar_file_name{i};
        data = readtable(append(dir_polars,"/",file_name));
        FAST_object.Airfoil.Cm{i} = data.c_m;
    end
    
    %% change SnSlope
    FAST_object.Airfoil.CnSlope = scalar_info.CnSlope.';
    
    %% change StallAngle1
    FAST_object.Airfoil.StallAngle1 = scalar_info.StallAngle1.';
    
    %% change StallAngle2
    FAST_object.Airfoil.StallAngle2 = scalar_info.StallAngle2.';
    
    %% change CritCn1
    FAST_object.Airfoil.CritCn1= scalar_info.CritCn1.';
    
    %% change CritCn2
    FAST_object.Airfoil.CritCn2= scalar_info.CritCn2.';
end

