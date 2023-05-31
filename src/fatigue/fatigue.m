clear
close all
%% User input
dir_root = "../../data/fatigue";
pdf_wind_speed = load(append(dir_root, "/task_2/pdf_wind_speed_histogram.mat")).histogram_wind_hub;
dir_load_cases = append(dir_root, "/FASTtool_load_cases");
file_base_simulation = append(dir_load_cases, "/simulation");

seeds = [1 2];
wind_speed_range = 1:27;
blades = [1 2 3];
load_types = ["RootMEdg" "RootMFlp"]; 
t_loads = 600; % seconds
life_time = 20; % years
m = 9;
load_bin_width = 0.01;
equivalent_load_cycles = [5e6 1e7 20*365*24*3600];
material_max_stress = 600e6; % in Pa=Nm^-2
fos_load = 1;
fos_material = 1.7;
fos_severity = 1.3;
root_chord = 4.5; % m
second_MoA_root = 118585511936.00403/21.79e9; % stiffness/E in m^4

%% pre structure
loads = strings(1, size(load_types, 2)*size(blades, 2));
n_blades = size(blades, 2);
n_load_types = size(load_types, 2);
for i_blade = 0:n_blades-1
    for i_load_type = 1:n_load_types
        loads(i_blade*n_load_types+i_load_type) = append(load_types(i_load_type), ...
            num2str(i_blade+1));
    end
end
loads_str = regexprep(mat2str(load_types), '"', '');

%% Get load results
concatenated_file = append(dir_root, "/loads_", loads_str, "_blades_", mat2str(blades), ...
    "_seeds_", mat2str(seeds), "_U_range_", num2str(wind_speed_range(1)), "_", ...
    num2str(wind_speed_range(end)), ".mat");

if isfile(concatenated_file)
    concatenated_data = load(concatenated_file);
    disp("These loads have been concatenated before. Loading saved file.")
else
    % build data structure
    disp("Concatenating and saving results.")
    file_example = append(file_base_simulation, "_U=", num2str(wind_speed_range(1)), ...
        ".00_seed=", num2str(seeds(1)), ".mat");
    example_results = load(file_example);
    load_matrix = zeros(size(wind_speed_range, 2), ...
        size(seeds, 2)*size(example_results.Time, 1));
    t_extrapolate = zeros(size(wind_speed_range, 2), 1);
    base_struct = struct("load_matrix", load_matrix, "t_extrapolate", t_extrapolate);
    concatenated_data = struct;
    for i = 1:size(loads, 2)
        concatenated_data.(loads(i)) = base_struct;
    end
    bin_centres = (pdf_wind_speed.BinEdges(1:end-1)+pdf_wind_speed.BinEdges(2:end))/2;
    wind_speed_prob = dictionary(bin_centres, pdf_wind_speed.Values);
    
    % concatenate results
    for l = 1:size(loads, 2)
        disp(append("Working on ", loads(l)))
        load_type = loads(l);
        for i = 1:size(wind_speed_range, 2)
            load_time_series = [];
            wind_speed = wind_speed_range(i);
            for n = 1:size(seeds)+1
                seed = seeds(n);
                simulation_file = append(file_base_simulation, "_U=", num2str(wind_speed), ...
                    ".00_seed=", num2str(seed), ".mat");
                sim_results = load(simulation_file);
                load_time_series = [load_time_series sim_results.(load_type)'*1e3];
            end
            concatenated_data.(load_type).load_matrix(i, :) = load_time_series;
            load_time_for_lifetime = life_time*365*24*3600*wind_speed_prob(wind_speed);
            concatenated_data.(load_type).t_extrapolate(i) = load_time_for_lifetime;
        end
    end
    disp("Saving concatenated data.")
    save(concatenated_file, "-struct", "concatenated_data")
end

%% Calculation of damage and damage equivalent loads
disp("Calculating the absolute damage, the relative damage until failure, and " + ...
    "equivalent moment ranges.")
combined_loads_file = append(dir_root, "/combined_loads_", loads_str, "_blades_", ...
    mat2str(blades), "_seeds_", mat2str(seeds), "_U_range_", ...
    num2str(wind_speed_range(1)), "_", num2str(wind_speed_range(end)), ".mat");
if isfile(combined_loads_file)
    combined_loads = load(combined_loads_file);
    disp("These loads have been combined before. Loading saved file")
else
    r_root = root_chord/2;
    load_range = 2*material_max_stress*second_MoA_root/r_root;
    failure = struct("n", 1, "load_range", load_range);
    safety_factor = fos_severity*fos_material*fos_load;
    combined_loads = struct;
    t_combined_loads = t_loads*max(size(seeds));
    for i = 1:size(loads, 2)
        disp(append("Working on load ", loads(i)))
        load_type = loads(i);
        plot_data = concatenated_data.(load_type);
        
        [absolute_damage, damage, equivalent_load_ranges] = combine_loads(plot_data.load_matrix, ...
            t_combined_loads, plot_data.t_extrapolate, load_bin_width, equivalent_load_cycles, ...
            m, failure, safety_factor);
        combined_loads.(loads(i)) = struct("absolute_damage", absolute_damage, "damage", ...
            damage, "equivalent_load_ranges", equivalent_load_ranges);
    end
    disp("Saving combined loads")
    save(combined_loads_file, "-struct", "combined_loads")
end

%% Exporting and/or plotting data. A basis for plotting here in MATLAB is provided.
% axes = cell(size(load_types, 2)+1, 1);
plot_data_types = cat(2, "equivalent_load_range", load_types);
plot_data = cell(size(plot_data_types,  2), 1);
for i_ax = 1:size(plot_data_types, 2)
    % figure;
    % axes{i_ax} = gca;
    % hold(axes{i_ax});
    if plot_data_types(i_ax) == "equivalent_load_range"
        N_eq = combined_loads.(loads(1)).equivalent_load_ranges.keys;
        cells = cell(size(N_eq, 1), size(loads, 2) + 1);
        for N = 1:size(N_eq, 1)
            cells{N, 1} = N_eq(N);
        end
        header = cat(2, "equivalent_cycles", loads);
    else
        cells = cell(size(wind_speed_range, 2), size(blades, 2) + 1);
        for i = 1:size(wind_speed_range, 2)
            cells{i, 1} = wind_speed_range(i);
        end
        header = cat(2, "wind_speeds", "blade"+blades);
    end
    plot_data{i_ax} = cell2table(cells, "VariableNames", header);
end

i_ax = repmat(2:size(load_types, 2)+1, 1, size(blades, 2));
i_blade = repmat(1:size(blades, 2), 1, size(load_types, 2));
for i_load_type = 1:size(loads, 2)
    blade_no = i_blade(i_load_type);
    load_type = loads(i_load_type);
    load_data = combined_loads.(load_type);
    plot_data{1}.(load_type) = load_data.equivalent_load_ranges.values;
    plot_data{i_ax(i_load_type)}.("blade"+num2str(blade_no)) = load_data.absolute_damage;

    % loglog(axes{1}, N_eq, load_data.equivalent_load_ranges.values, "-x", "DisplayName", ...
    %     load_type);
    % plot(axes{i_ax(i_load_type)}, wind_speed_range, load_data.absolute_damage, "-x", "DisplayName", load_type);
end
% hold(axes{1}, "off")
% set(axes{1}, "YScale", "log", "XScale", "log")

% for i_ax = 1:size(load_types, 2)+1
%     legend(axes{i_ax});
%     grid(axes{i_ax});
% end

for i = 1:size(plot_data_types, 2)
    file_name = append(dir_root, "/plot_data/", plot_data_types(i), ".dat");
    writetable(plot_data{i, 1}, file_name);
end





