function [absolute_damage, damage, equivalent_loads] = combine_loads(loads, t_loads, t_extrapolate, ...
    bin_width, eq_load_cycles, m, failure, safety_factor)
    % Calculates the equivalent load of the combined load of single or mutliple load cases.
 
    % loads: matrix of load time series. Each row is a different load case

    % t_loads: vector of the time span of the load cases. Each row
    %   corresponds to the row of 'loads'. If given as scalar, this value is
    %   used for all load cases

    % t_extrapolate: vector of the time spans for which each load is
    %   extrapolated. Each row corresonds to the row of 'loads'

    % bin_width: bin width of the histogram of the load ranges

    % eq_load_cycles: array of the number of load cycles for the equivalent
    % loads

    % m: Wöhler exponent

    % failure: a structure with field 'n' (scalar) and 'load_range' (scalar) for which the
    % material fails.
    
    % safety factor: factor that is multiplied to the load ranges
    if max(size(t_loads)) == 1
        t_loads = t_loads*ones(size(loads,1));
    end

    if size(loads, 1) ~= size(t_loads, 1) || size(loads, 1) ~= size(t_extrapolate, 1)
        error("'loads', 't_loads', and 't_extrapolate' are not in accord regarding the number of load cases.")
    end
    
    absolute_damage = zeros(size(loads,1), 1);
    damage = zeros(size(loads,1), 1);
    for i_load = 1:size(loads, 1)
        [range, mean] = rainflow(loads(i_load, :));
        range = safety_factor*range; % FASTtool returns kNm
        [cycles_t_load, edges] = histcounts(range, "BinWidth", bin_width);
        load_ranges = (edges(1:end-1)+edges(2:end))/2;
        load_cycles = cycles_t_load/t_loads(i_load)*t_extrapolate(i_load);
        absolute_damage(i_load) = load_cycles*load_ranges'.^m;
        damage(i_load) = absolute_damage(i_load)/(failure.n*failure.load_range^m);
    end

    equivalent_loads = dictionary(eq_load_cycles, zeros(size(eq_load_cycles)));
    total_absolute_damage = sum(absolute_damage);
    for i_cycles = 1:max(size(eq_load_cycles))
        equivalent_loads(eq_load_cycles(i_cycles)) = (total_absolute_damage/eq_load_cycles(i_cycles))^(1/m);
    end
end

