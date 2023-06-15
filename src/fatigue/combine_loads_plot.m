function [absolute_damage, damage, equivalent_loads] = combine_loads_plot(loads, t_loads, t_extrapolate, ...
    bin_width, eq_load_cycles, m, failure, safety_factor, threshold_cycles, max_loading_fac, ...
    U, threshold)
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
        [n_cycles, load_bin_edges] = histcounts(range, "BinWidth", bin_width);
        load_bin_centre = (load_bin_edges(1:end-1)+load_bin_edges(2:end))/2;
        load_bin_centre = load_bin_centre(n_cycles ~= 0);
        n_cycles = nonzeros(n_cycles);
        total_cycles = sum(n_cycles);

        ids_above_threshold = cumsum(n_cycles)/total_cycles>threshold_cycles;
        load_bin_centre_above_th = load_bin_centre(ids_above_threshold);
        n_cycles_above_th = n_cycles(ids_above_threshold);
        total_cycles_above_th = sum(n_cycles_above_th);
        load_centres_above_th = zeros(total_cycles_above_th, 1);
        where = 1;
        for i = 1:size(load_bin_centre_above_th, 2)
            for n = 1:n_cycles_above_th(i)
                load_centres_above_th(where) = load_bin_centre_above_th(i);
                where = where + 1;
            end
        end
        pd_loads_above_th = fitdist(load_centres_above_th, 'Weibull'); 
        above_th_range = load_bin_centre_above_th(end)-load_bin_centre_above_th(1);
        extrapolated_loads = load_bin_centre_above_th(1):bin_width:max_loading_fac*above_th_range+load_bin_centre_above_th(1);
        exceedence_prob = cdf(pd_loads_above_th, extrapolated_loads);
        % plot(extrapolated_loads, exceedence_prob)

        plot(extrapolated_loads, (1-(total_cycles-total_cycles_above_th)/total_cycles)*(1-exceedence_prob), "LineWidth", 2, "DisplayName", "Weibull fit")
        set(gca, "YScale", "log")
        hold on

        exceedance_probability = [];
        load_cycles = sum(n_cycles);
        for i = 1:size(load_bin_centre, 2)
            exceedance_probability(i) = sum(n_cycles(i:end)/load_cycles);
        end
        plot(load_bin_centre, exceedance_probability, "LineWidth", 2, "DisplayName", "Simulation below threshold")
        plot(load_bin_centre(ids_above_threshold), exceedance_probability(ids_above_threshold), "LineWidth", 2, "DisplayName", "Simulation above threshold")
        ax = gca;
        ax.Title.String = "U="+ U +"m/s, threshold=" + threshold;
        ax.XLabel.String = "Load range (Nm)";
        ax.YLabel.String = "Probability of exceeding the load range";
        ax.FontSize = 12;
        grid(ax, "on")
        legend(ax, "Location","best")
        saveas(gcf, "../../data/fatigue/plots/fit_U_"+U+"_th_"+threshold+".png")
        close(gcf)

        % above_threshold_probability = dictionary;
        % for el = 1:size(extrapolated_loads)
        %     above_threshold_probability(extrapolated_loads(i)) = cdf(pd_wbl_loads_abvove_th, ...
        %         extrapolated_loads(i));
        % end
        % above_threshold_probability



        % 
        % load_ranges = (edges(1:end-1)+edges(2:end))/2;
        % load_cycles = cycles_t_load/t_loads(i_load)*t_extrapolate(i_load);
        % absolute_damage(i_load) = load_cycles*load_ranges'.^m;
        % damage(i_load) = absolute_damage(i_load)/(failure.n*failure.load_range^m);
    end

    equivalent_loads = dictionary(eq_load_cycles, zeros(size(eq_load_cycles)));
    total_absolute_damage = sum(absolute_damage);
    for i_cycles = 1:max(size(eq_load_cycles))
        equivalent_loads(eq_load_cycles(i_cycles)) = (total_absolute_damage/eq_load_cycles(i_cycles))^(1/m);
    end
end

