clear
load_range_FA = 1;
load_range_SS = 0.5;
n_points = 12;

total_eq_load_range = zeros(n_points,1);
FA_eq_load_range = zeros(n_points, 1);
SS_eq_load_range = zeros(n_points, 1);
angles = zeros(n_points, 1);
for n = 0:n_points-1
    angle = 2*pi*n/n_points;
    angles(n+1) = rad2deg(angle);
    FA_eq_load_range(n+1) = load_range_FA*abs(cos(angle));
    SS_eq_load_range(n+1) = load_range_SS*abs(sin(angle));
    total_eq_load_range(n+1) = sqrt(FA_eq_load_range(n+1)^2+SS_eq_load_range(n+1)^2);
end

