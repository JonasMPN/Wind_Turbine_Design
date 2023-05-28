clear
close all
f1p = 1/4;
M1p = 2.2;
f1e = 0.9;
M1e = 0.8;
t_duration = 600;
N_eq = 5e6;
m = 9;

t = 0:0.01:t_duration;
load_p1 = load(M1p, f1p, t);
load_p1e = load_p1 + load(M1e, f1e, t);
t_1e = 1630*3600;
t_1 = 20*365*24*3600*0.906-t_1e;

loads = [load_p1; load_p1e];
t_loads = t_duration;
t_extrapolate = [t_1; t_1e];

eq_load_range = equivalent_load_range(loads, t_loads, t_extrapolate, 0.01, N_eq, m);


function load = load(load_range, frequency, times)
    load = load_range/2*sin(times*2*pi*frequency);
end