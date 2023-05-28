clear
data = readtable("../data/results/fatigue_wind.dat");
N_eq = 5e6;
m = 4;
%%
damage = data.n.*data.M_range.^m;
total_damage = sum(damage);
rel_damage = damage/total_damage;
eq_M_range = (sum(damage)/N_eq)^(1/m);


