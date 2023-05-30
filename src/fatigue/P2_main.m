clear all
close
clc

dir_data = "../../data/fatigue/task_2/";

ncdisp(append(dir_data, "/",'21_data.nc'));
data_latitude_1  = ncread(append(dir_data, "/",'21_data.nc'),'latitude');
data_longitude_1  = ncread(append(dir_data, "/", '21_data.nc'),'longitude');
data_time_1 = ncread(append(dir_data, "/", '21_data.nc'),'time');
data_u_1 = ncread(append(dir_data, "/", '21_data.nc'),'u100'); % wind speed (longitudeXlatitudeXtime)
data_v_1 = ncread(append(dir_data, "/", '21_data.nc'),'v100'); 
wind_1 = zeros(1,length(data_time_1));

data_latitude_2  = ncread(append(dir_data, "/", '121314_data.nc'),'latitude');
data_longitude_2  = ncread(append(dir_data, "/", '121314_data.nc'),'longitude');
data_time_2 = ncread(append(dir_data, "/", '121314_data.nc'),'time');
data_u_2 = ncread(append(dir_data, "/", '121314_data.nc'),'u100'); % wind speed (longitudeXlatitudeXtime)
data_v_2 = ncread(append(dir_data, "/", '121314_data.nc'),'v100'); 
wind_2 = zeros(1,length(data_time_2));
    
data_latitude_3  = ncread(append(dir_data, "/", '151617_data.nc'),'latitude');
data_longitude_3  = ncread(append(dir_data, "/", '151617_data.nc'),'longitude');
data_time_3 = ncread(append(dir_data, "/", '151617_data.nc'),'time');
data_u_3 = ncread(append(dir_data, "/", '151617_data.nc'),'u100'); % wind speed (longitudeXlatitudeXtime)
data_v_3 = ncread(append(dir_data, "/", '151617_data.nc'),'v100'); 
wind_3 = zeros(1,length(data_time_3));

data_latitude_4  = ncread(append(dir_data, "/", '181920_data.nc'),'latitude');
data_longitude_4  = ncread(append(dir_data, "/", '181920_data.nc'),'longitude');
data_time_4 = ncread(append(dir_data, "/", '181920_data.nc'),'time');
data_u_4 = ncread(append(dir_data, "/", '181920_data.nc'),'u100'); % wind speed (longitudeXlatitudeXtime)
data_v_4 = ncread(append(dir_data, "/", '181920_data.nc'),'v100'); 
wind_4 = zeros(1,length(data_time_4));
%% Compute wind speed magnitude 
        for k = 1:length(data_time_1)
              wind_1(k) = sqrt(data_u_1(k)^2+data_v_1(k)^2);
        end
        for k = 1:length(data_time_2)
              wind_2(k) = sqrt(data_u_2(k)^2+data_v_2(k)^2);
        end
        for k = 1:length(data_time_3)
              wind_3(k) = sqrt(data_u_3(k)^2+data_v_3(k)^2);
        end
        for k = 1:length(data_time_4)
              wind_4(k) = sqrt(data_u_4(k)^2+data_v_4(k)^2);
        end
%% Computing the mean wind speed
% wind_1 = wind_1+2;
% wind_2 = wind_2+2;
% wind_3 = wind_3+2;
% h1 = 60;
% alpha_sheer = 0.11;
% wind_10 = [wind_1,wind_2,wind_3,wind_4];
% wind_10_mean = mean([wind_1,wind_2,wind_3,wind_4]);
% wind_60 = [wind_1,wind_2,wind_3,wind_4]*(h1/10)^alpha_sheer;
% h2 = 100;
% z_0 = 0.0002;
% wind_hub = wind_60*(log(h2/z_0)/log(h1/z_0));
wind_hub = [wind_1,wind_2,wind_3,wind_4];
wind_hub_mean = mean(wind_hub);
% wind_hub_mean = mean([wind_1,wind_2,wind_3,wind_4]*(z_hub/10)^alpha_sheer);
wind_mean_hub_max = ceil(max(wind_hub));

% param_wbl_wind_10 = wblfit(wind_10);
param_wbl_wind_hub = wblfit(wind_hub);