clear all
close
clc

run('P2_main.m');

% %% Plot 10m
% 
% datapoints_sum = length(data_time_1)+length(data_time_2)+length(data_time_3)+length(data_time_4); % Total number of data points
% 
% wind_x_10 = linspace(0,23,24);
% pd_x_10 = linspace(0,23,100);
% 
% pd_wbl_wind=fitdist(wind_10','Weibull'); 
% wind_y = pdf(pd_wbl_wind,pd_x_10);
% 
% histogram_wind_10 = histogram(wind_10,'BinEdges',wind_x_10,'Normalization','probability');
% hold on
% 
% plot(pd_x_10,wind_y,'r','LineWidth',2);
% xlabel('Wind speed (m/s)');
% ylabel('Probability density');
% title('Weibull Distribution of Mean Wind Speed (10m)');
% legend('Wind Data Histogram','Weibull Distribution');
% xlim([0 30]);
% hold off
%% Plot hub height (110 m)

figure();
wind_x_hub = linspace(0.5,wind_mean_hub_max+0.5,wind_mean_hub_max+1);
pd_x_hub = linspace(0,wind_mean_hub_max,100);

pd_wbl_wind_hub=fitdist(wind_hub','Weibull'); 
wind_y_hub = pdf(pd_wbl_wind_hub,pd_x_hub);

histogram_wind_hub = histogram(wind_hub,'BinEdges',wind_x_hub,'Normalization','probability');
save("../../data/fatigue//task_2/pdf_wind_speed_histogram.mat", "histogram_wind_hub")
hold on

plot(pd_x_hub,wind_y_hub,'r','LineWidth',2);
xlabel('Wind speed (m/s)');
ylabel('Probability density');
title('Weibull Distribution of Wind Speed (100m)');
legend('Wind Data Histogram','Weibull Distribution');
xlim([0 30]);
hold off

figure()
histfit(wind_hub,wind_mean_hub_max+1,'Weibull');