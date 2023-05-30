clear
close
clc

global para
run('P2_main.m');

para.rho = 1.225;

para.R1 = 141/2;
para.V1_rated = 11.5;
para.P1_rated = 7e6;
para.A1 = pi * para.R1^2;
para.Cp1= para.P1_rated/(0.5*para.rho*para.V1_rated^3*para.A1);


V = linspace(0,25,1000);
V_AEP = 0.5:1:24.5;

P1_temp = 0.5*para.Cp1 * para.rho * V.^3 * para.A1;
P1 = zeros(1,length(V)); % Initialize

P1_AEP_temp = 0.5*para.Cp1 * para.rho * V_AEP.^3 * para.A1;
P1_AEP = zeros(1,length(V_AEP)); % Initialize

for i = 1:length(V)
   if (P1_temp(i)<para.P1_rated)
       P1(i) = 0.5*para.Cp1 * para.rho * V(i)^3 * para.A1;
   else
       P1(i) = para.P1_rated; 
   end
end


for i = 1:length(V_AEP)
   if (P1_AEP_temp(i)<para.P1_rated)
       P1_AEP(i) = 0.5*para.Cp1 * para.rho * V_AEP(i)^3 * para.A1;
   else
       P1_AEP(i) = para.P1_rated; 
   end
end




AEP_1_inst = zeros(1,length(V_AEP));
for i = 1:length(V_AEP) % wind_mean_hub_max
    AEP_1_inst(i) = (exp(-((i-1)/param_wbl_wind_hub(1))^param_wbl_wind_hub(2))-exp(-((i)/param_wbl_wind_hub(1))^param_wbl_wind_hub(2)))*P1_AEP(i);
end


AEP_1 = 8766*3600*sum(AEP_1_inst);

AEP_1_MWh = AEP_1/3e6;
%% Plot


plot(V, P1);
xlabel('Wind Speed (m/s)');
ylabel('Generated Power (W)');
grid on
grid minor
legend('7 MW Wind Turbine');
title('Power Curve');
ylim([0 8e6]);
hold off