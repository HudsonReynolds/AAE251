%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE251 Fall 2024
% Homework 9
% AAE251_HW9_Q1ab
% Author: Preston Wright and Hudson Reynolds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

alt = linspace(0,30,30/100+1);  % Altitude array [km]
powerAvailable = [];            % Initialized available power array [kW]
powerRequired = [];             % Initialized required power array [kN]

rhoSea = 1.2250;        % Density at sea level [kg/m^3]
g = 9.81;               % Gravitational acceleration [m/s^2]
m = 1315;               % Mass of the aircraft [kg]
mAD = 0.9;              % Air density exponent
powerMax = ;            % Maximum available power at sea level [kW]
propEff = 0.8;          % Propeller efficiency
area = 16.3;            % Wing area [m]
e = ;                   % Oswald effeciency factor
aspectRatio = ;         % Aspect ratio of wing
minimumDrag = ;         % Parasitic drag


%% Calculations

K = 1/(pi*e*aspectRatio);

for i = 1:length(alt)
    [~,~,~,rhoAlt] = atmosisa(alt(i)*1000);

    powerAvailable(i) = propEff*((rhoAlt/rhoSea)^mAD)*powerMax;
    powerRequired(i) = (4/3)*(((2*(m*g)^3)/(rho*area))*sqrt(3*(K^3)*minimumDrag));
end

%% Graphing

figure(1)
plot(alt,powerRequired)
hold on
plot(alt,powerAvailable)
grid minor
title("Available and Required Power Vs. Altitude")
xlabel("Altitude [km]")
ylabel("Power [kW]")
legend("Required Power", "Available Power")

