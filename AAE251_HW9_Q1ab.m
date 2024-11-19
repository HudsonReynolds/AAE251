%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE251 Fall 2024
% Homework 9
% AAE251_HW9_Q1ab
% Author: Preston Wright and Hudson Reynolds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

alt = linspace(0,30,30*100+1);  % Altitude array [km]
powerAvailable = [];            % Initialized available power array [kW]
powerRequired = [];             % Initialized required power array [kN]

rhoSea = 1.2250;        % Density at sea level [kg/m^3]
g = 9.81;               % Gravitational acceleration [m/s^2]
m = 1315;               % Mass of the aircraft [kg]
mAD = 0.6;              % Air density exponent
K = 0.054;              % Span efficiency
powerMax = 216;         % Maximum available power at sea level [kW]
propEff = 0.8;          % Propeller efficiency
area = 16.3;            % Wing area [m]
parasiteDrag = 0.026;   % Parasitic drag


%% Calculations

% Loop through the altitudes calculating the available and required power
for i = 1:length(alt)
    [~,~,~,rhoAlt] = atmosisa(alt(i)*1000,extended="on");

    powerAvailable(i) = propEff*((rhoAlt/rhoSea)^mAD)*powerMax;
    powerRequired(i) = (4/3)*sqrt(((2*(m*g)^3)/(rhoAlt*area))*sqrt(3*(K^3)*parasiteDrag))/1000;
end

%% Graphing

% Output the required and available power with respect to altitude
figure(1)
plot(alt,powerRequired)
hold on
plot(alt,powerAvailable)
grid minor
title("Available and Required Power Vs. Altitude")
xlabel("Altitude [km]")
ylabel("Power [kW]")
legend("Required Power", "Available Power")

