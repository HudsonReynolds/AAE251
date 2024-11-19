%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE251 Fall 2024
% Homework 9
% AAE251_HW9_Q1cd
% Author: Preston Wright and Hudson Reynolds
% Description: Sets up and calculates the available and required thrust with
% respect to altitude for a given jet aircraft, plotting those values versus
% the altitude used to calculate them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

alt = linspace(0,30,30*100+1);   % Altitude array [km]
thrustAvailable = [];            % Initialized available power array [kW]
thrustRequired = [];             % Initialized required power array [kN]

rhoSea = 1.2250;        % Density at sea level [kg/m^3]
g = 9.81;               % Gravitational acceleration [m/s^2]
m = 33100;              % Mass of the aircraft [kg]
mAD = 0.6;              % Air density exponent
thrustMax = 55620;      % Maximum available thrust at sea level [kN]
K = 0.05;               % Wingspan efficiency
parasiteDrag = 0.015;   % Parasitic drag


%% Calculations

% Loop through the altitudes calculating the available and required thrust
for i = 1:length(alt)
    [~,~,~,rhoAlt] = atmosisa(alt(i)*1000, extended="on");
    thrustAvailable(i) = ((rhoAlt/rhoSea)^mAD)*thrustMax;
    thrustRequired(i) = 2*m*g*sqrt(K*parasiteDrag);
end

%% Graphing

% Output the required and available thrust with respect to altitude
figure(1)
plot(alt,thrustRequired)
hold on
plot(alt,thrustAvailable)
grid minor
title("Available and Required Thrust Vs. Altitude")
xlabel("Altitude [km]")
ylabel("Thrust [kN]")
legend("Required Thrust", "Available Thrust", location="northeast")

