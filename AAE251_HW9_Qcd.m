%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE251 Fall 2024
% Homework 9
% AAE251_HW9_Q1cd
% Author: Preston Wright and Hudson Reynolds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

alt = linspace(0,30,30/100+1);  % Altitude array [km]
thrustAvailable = [];            % Initialized available power array [kW]
thrustRequired = [];             % Initialized required power array [kN]

rhoSea = 1.2250;        % Density at sea level [kg/m^3]
g = 9.81;               % Gravitational acceleration [m/s^2]
m = ;                   % Mass of the aircraft [kg]
mAD = ;                 % Air density exponent
thrustMax = ;           % Maximum available thrust at sea level [kN]
area = ;                % Wing area [m]
e = ;                   % Oswald effeciency factor
aspectRatio = ;         % Aspect ratio of wing
minimumDrag = ;         % Parasitic drag


%% Calculations

K = 1/(pi*e*aspectRatio);

for i = 1:length(alt)
    [~,~,~,rhoAlt] = atmosisa(alt(i)*1000);

    thrustAvailable(i) = ((rhoAlt/rhoSea)^mAD)*thrustMax;
    thrustRequired(i) = 2*m*g*sqrt(K*minimumDrag);
end

%% Graphing

figure(1)
plot(alt,thrustRequired)
hold on
plot(alt,thrustAvailable)
grid minor
title("Available and Required Thrust Vs. Altitude")
xlabel("Altitude [km]")
ylabel("Thrust [kN]")
legend("Required Thrust", "Available Thrust")

