function [thrust] = ThrustSLUFFunc(rho, vel, area, cD0, spanEfficiency, mass)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Preston Wright, Hudson Reynolds
% Description: function that finds results from the lift and drag equation
% assuming SLUF flight conditions
% and outputs the forces
% area - wing / reference area [m^2]
% rho - density of air [kg/m^3]
% cD0 - coefficient of drag at zero AoA
% vel - velocity [m/s]
% mass - mass of aircraft [kg]
%
% Outputs:
% thrust - calculated thrust for given inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializations
g = 9.81;

%% Calculations 
thrust = ((1/2) * rho * vel^2 * area * cD0) + ...
    2 * spanEfficiency * ((mass * g)^2/(rho * area * vel^2));