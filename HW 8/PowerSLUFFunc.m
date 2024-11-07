function [power] = PowerSLUFFunc(rho, vel, area, cD0, spanEfficiency, mass)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Preston Wright, Hudson Reynolds
% Description: function that finds the thrust at SLUF for given conditions
% and outputs the calculated thrust
%
% Inputs:
% area - wing / reference area [m^2]
% rho - density of air [kg/m^3]
% cD0 - coefficient of drag at zero AoA
% vel - velocity [m/s]
% mass - mass of aircraft [kg]
%
% Outputs:
% power - calculated thrust for given inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations
g = 9.81;

%% Calculations
power = ((1/2) * rho * vel^3 * area * cD0) + ...
    2 * spanEfficiency * ((mass * g)^2/(rho * area * vel));