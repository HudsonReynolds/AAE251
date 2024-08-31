%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 251 
% Program Description 
% This program calculates temperature, pressure, density, and speed of
% sound for a given input height. Data from Anderson Introduction to Flight
% 8th Editon.
%
% Assignment Information
%   Assignment:     HW1, Section B
%   Author:         Hudson Reynolds, reyno250@purdue.edu
%   Team ID:        R406
%   Academic Integrity:
%     [X] I worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers I worked with: Preston Wright, wrigh735@purdue.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init:
% create an array of the heights
h = linspace(0,105000,106);

% get the values from the function
[T] = AtmoModel(h,0);

%% Plotting
figure(1)
plot(T,h)


%% FUNCTION DEFINITIONS:
function[T,P,rho,A] = AtmoModel(height, units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function Description:
% This function calculates temperature, pressure, density, and speed of
% sound for a given input height. Data from Anderson Introduction to Flight
% 8th Editon.
%
% Inputs:
% height = desired height [km].
% units = 0 => SI, 1 => English
%
% Outputs: 
%
% Outputs shown [SI / English]
%
% T= temperature [K / F]
% P = pressure [Pa / psi]
% rho = density [kg * m^-3 / slug * ft^-3]
% A = speed of sound [m * s^-1 / mph]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start by calculating everything in SI:

% Start with temperature, as Anderson does. This is defined from
% experimental data so we have no derivation, just input:

% temperatures in K
T0 = 288.16;
T1 = 216.66;
T4 = 282.66;
T6 = 165.66;
T7 = 225.66;

% Altitudes in km:
A0 = 0;
A1 = 11000;
A2 = 25000;
A3 = 47000;
A4 = 53000;
A5 = 79000;
A6 = 90000;
A7 = 105000;

% slopes of lines:
m1 = -6.5e-3;
m2 = 3e-3;
m3 = -4.5e-3;
m4 = 4e-3;

% pressure

% input into for loop to find all values:
for i = 0:length(height)-1
    i = i + 1;
    if height(i) < 0
        fprintf('Function Error! Exceeded Lower Bounds!\n');
    elseif height(i) > A0 & height(i) < A1 
        T(i) = T0 + m1 * height(i);
    elseif height(i) >= A1 & height(i) <= A2
        T(i) = T1;
    elseif height(i) > A2 & height(i) < A3
        T(i) = T1 + m2 * (height(i) - A2);
    elseif height(i) >= A3 & height(i) <= A4
        T(i) = T4;
    elseif height(i) > A4 & height(i) < A5
        T(i) = T4 + m3 * (height(i) - A4);
    elseif height(i) >= A5 & height(i) <= A6
        T(i) = T6;
    elseif height(i) > A6 & height(i) < A7
        T(i) = T6 + m4 * (height(i) - A6);
    end

end

% Set initial and final values:
T(1) = T0;
T(106) = T7;

end








