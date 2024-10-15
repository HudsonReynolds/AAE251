%Homework 5
%AAE251_HW5
%Authors: Hudson Reynolds & Preston Wright

function AAE251_HW5(NACA)
% This function takes a 4-digit NACA profile as an input and returns the
% associated geometry of the profile in a plot.
% 
% Inputs: NACA - 4 digit NACA profile
%
% Outputs: Plot
%
%AAE 251 Fall 2024





%% NACA 2415 profile:

% NACA follows m p tt convention with numbers

NACA = 2415;

NACAString = num2str(NACA);

% m = maximum camber percentage
m = str2num(NACAString(1)) / 100;
% p = location of max camber out of 10
p = str2num(NACAString(2)) / 10;
% tt = thickness to chord ratio percentage
t = str2num(NACAString(3:4)) / 100;


%camber line equation

% define chord length as 1
x = linspace(0, 1, 401);

% find the equation of the line of camber:
for index = 1:length(x)
    if x < p
        camberLine(index) = m * (x(index) / p^2) * (2 * p - x(index));
    else
        camberLine(index) = m * ((1 - x(index)) / (1 - p)^2) * (1 + x(index) - 2 * p);
    end
end

%take the derivative wrt to x to find the angle at any point:

for index = 1:length(x)
    if x < p
        camberLineAngle(index) = (2 * m / p^2) * (2 * p - x(index));
    else
        camberLineAngle(index) = (2 * m / (1 - p)^2) * (p - x(index));
    end
end

theta = atan(camberLineAngle);

% find the thickness at any given point:

thickness = 5 * t * (0.2969 * x.^(1/2) - 0.1260 .* x - 0.3516 * x.^2 + 0.2843 .* x.^3 - 0.1015 * x.^4);

% find the equations of the surfaces:

xTop = x - thickness .* sin(theta);
xLow = x + thickness .* sin(theta);

yTop = camberLine + thickness .* cos(theta);
yLow = camberLine - thickness .* cos(theta);

% plot the wing 

figure(1)
hold on
plot(x,camberLine, '--', 'Color', 'red');
plot(xTop, yTop, 'Color', 'blue', 'LineWidth', 1)
plot(xLow, yLow, 'Color', 'blue', 'LineWidth', 1)
title(sprintf('Airfoil Geometry for NACA %s',  NACAString))
ylabel('Thickness')
xlabel('Chord %');
axis equal;
legend('Mean Camber Line', 'Airfoil')



