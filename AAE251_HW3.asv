%AAE 251 Fall 2024
%Homework 3
%AAE251_HW3
%Authors: Hudson Reynolds

%% Initializations:
% Constants Definitions:
muVenus = 3.2486e5; % G * Mass of venus [km^3 / s^2]
rVenus = 6.052e3;    % radius of venus [km]
trueAnomoly = linspace(0,2 * pi, 1001);

% Parent orbit definitions:
aParent = rVenus + 325; % semi-major axis of parent orbit

% Child orbit definitions:
aChild = 6.5645e3;      % semi-major axis of child orbit
apRad = rVenus + 700;   % apoapsis of child orbit
e = 0.0285627;          % eccentricity of child orbit

%% Calculations;

% Parent orbit calc:
rParent = aParent;

% Child orbit calc:
p = aChild * (1-e^2);
rChild = p ./ (1 + e * cos(trueAnomoly));

%% Plotting:
% convert to polar coords

%venus surface
xVenus = rVenus * cos(trueAnomoly);
yVenus = rVenus * sin(trueAnomoly);

%parent orbit:
xParent = rParent * cos(trueAnomoly);
yParent = rParent * sin(trueAnomoly);

%child orbit:
for index = 1:length(trueAnomoly)
    xChild(index) = rChild(index) * cos(trueAnomoly(index));
    yChild(index) = rChild(index) * sin(trueAnomoly(index));
end

% make figure

figure(1)
hold on
plot(xVenus, yVenus, '-', 'LineWidth', 1, 'Color', '#dbb40c')
plot(xParent, yParent, 'LineWidth', 0.75,'Color', 'red')
plot(xChild, yChild, 'LineWidth', 0.75, 'Color', 'blue')
grid on
axis square
title("Parent and Child Orbits around Venus")
xlabel("Position in x direction [km]")
ylabel("Position in y direction [km]")
legend("Venus Surface", "Parent Orbit", "Child Orbit", location='southwest')


