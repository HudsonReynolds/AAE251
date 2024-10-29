% AAE251 Fall 2024
% Homework 7
% AAE251_HW7_Q2a
% Author: Preston Wright and Hudson Reynolds
%% Initializations

fInert = 1/12;              % Inert mass fraction
g = 9.81;                   % Gravity
specImpOne = 295;           % Specific impulse of stage one
specImpTwo = 355;           % Specific impulse of stage two
deltaVTot = 11.85*1000;     % Total delta V required
mPayload = 99;              % Payload mass
tests = 10000;              % Total iterations
i = 1;                      % Iteration counter

%% Calculations

% Iterate from .3 to .6 for alpha one
for k = .3*tests:tests-.4*tests

    % Find current value of alpha one and two
    alphaOne = k/tests;
    alphaOneArray(i) = alphaOne;
    alphaTwo = 1 - alphaOne;

    % Calculate delta V and the mass fraction
    deltaVOne = alphaOne * deltaVTot;
    deltaVTwo = alphaTwo * deltaVTot;
    mFracOne = exp(deltaVOne/(specImpOne*g));
    mFracTwo = exp(deltaVTwo/(specImpTwo*g));

    % Calculate initial mass
    mInit(i) = mPayload * ((mFracOne*(1-fInert))/(1-fInert*mFracOne)) * ((mFracTwo*(1-fInert))/(1-fInert*mFracTwo));

    % Update array index
    i = i + 1;
end

% Find the location of the minimum initial mass and find the corresponding
% alpha one value
[idealMInit,idealAlphaOneLoc] = min(mInit);
idealAlphaOne = alphaOneArray(idealAlphaOneLoc);

% Calculate corresponding alpha two, final mass, and propellant mass
idealAlphaTwo = 1 - idealAlphaOne;
idealMFinal = idealMInit / mFracOne;
idealMProp = idealMInit - idealMFinal;

%% Graphing

% Plot GLOW vs. Alpha One, with a marker at the minimum value
figure(1)
title("GLOW vs. Alpha One")
plot(alphaOneArray, mInit)
hold on
plot(idealAlphaOne, idealMInit, Marker="x", MarkerSize=10)
xlabel("Alpha One")
ylabel("GLOW")
grid on

%% Outputs

% Output necessary information
fprintf("The alpha one value that leads to a minimum GLOW value is: %.4f\n", idealAlphaOne)
fprintf("The corresponding alpha two value is: %.4f\n", idealAlphaTwo)
fprintf("The corresponding propellant mass is: %.4f\n", idealMProp)

