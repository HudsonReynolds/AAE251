% AAE251 Fall 2024
% Homework 7
% AAE251_HW7_Q2a
% Author: Preston Wright and Hudson Reynolds

%% Initializations 

isp2List = [356, 372, 383, 399, 445]; % list of possible isps [s]
alphaList = 0.30:0.001:0.5;           % list of alpha values
fInert = 1/12;                        % Inert mass fraction
g = 9.81;                             % Gravity
isp1 = 295;                           % Specific impulse of stage one
deltaVTot = 11.85*1000;               % Total delta V required
mPayload = 99;                        % Payload mass
i = 1;                                % Iteration counter

%% Calculations

% The initial Gross Lift Off Weight (GLOW) is a function of the delta-v
% proportion for stage 1 and stage 2, alpha 1 and alpha 2, respectively.
% These proportions add to one. To find the ideal value, the value of alpha
% 1 is iterated over in the span from 0.5 to 0.67.

% Initial mass is determined from a variation of the Tsiolkovsky Rocket
% equation for two stages

for k = alphaList

    % Find current value of alpha one and two
    alpha1 = alphaList(i);
    alphaOneArray(i) = alphaList(i);
    alpha2 = 1 - alpha1;

    % Calculate delta V and the mass fraction
    deltaV1 = alpha1 * deltaVTot;
    deltaV2 = alpha2 * deltaVTot;
    mFracOne = exp(deltaV1/(isp1*g));
    for j = 1:length(isp2List)
        mFracTwo(j) = exp(deltaV2/(isp2List(j)*g));
        % Calculate initial mass
        mInit(i, j) = mPayload * ((mFracOne*(1-fInert))/(1-fInert*mFracOne)) * ((mFracTwo*(1-fInert))/(1-fInert*mFracTwo));
    end

    % Update array index
    i = i + 1;
end

figure(1)

plot(alphaList, mInit)
xlim([0.33 0.5])

legend("356 s", "372 s", "383 s", "399 s", "445 s")
