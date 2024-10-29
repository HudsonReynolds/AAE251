%AAE 251 Fall 2024
%Homework 4
%AAE251_HW4_Q2
%Authors: Hudson Reynolds and Preston Wright

%% Initializations:
% constants:
g0 = 9.80665; % gravitational constant [m s^-2]

% write the initializations for the two stage rocket
dVTot = 11850; % total delta-V [m/s]
mPay = 99;     % payload mass [kg]
fInert = 1/12; % inert mass fraction for both stages

% define parameters for the first stage:
IspStage1 = 295;     % Isp of first stage [s]
C1 = IspStage1 * g0; % exhaust velocity of stage 1 [m/s]

% define parameters for the second stage:
IspStage2 = 355;     % Isp of second stage [s]
C2 = IspStage2 * g0; % exhaust velocity of stage 2 [m/s]

%% Calculations:

% The initial Gross Lift Off Weight (GLOW) is a function of the delta-v
% proportion for stage 1 and stage 2, alpha 1 and alpha 2, respectively.
% These proportions add to one. To find the ideal value, the value of alpha
% 1 is iterated over in the span from 0.1 to 0.8.

% Initial mass is determined from a variation of the Tsiolkovsky Rocket
% equation for two stages

% create the for loop for iteration:

i = 0;

alphaList = 0.3:0.0001:0.6;

for alpha1 = alphaList

    % create an array index
    i = i + 1;

    % write the delta v of each stage in terms of the alpha variable
    dV1 = alpha1 * dVTot;
    dV2 = (1 - alpha1) * dVTot;
    
    % calculate the mass fraction of each stage:
    MR1 = exp(dV1 / C1);
    MR2 = exp(dV2 / C2);


    % figure out the initial mass:

    mI(i) = mPay * ((MR1 * (1-fInert)) / (1-fInert*MR1)) * ((MR2*(1-fInert)) / (1-fInert*MR2));

    % mI(i) = mPay * (MR2 * (1 - fInert) * MR1 * (1 - fInert)) / ((1 - fInert * MR2) * (1 - fInert * MR1));

end

figure(1)
plot(alphaList, mI)


