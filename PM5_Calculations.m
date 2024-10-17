%AAE 251 Fall 2024
%PM5
%PM5_Calculations
%Authors: Hudson Reynolds

%% PART 1:

%% Initializations:
% values for earth and orbit:
OmegaEarth = 7.292e-5; % sidereal rate [rad s^-1]
radEarth = 6371;       % radius of earth [km]
muEarth = 3.986e5;     % gravitational parameter [km^3 s^-2]
dVLoss = 1.65;         % delta-V loss [km s^-1]
radOrbit = 200;        % orbital altitude [km]

% KSC inputs:
latKSC = [28,31,27];      % latitude [degree, min, sec]
latKSC = latKSC(1) + (latKSC(2)/60) + (latKSC(3)/3600);
Az = 90;

dVEarthHelp = OmegaEarth * radEarth * cosd(latKSC) * sind(Az)

dVLEO = sqrt(muEarth / (radOrbit + radEarth))

dVTot = dVLEO + dVLoss - dVEarthHelp

%% PART 2:

%% Orbits start on March 23, 2025 in an inferior conjuction.    

% inits:

muSun = 1.33e11;      % mu value of sun [km^3 s^-2]

aEarthOrbit = 1.5e8;  % radius of earth orbit [km]

aVenusOrbit = 1.08e8; % radius of venus orbit [km]

aTransferOrbit = (aEarthOrbit + aVenusOrbit) / 2;

periodEarth = 2 * pi * sqrt(aEarthOrbit^3 / muSun);

periodEarthYrs = periodEarth / (60 * 60 * 24 * 365);

periodVenus = 2 * pi * sqrt(aVenusOrbit^3 / muSun);

transferTime = pi * sqrt(aTransferOrbit^3 / muSun);

transferTimeDays = transferTime / (60 * 60 * 24);

meanMotionEarth = 2 * pi / periodEarth;

trueAnomolyEarth = meanMotionEarth * t;



