%AAE 251 Fall 2024
%Homework 4
%AAE251_HW4
%Authors: Hudson Reynolds and Preston Wright

%% Initializations:
% values for earth and orbit:
OmegaEarth = 7.292e-5; % sidereal rate [rad s^-1]
radEarth = 6371;    % radius of earth [km]
muEarth = 3.986e5;     % gravitational parameter [km^3 s^-2]
dVLoss = 1.7;          % delta-V loss [km s^-1]
radOrbit = 290;        % orbital altitude [km]

%kourou inputs:
latKourou = [5,3,27];     % latitude [degree, min, sec]
azimuthsKourou = [-10.5,93.5]; % permissable launch azi [min, max]

% KSC inputs:
latKSC = [28,31,27];      % latitude [degree, min, sec]
azimuthsKSC = [35, 120];       % permissable launch azi [min, max]


%% Part A Calculations:

% find latitudes in degrees:
latKourou = latKourou(1) + (latKourou(2)/60) + (latKourou(3)/3600);

latKSC = latKSC(1) + (latKSC(2)/60) + (latKSC(3)/3600);

% delta-v added by the earth for azimuths -20 to 130 degrees:

n = 1;

for Az = -20:130
    dVEarthHelpKourou(n) = OmegaEarth * radEarth * cosd(latKourou) * sind(Az);
    
    dVEarthHelpKSC(n) = OmegaEarth * radEarth * cosd(latKSC) * sind(Az);

    AzList(n) = Az;
    n = n + 1;
end

dVRangeKourou = OmegaEarth * radEarth * cosd(latKourou) * sind(90) ...
    - OmegaEarth * radEarth * cosd(latKourou) * sind(azimuthsKourou(1));

dVRangeKSC = OmegaEarth * radEarth * cosd(latKSC) * sind(90) ...
    - OmegaEarth * radEarth * cosd(latKSC) * sind(azimuthsKSC(1));

dVMax = OmegaEarth * radEarth * cosd(latKourou) * sind(90);

% plot figure:

figure(1)
subplot(2,1,1)
plot(AzList, dVEarthHelpKSC, 'Color', "red", "LineWidth", 1)
xlabel('Launch Azimuth Angle (deg)')
xline(azimuthsKSC(1), '--', 'Minimum Launch Azimuth', 'LabelVerticalAlignment','bottom')
xline(azimuthsKSC(2), '--', 'Maximum Launch Azimuth', 'LabelVerticalAlignment','bottom')
ylabel('Delta-V (km/s)')
title('Delta-V Imparted by Earth at KSC by Launch Azimuth')
grid on

subplot(2,1,2)
plot(AzList, dVEarthHelpKourou, "Color", "blue", "LineWidth", 1)
xlabel('Launch Azimuth Angle [deg]')
xline(azimuthsKourou(1), '--', 'Minimum Launch Azimuth', 'LabelVerticalAlignment','top')
xline(azimuthsKourou(2), '--', 'Maximum Launch Azimuth', 'LabelVerticalAlignment','bottom')
ylabel('ΔV Imparted by Earth [km/s]')
title('ΔV Imparted by Earth at Kourou by Launch Azimuth')
grid on;

% print statements:

fprintf("\nThe Range of ΔV Earth Help from Kourou is %.2f km/s.\n", dVRangeKourou);
fprintf("The Range of ΔV Earth Help from KSC is %.2f km/s.\n", dVRangeKSC);

fprintf("\nThe Maximum ΔV Earth Help of %.2f km/s occurs at Kourou at 90 degrees launch azimuth.\n", dVMax);

% sin(Az) = cos(i) / cos(La)
% i >= La

%% Part B Calculations:

% final inclinations:
InitAziKourou = [-10,35,90];
InitAziKSC = [40,90,110];

% v needed for low earth orbit:
dVLEO = sqrt(muEarth / (radOrbit + radEarth));

%calculate the total energy to get into orbit from certain azimuth:
dVTotKSC = dVLEO + dVLoss - (OmegaEarth * radEarth * cosd(latKSC) * sind(InitAziKSC))
dVTotKourou = dVLEO + dVLoss - (OmegaEarth * radEarth * cosd(latKourou) * sind(InitAziKourou))
dVTotCombined = [dVTotKSC, dVTotKourou];

% calculate the final inclination based on the launch location and azimuth:
orbitalIncKSC = acosd(sind(InitAziKSC) * cosd(latKSC));
orbitalIncKourou = acosd(sind(InitAziKourou) * cosd(latKourou));

orbitalIncCombined = [orbitalIncKSC, orbitalIncKourou];

n = 1;

for index = 1:length(orbitalIncCombined)
    for finalAzimuth = 0:0.1:90
        nu = abs(orbitalIncCombined(index) - finalAzimuth);
        dVIncChange = 2 * dVLEO * sind(nu);
        finalAzimuthList(n) = finalAzimuth;
        dVTot(index,n) = dVIncChange + dVTotCombined(index);
        n = n + 1;
    end
    n = 1;
end


figure(2)
plot(finalAzimuthList, dVTot, 'LineWidth', 1)
title('ΔV requirements for Inclination Changes 0 to 90 degrees from KSC and Kourou');
xlabel('Final Orbit Inclination [deg]');
ylabel('ΔV Requirement [km/s]')
grid on
legend('KSC 40 deg Azimuth', 'KSC 90 deg Azimuth', 'KSC 110 deg Azimuth',...
    'Kourou -10 deg Azimuth', 'Kourou 35 deg Azimuth', 'Kourou 90 deg Azimuth', Location='best')



 