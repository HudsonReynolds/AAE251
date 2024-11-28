%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: AAE251 Porkchop
% Author: Hudson Reynolds - Created: 11/27/2024
% Last modified: 11/27/2024
%
% Description: This is the script that performs the calculations for the
% tranfer between planets. This is done by solving Lambert's problem
% between Earth and Venus using a patched conics model. This script outputs
% the delta-V's for each burn and the 

% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization:
% clear console and figures
clear;
clc;
close all;

% turn plots on and off

plotnum1 = 0;
plotnum2 = 1;
plotnum3 = 0;

% conversions

AU2km = 1.496e8; % 1 AU in kilometers
s2yr = 3.1536e7; % number of seconds in year
days2s = 86400;
G = 6.6743e-11;

% simulation timestep parameters
time = 0;      % simulation starting time [s]
dt = 100;      % time step [s]
timespan = (s2yr/365) * 1; % span of time to integrate over [s]

% Sun Inits: Initialize with no velocity, act as zero velocity reference
% x, y, z initial coordinates of sun:
sunPos = [0,0,0];
% velocity init:
sunVel = [0,0,0];
% gravitational constant for planet
sunMu = 1.3271e11;
sunMass = 1.989e30;

% Earth Inits, pulled from J2000 on Nov 15, 2024:
t1Earth = datetime('2024-November-15')

earthA = 1 * AU2km;     % semi major axis
earthE = 0.0167;        % eccentricity
earthI = 0.00005;       % inclination
earthRAAN = -11.26064;  % longitude of ascending node [deg]
earthAOP = 102.94719;   % argument of periapsis [deg]
earthMeanLong = 100.46435;
earthMu = 3.986e5;      % gravitational constant [km^3/s^2]
earthMass = 5.9722e24;  % mass of Earth [kg]

earthPeriod = s2yr;     % earth period in seconds

% calculate the cartesian state for the given inputs:
[earthPos, earthVel] = orbitalElements(earthA, earthE, earthI, earthRAAN, earthAOP, earthMeanLong, sunMu, sunPos, 0);

% Venus Inits, pulled from J2000 on Jan 11, 2024:

t1Venus = datetime('2024-January-11');

venusA = 0.723 * AU2km; % semi major axis
venusE = 0.00677;       % eccentricity
venusI = 3.39471;       % inclination
venusRAAN = 76.68069;   % longitude of ascending node [deg]
venusAOP = 131.53298;   % argument of periapsis [deg]
venusMeanLong = 181.97973;
venusMu = 3.248e5;  % gravitational constant 

venusPeroid = s2yr * (224.7/365.25);

% calculate the cartesian state for the given inputs:
[venusPos, venusVel] = orbitalElements(venusA, venusE, venusI, venusRAAN, venusAOP, venusMeanLong, sunMu, sunPos, 0);


t2 = datetime('2025-March-12')

timeDiffEarth = between(t1Earth,t2);

TransferTime = 147 * days2s;

r0 = earthPos

rF = venusPos


[v0, vF] = lambertProblem(earthPos, venusPos, TransferTime, 1, sunMu, 1e-6, 500, 0, 4*pi^2, -4*pi)


% perform the lambert calcs

function [v0, vF] = lambertProblem(r0, rF, dT, tm, mu, tolerance, maxSteps, psi, psi_ub, psi_lb)

    sqrtMu = sqrt(mu);
    r0Norm = norm(r0);
    rFNorm = norm(rF);

    gamma = dot(r0, rF) / (r0Norm * rFNorm);

    beta = tm * sqrt(1 - gamma^2);

    A = tm * sqrt(r0Norm* rFNorm * (1 + gamma));

    if A == 0
        A = [0,0,0];
    end

    c2 = 0.5;
    c3 = 1/6;

    for i = 1:maxSteps
        B = r0Norm + rFNorm * (psi * c3 - 1) / sqrt(c2);

        if A > 0 & B < 0
            psi_lb = psi_lb + pi;

            B = B * -1;
        end

        chi3 = sqrt(B / c2)^3;

        deltaT = (chi3*c3+A*sqrt(B)) / sqrtMu;

        if abs(dT - deltaT) < tolerance

            break
        end

        if deltaT <= dT

            psi_lb = psi;
        else
            psi_ub = psi;
        end

        psi = (psi_ub + psi_lb) / 2;

        c2 = (1 - cos(sqrt(psi))) / psi;

        c3 = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi));

        f = 1 - B / r0Norm;

        g = A * sqrt(B/mu);

        gdot = 1 - B/rFNorm;

        v0 = (rF - f*r0)/g;
        vF = (gdot*rF-r0)/g;

    end
end

%% Orbital Elements Function
% function converting orbital elements to initial cartesian states
% inputs:
% a - semi major axis
% e - eccentricity
% i - inclination (relative to x-y plane)
% laan - longitude of ascending node [deg]
% aop - argument of periapsis [deg]
% mu - gravitational constant of parent body
% posInit - initial position of parent body 
% outputs:
% pos - initial position vector
% vel - initial velocity vector
function [pos, vel] = orbitalElements(a, e, i, raan, lop, meanLong, mu, posInit, velInit)

% convert inputs to radians:
raan = deg2rad(raan);
i = deg2rad(i);
lop = deg2rad(lop);
meanLong = deg2rad(meanLong);

aop = lop - raan;

meanAnomoly = meanLong - raan - aop;

% find the rotation matrix wrt to the planet centered frame:
vecRotMat = eul2rotm([raan,i,aop], 'ZXZ');

vecRotMat = eul2rotm([raan, i, aop + meanAnomoly], 'ZXZ');

% find the vector pointing towards the periapsis
vec = vecRotMat * [1;0;0];

% find the length of the periapsis from orbital elements:
periapsis = a * (1 - e);

r = a * (1 - e^2) / (1 + e * cos(meanAnomoly));

% find the initial position (add parent body origin to get absolute origin)
pos = vec.' * periapsis + posInit;

pos = vec.' * r + posInit;

% get velocity from vis-viva
velMag = sqrt(mu * ((2/periapsis) - (1 / a)));

velMag = sqrt(mu * ((2/r) - (1 / a)));

% find the velocity vector:
% find another vector 90 degrees along orbit plane:
vector2RotMat = eul2rotm([raan,i,aop + pi/2], 'ZXZ');
vec2 = vector2RotMat * [1;0;0];

% find orbit normal (z-axis of orbit):
orbitNorm = cross(vec2, vec) / norm(cross(vec2, vec));

% find orbit velocity direction w/ cross product:
velDir = cross(orbitNorm, vec) / norm(cross(orbitNorm, vec));
vel = (velMag * velDir).' + velInit;
end

