%% Orbital Elements Function
% Author: Hudson Reynolds, Last Modified: 11/28/2024
%
% function converting orbital elements to initial cartesian states
%
% Inputs:
% a - semi major axis
% e - eccentricity
% i - inclination (relative to x-y plane)
% laan - longitude of ascending node [deg]
% aop - argument of periapsis [deg]
% mu - gravitational constant of parent body
% posInit - initial position of parent body 
%
% Outputs:
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