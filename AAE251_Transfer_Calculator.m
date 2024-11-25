%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: Transfer Orbit Calculations
% Author: Hudson Reynolds - Created: 11/24/2024
% Last modified: 11/25/2024
%
% Description: This is the overarching script that runs the satellite trajectory.
% The overarching simulation structure uses a custom RK4 structure to allow
% for independent calcuation of different bodies. The simulation uses real world 
% values for the orbits of the planets to simulate realistic orbits with full 
% 6DoF dynamics and N-body physics. The orbits of Venus, Earth, and the
% satellite are simulated.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization:
% clear console and figures
clear;
clc;
close all;

% conversions

AU2km = 1.496e8; % 1 AU in kilometers
s2yr = 3.1536e7; % number of seconds in year

% simulation timestep parameters
time = 0;      % simulation starting time [s]
dt = 3600;      % time step [s]
timespan = s2yr; % span of time to integrate over [s]

% Sun Inits: Initialize with no velocity, act as zero velocity reference
% x, y, z initial coordinates of sun:
sunPos = [0,0,0];
% velocity init:
sunVel = [0,0,0];
% gravitational constant for planet
sunMu = 1.3271e11;

% Earth Inits, pulled from J2000:
earthA = 1 * AU2km;     % semi major axis
earthE = 0.0167;        % eccentricity
earthI = 0.00005;       % inclination
earthRAAN = -11.26064;  % longitude of ascending node [deg]
earthAOP = 102.94719;   % argument of periapsis [deg]
earthMeanLong = 100.46435;
earthMu = 3.986e5;      % gravitational constant [km^3/s^2]

% calculate the cartesian state for the given inputs:
[earthPos, earthVel] = orbitalElements(earthA, earthE, earthI, earthRAAN, earthAOP, earthMeanLong, sunMu, sunPos);

% Venus Inits, pulled from J2000:
venusA = 0.723 * AU2km; % semi major axis
venusE = 0.00677;   % eccentricity
venusI = 3.39471;    % inclination
venusRAAN = 76.68069; % longitude of ascending node [deg]
venusAOP = 131.53298;  % argument of periapsis [deg]
venusMeanLong = 181.97973;
venusMu = 0.1;  % gravitational constant 
% calculate the cartesian state for the given inputs:
[venusPos, venusVel] = orbitalElements(venusA, venusE, venusI, venusRAAN, venusAOP, venusMeanLong, sunMu, sunPos);

% Satellite Inits:
satA = 6600;     % semi major axis
satE = 0;    % eccentricity
satI = 28.5;     % inclination
satRAAN = 50;  % longitude of ascending node [deg]
satAOP = 90;   % argument of periapsis [deg]
satMeanLong = 0;
% calculate the cartesian state for the given inputs:
[satPos, satVel] = orbitalElements(satA, satE, satI, satRAAN, satAOP, satMeanLong, sunMu, sunPos);

% manuver inits:
burnStartTime = 3.5; % start time of burn [s]
burnDuration = .2;   % duration of burn [s]
deltaV = .5;         % change in velocity of burn [-/s]

% initialize arrays:
sunPosArray = sunPos;
satPosArray = satPos;
earthPosArray = earthPos;
venusPosArray = venusPos;


%% Main simulation loop:

endTime = timespan / dt;

for i = 1:endTime
    t = dt * i;
    tArray(i) = t;

    % Planet Euler Method:
    sunPos = sunPos + sunVel * dt;
    sunPosArray = vertcat(sunPosArray, sunPos);

    %earth rk4 method:
    rEarth = sunPos - earthPos;
    % rk4 integrate:
    earthVel = rk4(@(t,y)accel(t,earthVel,rEarth, 0, 0, 0, sunMu, 0), dt, t, earthVel);
    earthPos = rk4(@(t,y)vel(t,earthVel), dt, t, earthPos);
    earthPosArray = vertcat(earthPosArray,earthPos);

    %venus rk4 method:
    rVenus = sunPos - venusPos;
    % rk4 integrate:
    venusVel = rk4(@(t,y)accel(t,venusVel,rVenus, 0, 0, 0, sunMu, 0), dt, t, venusVel);
    venusPos = rk4(@(t,y)vel(t,venusVel), dt, t, venusPos);
    venusPosArray = vertcat(venusPosArray,venusPos);

    %sat rk4 method:
    rSatPlanet = sunPos - satPos;
    rSatMoon = earthPos - satPos;
    % define vectors for multibody dynamics:
    rList = [rSatPlanet, rSatMoon];
    muList = [sunMu, earthMu];
    % rk4 integrate:
    satVel = rk4(@(t,y)accel(t,satVel,rList, burnStartTime, burnDuration, deltaV, muList, 1), dt, t, satVel);
    satPos = rk4(@(t,y)vel(t,satVel), dt, t, satPos);
    satPosArray = vertcat(satPosArray, satPos);
end


% Trajectory figure
figure(2)
plot3(sunPosArray(:,1), sunPosArray(:,2), sunPosArray(:,3), '.', 'MarkerSize', 20, 'Color', '#FFA500')
hold on
%plot3(satPosArray(:,1), satPosArray(:,2), satPosArray(:,3), 'LineWidth', 1, 'Color','#f97306')
plot3(earthPosArray(:,1), earthPosArray(:,2), earthPosArray(:,3), 'LineWidth', 1.5, 'Color', '#3b3b3b')
plot3(venusPosArray(:,1), venusPosArray(:,2), venusPosArray(:,3), 'LineWidth', 1.5, 'Color', 'red')
%plot3(satPosArray(burnIndex,1), satPosArray(burnIndex,2), satPosArray(burnIndex, 3), 'r*')
title('Planet Satellite Model 3D')
xlabel('Displacement-X')
ylabel('Displacement-Y')
zlabel('Displacement-Z')
legend('Sun', 'Earth', 'Venus')
hold off

xl = xlim;
yl = ylim;
zl = zlim;

axis equal

% Trajectory Animation:
figure(1)
curve = animatedline('LineWidth', 1.5, 'Color','#069af3');
curve1 = animatedline('LineWidth', 1, 'Color','#f97306');
curve2 = animatedline('LineWidth', 1.5, 'Color', '#3b3b3b');
xlabel('Displacement-X')
ylabel('Displacement-Y')
zlabel('Displacement-Z')
legend('Planet','Satellite', 'Moon')
axis equal
view(43,24);

output = 0;
playbackSpeed = 2;

xlim(xl);
ylim(yl);
zlim(zl);

for i=1:endTime
    addpoints(curve, sunPosArray(i,1), sunPosArray(i,2), sunPosArray(i,3));
    addpoints(curve1,earthPosArray(i,1), earthPosArray(i,2), earthPosArray(i,3));
    %addpoints(curve2,earthPosArray(i,1), earthPosArray(i,2), earthPosArray(i,3));
    title(sprintf("Satellite Trajectory at time %.2f s", tArray(i)));
    drawnow;
    %pause(dt / playbackSpeed);
    grid on;

    if output == 1
        frame = getframe(gcf);
        img =  frame2im(frame);
        [img,cmap] = rgb2ind(img,256);
        if i == 1
            imwrite(img,cmap,'TrajAnimation.gif','gif','LoopCount',Inf,'DelayTime',dt/playbackSpeed);
        else
            imwrite(img,cmap,'TrajAnimation.gif','gif','WriteMode','append','DelayTime',dt/playbackSpeed);
        end
    end
end






%% FUNCTIONS:

%% RK4 Integrator
% numerically integrates a first order initial value problem
% inputs:
% fun - function of integration
% dt - time step [s]
% tIn - input time [s]
% xIn - initial value input vector
% outputs:
% out - final value output vector
function out = rk4(fun, dt, tIn, xIn)
    f1 = fun(tIn,xIn);
    f2 = fun(tIn + dt/2, xIn + (dt/2) .* f1);
    f3 = fun(tIn + dt/2, xIn + (dt/2) .* f2);
    f4 = fun(tIn + dt, xIn + dt*f3);
    
    out = xIn + (dt / 6)*(f1 + 2*f2 + 2*f3+f4);
end
   
%% Acceleration Function
% function describing the translational accelerations on the satellite
% inputs:
% t- input time [s]
% v - current velocity [-/s]
% r - distance from bodies
% burnStart - start time of burn [s]
% burnDuration - duration of burn [s]
% dV - total delta V of the burn [-/s]
% mu - gravitational constant of bodies (must be in same order as r)
% burn - 0 => no burn compute, 1 => burn compute
% outputs:
% dv - acceleration vector for integration
function dv = accel(t, v, r, burnStart, burnDuration, dV, mu, burn)
    % gravity forces:
    gravAccel = 0;
    for i = 1:length(mu)
        r0 = (i - 1) * 3 + 1;
        rf =  i * 3;
        gravAccel = gravAccel + (mu(i) / (norm(r(r0:rf))^2) * (r(r0:rf) / norm(r(r0:rf))));
    end
    
    thrustAccel = 0;

    if burn == 1
        %thrust force:
        if t >= burnStart && t <= (burnStart + burnDuration)
            % apply force along the velocity vector
            thrustAccel = (v / norm(v)) * (dV / burnDuration);
        end
    end
    dv = gravAccel + thrustAccel;
end

%% Velocity Function
% function describing the position-velocity on the satellite
% inputs:
% t- input time [s]
% v - current velocity [-/s]
% outputs:
% dr - velocity vector for integration
function dr = vel(t, v)
dr = v;
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
function [pos, vel] = orbitalElements(a, e, i, raan, aop, meanLong, mu, posInit)

% convert inputs to radians:
raan = deg2rad(raan);
i = deg2rad(i);
aop = deg2rad(aop);

% find the rotation matrix wrt to the planet centered frame:
vecRotMat = eul2rotm([raan,i,aop], 'ZXZ');

% find the vector pointing towards the periapsis
vec = vecRotMat * [1;0;0];

% find the length of the periapsis from orbital elements:
periapsis = a * (1 - e);

% find the initial position (add parent body origin to get absolute origin)
pos = vec.' * periapsis + posInit;

% get velocity from vis-viva
velMag = sqrt(mu * ((2/periapsis) - (1 / a)));

% find the velocity vector:
% find another vector 90 degrees along orbit plane:
vector2RotMat = eul2rotm([raan,i,aop + pi/2], 'ZXZ');
vec2 = vector2RotMat * [1;0;0];

% find orbit normal (z-axis of orbit):
orbitNorm = cross(vec2, vec) / norm(cross(vec2, vec));

% find orbit velocity direction w/ cross product:
velDir = cross(orbitNorm, vec) / norm(cross(orbitNorm, vec));
vel = (velMag * velDir).';
end