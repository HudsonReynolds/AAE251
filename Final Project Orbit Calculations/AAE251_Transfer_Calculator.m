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

% turn plots on and off

plotnum1 = 0;
plotnum2 = 1;
plotnum3 = 0;
plotnum4 = 0;

% conversions

AU2km = 1.496e8; % 1 AU in kilometers
yr2s = 3.1536e7; % number of seconds in year
G = 6.6743e-11;  % gravitational constant
days2s = 86400;  % number of seconds in a day


% simulation time parameters
time = 0;      % simulation starting time [s]
dt = 600;      % time step [s]

transferTime = 147 * days2s; % amount of time for the transfer

t2 = datetime('15-11-2025', 'InputFormat','dd-MM-yyyy'); % date at which the simulation starts

tInt = t2 + seconds(transferTime);

timespan = transferTime; % span of time to integrate over [s]

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

% use the starting mean longitude to find the location at the start date:
earthMeanLong0 = 100.46435; 
earthPeriod = 365.256 * days2s;
earthMeanMotion = 2*pi / earthPeriod;

earthMeanLong = meanLongSolver(t1Earth, t2, earthMeanMotion, earthMeanLong0)

earthA = 1 * AU2km;     % semi major axis
earthE = 0.0167;        % eccentricity
earthI = 0.00005;       % inclination
earthRAAN = -11.26064;  % longitude of ascending node [deg]
earthAOP = 102.94719;   % argument of periapsis [deg]
earthMu = 3.986e5;      % gravitational constant [km^3/s^2]
earthMass = 5.9722e24;  % mass of Earth [kg]

% calculate the cartesian state for the given inputs:
[earthPos, earthVel] = orbitalElements(earthA, earthE, earthI, earthRAAN, earthAOP, earthMeanLong, sunMu, sunPos, 0)

[earthPos, earthVel] = planetEphemeris(juliandate(t2), 'Sun', 'Earth')

% Venus Inits, pulled from J2000 on Jan 11, 2024:
t1Venus = datetime('2024-January-11');
venusMeanLong0 = 181.97973;
venusPeriod = 224.7 * days2s;
venusMeanMotion = 2*pi / venusPeriod;

venusMeanLong = meanLongSolver(t1Venus, t2, venusMeanMotion, venusMeanLong0)
venusMeanLongIntercept = meanLongSolver(t1Venus, tInt, venusMeanMotion, venusMeanLong0)

venusA = 0.723 * AU2km; % semi major axis
venusE = 0.00677;       % eccentricity
venusI = 3.39471;       % inclination
venusRAAN = 76.68069;   % longitude of ascending node [deg]
venusAOP = 131.53298;   % argument of periapsis [deg]
venusMu = 3.248e5;  % gravitational constant 
% calculate the cartesian state for the given inputs:
[venusPos, venusVel] = orbitalElements(venusA, venusE, venusI, venusRAAN, venusAOP, venusMeanLong, sunMu, sunPos, 0);

[venusPos, venusVel] = planetEphemeris(juliandate(t2), 'Sun', 'Venus')

endTime = timespan / dt;

%prealloc arrays for speed:

sunPosArray = zeros(timespan / dt, 3);
earthPosArray = zeros(timespan / dt, 3);
earthVelArray = zeros(timespan / dt, 3);
venusPosArray = zeros(timespan / dt, 3);
satPosArray = zeros(timespan / dt, 3);
satVelArray = zeros(timespan / dt, 3);
tArray = zeros(timespan/dt , 1);

venusPosIntercept = planetEphemeris(juliandate(t2) + transferTime / days2s, 'Sun', 'Venus')

% Satellite Inits:
% satA = 6600;     % semi major axis
% satE = 0;    % eccentricity
% satI = 28;     % inclination
% satRAAN = 30;  % longitude of ascending node [deg]
% satAOP = 90;   % argument of periapsis [deg]
% satMeanLong = 10;
% % calculate the cartesian state for the given inputs:
% [satPos, satVel] = orbitalElements(satA, satE, satI, satRAAN, satAOP, satMeanLong, earthMu, earthPos, earthVel);
% 
% % manuver inits:
% burn1StartTime = 0; % start time of burn [s]
% dV1 = 3.276;         % change in velocity of burn [km/s]


[v0, vF] = lambertProblem(earthPos, venusPosIntercept, transferTime, 1, sunMu, 1e-6, 500, 0, 4*pi^2, -4*pi)

satPos = earthPos;
satVel =  v0;



%% Main simulation loop:

for i = 1:endTime
    t = dt * i;
    tArray(i) = t;

    % Sun Euler Method:
    sunPos = sunPos + sunVel * dt;
    sunPosArray(i,:) = sunPos;

    %earth rk4 method:
    rEarth = sunPos - earthPos;
    % rk4 integrate:
    earthVel = rk4(@(t,y)accel(t,earthVel,rEarth,sunMu), dt, t, earthVel);
    earthPos = rk4(@(t,y)vel(t,earthVel), dt, t, earthPos);
    earthPosArray(i,:) = earthPos;
    earthVelArray(i,:) = earthVel;

    %venus rk4 method:
    rVenus = sunPos - venusPos;
    % rk4 integrate:
    venusVel = rk4(@(t,y)accel(t,venusVel,rVenus, sunMu), dt, t, venusVel);
    venusPos = rk4(@(t,y)vel(t,venusVel), dt, t, venusPos);
    %venusPosArray(i,:) = venusPos;

    %sat rk4 method:
    rSatEarth = earthPos - satPos;
    rSatVenus = venusPos - satPos;
    rSatSun = sunPos - satPos;
    % define vectors for multibody dynamics:
    % rList = [rSatEarth, rSatVenus, rSatSun];
    % muList = [earthMu, venusMu, sunMu];
    rList = [rSatVenus, rSatSun];
    muList = [venusMu, sunMu];
    % rk4 integrate:
    satVel = rk4(@(t,y)accel(t,satVel,rList, muList), dt, t, satVel);

    thrustAccel = 0;

    % %thrust force:
    % if t >= burn1StartTime && t < burn1StartTime + dt
    %     % apply force along the velocity vector
    %     deltaV1 = satVel / norm(satVel) * dV1
    % 
    %     satVel = satVel + deltaV1
    % end


    satPos = rk4(@(t,y)vel(t,satVel), dt, t, satPos);
    satPosArray(i,:) = satPos;
    satVelArray(i,:) = satVel;
end


% Energy of System
if plotnum1 == 1
    % find the kinetic energy, 1/2 m v^2
    earthKE = 1/2 * vecnorm(earthVelArray * 1000, 2, 2).^2 * earthMass;
    % find the potential energy, mg
    earthPE = G * earthMass * sunMass ./ (vecnorm(earthPosArray, 2, 2) * 1000);
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'Total Mechanical Energy';

    hold on   

    energy = (earthKE - earthPE);

    
    plot(tArray / yr2s, energy / 1e9, 'LineWidth', 2)

    title("Total Mechanical Energy for Earth")
    xlabel("Time [yrs]")
    ylabel("Total Energy [GJ]")
    
    
    picturewidth = 20; % set the width of image in cm
    hw_ratio = .6; % aspect ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust font size
    
    grid on
    
    set(findall(hfig,'-property','Box'),'Box','off') % turn off box
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    


    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-vector','-fillpage')
    
    print(hfig,fname,'-dpng','-r300')
end

% Trajectory figure

if plotnum2 == 1

    hfig = figure;  % save the figure handle in a variable
    fname = '3D Planets Plot';


    picturewidth = 20; % set the width of image in cm
    hw_ratio = .6; % aspect ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust font size


    plot3(sunPosArray(:,1), sunPosArray(:,2), sunPosArray(:,3), '.', 'MarkerSize', 20, 'Color', '#FFA500')
    hold on
    plot3(satPosArray(:,1), satPosArray(:,2), satPosArray(:,3), 'LineWidth', 1, 'Color','#f97306')
    plot3(earthPosArray(:,1), earthPosArray(:,2), earthPosArray(:,3), 'LineWidth', 1.5, 'Color', '#3b3b3b')
    plot3(venusPosArray(:,1), venusPosArray(:,2), venusPosArray(:,3), 'LineWidth', 1.5, 'Color', 'red')
    %plot3(satPosArray(burnIndex,1), satPosArray(burnIndex,2), satPosArray(burnIndex, 3), 'r*')
    title('Earth-Venus-Sun 3D')
    xlabel('Displacement-X')
    ylabel('Displacement-Y')
    zlabel('Displacement-Z')
    legend('Sun', 'Sat', 'Earth', 'Venus')
    hold off
    
    grid on
    
    set(findall(hfig,'-property','Box'),'Box','off') % turn off box
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    


    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-vector','-fillpage')
    
    print(hfig,fname,'-dpng','-r300')


    xl = xlim;
    yl = ylim;
    zl = zlim;
    
    %axis equal
end

% Trajectory figure 2

if plotnum3 == 1

    hfig = figure;  % save the figure handle in a variable
    fname = 'Sat Earth Plot';


    picturewidth = 20; % set the width of image in cm
    hw_ratio = .6; % aspect ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust font size
    plot3(earthPosArray(:,1), earthPosArray(:,2), earthPosArray(:,3), 'LineWidth', 1.5, 'Color', '#3b3b3b')
    hold on
    plot3(satPosArray(:,1), satPosArray(:,2), satPosArray(:,3), 'LineWidth', 1, 'Color','#f97306')
    %plot3(satPosArray(burnIndex,1), satPosArray(burnIndex,2), satPosArray(burnIndex, 3), 'r*')
    title('Earth-Sat 3D')
    xlabel('Displacement-X')
    ylabel('Displacement-Y')
    zlabel('Displacement-Z')
    legend('Earth', 'Sat')
    hold off
    
    grid on
    
    set(findall(hfig,'-property','Box'),'Box','off') % turn off box
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    


    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-vector','-fillpage')
    
    print(hfig,fname,'-dpng','-r300')


    
    %axis equal
end

% Trajectory Animation:
if plotnum4 == 1
    figure(3)
    curve = animatedline('LineWidth', 1.5, 'Color','#069af3');
    curve1 = animatedline('LineWidth', 1, 'Color','#f97306');
    curve2 = animatedline('LineWidth', 1.5, 'Color', '#3b3b3b');
    xlabel('Displacement-X')
    ylabel('Displacement-Y')
    zlabel('Displacement-Z')
    legend('Venus','Earth', 'Sat')
    axis equal
    view(43,24);
    
    output = 0;
    playbackSpeed = 100;

    xlim(xl);
    ylim(yl);
    zlim(zl);



    for i=1:endTime
        index = playbackSpeed * i;

        addpoints(curve, venusPosArray(index,1), venusPosArray(index,2), venusPosArray(index,3));
        addpoints(curve1,earthPosArray(index,1), earthPosArray(index,2), earthPosArray(index,3));
        addpoints(curve2,satPosArray(index,1), satPosArray(index,2), satPosArray(index,3));
        title(sprintf("Satellite Trajectory at time %.2f days", tArray(i) / (3600 * 24) * playbackSpeed));
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
end
