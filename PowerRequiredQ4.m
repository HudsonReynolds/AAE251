function [power, powerReserve, maxV] = PowerRequiredQ4(V, height, plotVal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hudson Reynolds, Preston Wright
% Description: function that finds power for prop aircraft based on the
% velocity
%
% Inputs:
% V - velocity [m/s]
%
% Outputs:
% thrust - the required thrust to maintain SLUF conditions [N]
% thrustReserve - the percentage of thrust remaining [N]
% plots - see outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

height = 4618.0248;
V = 50:1:175;
plotVal = 1;

A = 39.4; %wing area [m^2]
b = 19.78;
[~, ~, ~, rho0] = atmosisa(0); % density of air at sea level [kg/m^3]
[~, ~, ~, rho] = atmosisa(height); % density of air [kg/m^3]
W = 10500;      % weight [kg]
cD0 = 0.021;   % zero AoA cD
p0max = 1342.26 * 2 * 1000;   % sea level power [kW]
eta = 0.8;     % propeller efficiency 
e = 0.7;
AR = b^2 / A

cDa = 1 / (pi * e * AR);   % induced drag coefficient

power = 1/2 * rho * A * V.^3 * cD0 + 2 * (W * 9.81)^2 ./ (e * AR * pi * rho * A * V);

powerMax = eta * (rho / rho0)^0.6 * p0max;

powerReserve = 1 - (power / powerMax);

[~, minIndex] = min(abs(power - powerMax));

maxV = V(minIndex)

syms x

eqn = powerMax == 1/2 * rho * A * x^3 * cD0 + 2 * (W * 9.81)^2 / (e * AR * pi * rho * A * x);

sol = solve(eqn, x, real=true)

fprintf("%.2f\n %.2f\n", sol(1), sol(2))


if plotVal == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'Power v Velocity Graph Q4';

    hold on   
    
    plot(V, power / 1e3)
    yline(powerMax/ 1e3 ,'r--')
    title("Velocity v. Power at 15,151 ft")
    xlabel("Velocity [m/s]")
    ylabel("Power [kW]")
    

    
    
    picturewidth = 20; % set the width of image in cm
    hw_ratio = .6; % aspect ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust font size
    
    legend('Power Required', 'Max Power', 'FontSize', 12, 'location', 'northwest')
    
    
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