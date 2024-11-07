function [power, powerReserve] = PowerRequiredPropFunc(V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hudson Reynolds
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

height = 0;
V = 50:1:175;
plotVal = 1;

A = 16.3; %wing area [m^2]
[~, ~, ~, rho0] = atmosisa(0); % density of air at sea level [kg/m^3]
[~, ~, ~, rho] = atmosisa(height); % density of air [kg/m^3]
W = 1315;      % weight [kg]
cL0 = 0.02;    % zero AoA cL
cLa = 0.12;    % slope of cL / alpha
cD0 = 0.026;   % zero AoA cD
cDa = 0.054;   % induced drag coefficient
p0max = 216;   % sea level power [kW]
eta = 0.8;     % propeller efficiency 


[~, lift, drag] = LiftDragFunc(A, rho, cL0, cLa, cD0, cDa, V, W);

power = drag .* V;

powerMax = eta * (rho / rho0)^0.6 * p0max;

powerReserve = 1 - (power / powerMax);


if plotVal == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'Power v. Velocity Graph';

    hold on   
    
    plot(V, power / 1e3)
    title("Velocity v. Power")
    xlabel("Velocity [m/s]")
    ylabel("Power [kW]")
    
    
    
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