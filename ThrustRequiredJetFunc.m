%function [thrust, thrustReserve] = ThrustRequiredJetFunc(V, height, plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hudson Reynolds
% Description: function that finds thrust for jet based on velocity
%
% Inputs:
% V - velocity [m/s]
% height - current altitude of jet [m]
% plot - turn plotting on or off. Set to 1 to plot.
%
% Outputs:
% thrust - the required thrust to maintain SLUF conditions [N]
% thrustReserve - the percentage of thrust remaining [N]
% plots - see outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants:

height = 0
V = linspace(100,275,176)
plot = 1

A = 16.3; %wing area [m^2]
[~, ~, ~, rho0] = atmosisa(0); % density of air at sea level [kg/m^3]
[~, ~, ~, rho] = atmosisa(height); % density of air [kg/m^3]
W = 1315; % weight [kg]
cL0 = 0.02; % zero AoA cL
cLa = 0.12;  % slope of cL / alpha
cD0 = 0.026; % zero AoA cD
cDa = 0.054; % induced drag coefficient
t0max = 55620; % sea level thrust [N]


[~, lift, drag] = LiftDragFunc(A, rho, cL0, cLa, cD0, cDa, V, W);

thrust = drag;

thrustMax = (rho / rho0)^0.6 * t0max;

thrustReserve = 1 - thrust / thrustMax;


if plot == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable

    hold on
    
    fname = 'Thrust v. Velocity Graph';
    
    plot(V, thrust)
    title("Velocity v. Thrust")
    xlabel("Velocity [m/s]")
    ylabel("Thrust [kN]")
    
    
    
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
    print(hfig,fname,'-dpdf','-vector','-fillpage')
    
    %print(hfig,fname,'-dpng','-r300')
end




