function [thrust, thrustReserve] = ThrustRequiredJetFunc(V, height, plotVal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hudson Reynolds, Preston Wright
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

% included so script doesn't throw errors when publishing. Delete these to
% run it as a function
height = 0;
V = 100:1:275;
plotVal = 1;

% constants:
[~, ~, ~, rho0] = atmosisa(0);     % density of air at sea level [kg/m^3]
[~, ~, ~, rho] = atmosisa(height); % density of air [kg/m^3]
A = 88.2;                          %wing area [m^2]
W = 33100;                         % weight [kg]
cL0 = 0.02;                        % zero AoA cL
cLa = 0.12;                        % slope of cL / alpha 
cD0 = 0.015;                       % zero AoA cD
cDa = 0.05;                        % induced drag coefficient
t0max = 55620;                     % sea level thrust [N]

% calculations:
[~, lift, drag] = LiftDragFunc(A, rho, cL0, cLa, cD0, cDa, V, W);

thrust = drag;

thrustMax = (rho / rho0)^0.6 * t0max;

thrustReserve = 1 - (thrust / thrustMax);


%plots:
if plotVal == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'Thrust v Velocity Graph Jet';

    hold on
       
    plot(V, thrust / 1e3)
    title("Velocity v Thrust Jet Aircraft")
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
    %print(hfig,fname,'-dpdf','-vector','-fillpage')
    
    print(hfig,fname,'-dpng','-r300')
end




