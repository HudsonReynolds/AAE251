%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE251 Fall 2024
% Homework 9
% AAE251_HW9_Q1cd
% Author: Preston Wright and Hudson Reynolds
% Description: Sets up and calculates the available and required thrust with
% respect to altitude for a given jet aircraft, plotting those values versus
% the altitude used to calculate them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

alt = linspace(0,30,30*100+1);   % Altitude array [km]
thrustAvailable = [];            % Initialized available thrust array [N]
thrustRequired = [];             % Initialized required thrust array [N]

rhoSea = 1.2250;        % Density at sea level [kg/m^3]
g = 9.81;               % Gravitational acceleration [m/s^2]
m = 33100;              % Mass of the aircraft [kg]
mAD = 0.6;              % Air density exponent
thrustMax = 55620;      % Maximum available thrust at sea level [N]
K = 0.05;               % Wingspan efficiency
parasiteDrag = 0.015;   % Parasitic drag


%% Calculations

% Loop through the altitudes calculating the available and required thrust
for i = 1:length(alt)
    [~,~,~,rhoAlt] = atmosisa(alt(i)*1000, extended="on");
    thrustAvailable(i) = ((rhoAlt/rhoSea)^mAD)*thrustMax;
    thrustRequired(i) = 2*m*g*sqrt(K*parasiteDrag);
end

%% Graphing

% Output the required and available thrust with respect to altitude

hfig = figure;  % save the figure handle in a variable


fname = "Available and Required Thrust Vs Altitude";

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = .7; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust fontsize to your document

hold on

plot(alt,thrustRequired, 'linewidth', 1)
hold on
plot(alt,thrustAvailable, 'LineWidth',1)
grid minor
title("Available and Required Thrust Vs. Altitude")
xlabel("Altitude [km]")
ylabel("Thrust [N]")
legend("Required Thrust", "Available Thrust", location="northeast")


set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-r300')




