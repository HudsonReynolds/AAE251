function [thrust, thrustReserve] = ThrustRequiredJetFunc(V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hudson Reynolds
% Description: function that finds results from the lift and drag equation
% assuming SLUF flight conditions
% and outputs the forces
%
% Inputs:
% V - velocity [m/s]
%
% Outputs:
% thrust - the required thrust to maintain SLUF conditions [N]
% thrustReserve - the percentage of thrust remaining [N]
% plots - see outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 88.2;
[~, ~, ~, rho] = atmosisa(0);
W = 33100;
cL0 = 0.02;
cLa = 0.12;
cD0 = 0.015;
cDa = 0.05;


[~, lift, drag] = LiftDragFunc(A, rho, cL0, cLa, cD0, cDa, V, W);


hfig = figure;  % save the figure handle in a variable
hold on

fname = 'Lab 5 Shear Modulus Virt';

plot(V, drag / 1e3)
title("Velocity v. Thrust")
xlabel("Velocity [m/s]")
ylabel("Thrust [kN]")



picturewidth = 20; % set this parameter and keep it forever
hw_ratio = .6; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust fontsize to your document


grid on

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,fname,'-dpdf','-vector','-fillpage')

%print(hfig,fname,'-dpng','-r300')




