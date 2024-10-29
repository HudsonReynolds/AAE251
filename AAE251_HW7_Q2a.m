% AAE251 Fall 2024
% Homework 7
% AAE251_HW7_Q2a
% Author: Preston Wright and Hudson Reynolds
%% Initializations

fInert = 1/12;              % Inert mass fraction
g = 9.81;                   % Gravity
specImpOne = 295;           % Specific impulse of stage one
specImpTwo = 355;           % Specific impulse of stage two
deltaVTot = 11.85*1000;     % Total delta V required
mPayload = 99;              % Payload mass
tests = 10000;              % Total iterations
i = 1;                      % Iteration counter

%% Calculations

% The initial Gross Lift Off Weight (GLOW) is a function of the delta-v
% proportion for stage 1 and stage 2, alpha 1 and alpha 2, respectively.
% These proportions add to one. To find the ideal value, the value of alpha
% 1 is iterated over in the span from 0.1 to 0.8.

% Initial mass is determined from a variation of the Tsiolkovsky Rocket
% equation for two stages

% Iterate from .3 to .6 for alpha one
for k = .3*tests:tests-.4*tests

    % Find current value of alpha one and two
    alphaOne = k/tests;
    alphaOneArray(i) = alphaOne;
    alphaTwo = 1 - alphaOne;

    % Calculate delta V and the mass fraction
    deltaVOne = alphaOne * deltaVTot;
    deltaVTwo = alphaTwo * deltaVTot;
    mFracOne = exp(deltaVOne/(specImpOne*g));
    mFracTwo = exp(deltaVTwo/(specImpTwo*g));

    % Calculate initial mass
    mInit(i) = mPayload * ((mFracOne*(1-fInert))/(1-fInert*mFracOne)) * ((mFracTwo*(1-fInert))/(1-fInert*mFracTwo));

    % Update array index
    i = i + 1;
end

% Find the location of the minimum initial mass and find the corresponding
% alpha one value
[idealMInit,idealAlphaOneLoc] = min(mInit);
idealAlphaOne = alphaOneArray(idealAlphaOneLoc);

% Calculate corresponding alpha two, final mass, and propellant mass
idealAlphaTwo = 1 - idealAlphaOne;
idealMFinal = idealMInit / mFracOne;
idealMProp = idealMInit - idealMFinal;

%% Graphing

% Plot GLOW vs. Alpha One, with a marker at the minimum value

hfig = figure;  % save the figure handle in a variable

fname = 'HW7a';

title("GLOW vs. Alpha One")
plot(alphaOneArray, mInit, 'Linewidth', 1)
hold on
plot(idealAlphaOne, idealMInit, Marker="x", MarkerSize=15)
xlabel('$\alpha_1$')
ylabel('GLOW [kg]')
xlim([0.3 0.55])
title('Gross Lift Off Weight v. $\alpha_1$')

legend("GLOW", "Minimum GLOW")

grid on



picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.7; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng', '-r300')

%% Outputs

% Output necessary information
fprintf("The alpha one value that leads to a minimum GLOW value is: %.4f\n", idealAlphaOne)
fprintf("The corresponding alpha two value is: %.4f\n", idealAlphaTwo)
fprintf("The corresponding propellant mass is: %.4f\n", idealMProp)



