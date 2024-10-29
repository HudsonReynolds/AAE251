% AAE251 Fall 2024
% Homework 7
% AAE251_HW7_Q2a
% Author: Preston Wright and Hudson Reynolds

%% Initializations 

isp2List = [356, 372, 383, 399, 445]; % list of possible isps [s]
alphaList = 0.30:0.001:0.5;           % list of alpha values
fInert = 1/12;                        % Inert mass fraction
g = 9.81;                             % Gravity
isp1 = 295;                           % Specific impulse of stage one
deltaVTot = 11.85*1000;               % Total delta V required
mPayload = 99;                        % Payload mass
i = 1;                                % Iteration counter

%% Calculations

% The initial Gross Lift Off Weight (GLOW) is a function of the delta-v
% proportion for stage 1 and stage 2, alpha 1 and alpha 2, respectively.
% These proportions add to one. To find the ideal value, the value of alpha
% 1 is iterated over in the span from 0.5 to 0.67.

% Initial mass is determined from a variation of the Tsiolkovsky Rocket
% equation for two stages

for k = 1:length(alphaList)

    % Find current value of alpha one and two
    alpha1 = alphaList(k);
    alpha2 = 1 - alpha1;
    alpha2Array(k) = alpha2;

    % Calculate delta V and the mass fraction
    deltaV1 = alpha1 * deltaVTot;
    deltaV2 = alpha2 * deltaVTot;
    mFracOne = exp(deltaV1/(isp1*g));
    for j = 1:length(isp2List)
        mFracTwo(j) = exp(deltaV2/(isp2List(j)*g));
        % Calculate initial mass
        mInit(k, j) = mPayload * ((mFracOne*(1-fInert))/(1-fInert*mFracOne)) * ((mFracTwo(j)*(1-fInert))/(1-fInert*mFracTwo(j)));
    end

    % Update array index
    %i = i + 1;
end

hfig = figure;  % save the figure handle in a variable

fname = 'HW7b';


plot(alpha2Array, mInit)
title("GLOW vs. $\alpha_2$")

xlabel('$\alpha_2$')
ylabel('GLOW [kg]')

xlim([0.5 0.67])

legend("356 s", "372 s", "383 s", "399 s", "445 s", 'location', 'northwest')

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

