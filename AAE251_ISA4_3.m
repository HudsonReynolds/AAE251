%AAE 251 Fall 2024
%ISA 4
%AAE251_ISA4_3
%Author: Hudson Reynolds

%% Initializations:

radEarth = 6380;       % radius of earth [km]
muEarth = 3.986e5;     % gravitational parameter [km^3 s^-2]
hInit = 220;           % initial orbital altitude [km]
hFinal = 220:100:2020; % final height of orbit [km]
mPay = 4800;           %payload mass

plotVal = 0;           % turn plotting on / off

%% Calculations:
% find the initial orbital height
aInit = radEarth + hInit;
% semi major axis of transfer orbit
aTransfer = radEarth + (hInit + hFinal)/2;
% final semi major axis
aFinal = radEarth + hFinal;
% initial velocity in first orbit
v0 = sqrt(muEarth / aInit);
% velocity at the perigee of the transfer orbit
v1 = sqrt(muEarth * (2 / aInit - 1 ./ aTransfer));
% velocity at the apogee of the transfer orbit
v1_2 = sqrt(muEarth * (2 ./ aFinal - 1 ./ aTransfer));
% velocity of the final orbit
v2 = sqrt(muEarth ./ aFinal);
% delta V of the first burn
dV1 = v1-v0;
% delta V of the second burn
dV2 = v2-v1_2;

%% Delta V's:

%constants:
fInertStage1 = 0.12;
fInertStage2 = 0.145;
ispStage1 = 435;
ispStage2 = 448;
g = 9.81;

% calculations:

%calculate the exhaust velocities:
C1 = ispStage1 * g;
C2 = ispStage2 * g;
% calculate the mass ratios:
MR1 = exp(dV1 * 1000 / C1);
MR2 = exp(dV2 * 1000/ C2);
% calculate the initial mass using the formula with payload mass
mIStage1 = mPay .* MR2 .* (1 - fInertStage2) .* MR2 .* (1 - fInertStage1) ./ ((1 - fInertStage2 .* MR2) .* (1 - fInertStage1 .* MR1));
% calculate the final mass using the mass ratio
mF1 = mIStage1 ./ MR1;
% find the prop mass for first stage
mProp1 = mIStage1 - mF1;

% find prop mass for second stage
mProp2 = mPay .* (MR2-1) .* (1 - fInertStage2) ./ (1 - fInertStage2 * MR2);

%% Plotting
if plotVal == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'Stage 1 Prop Mass v Delta V';

    hold on
    
    plot(hFinal, mProp1, "LineWidth", 1)
    title("Final Orbit Height v. Propellant Mass (Stage 1)")
    xlabel("Final Orbit Height [km]")
    ylabel("Propellant Mass [kg]")

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

% plot 2

    hfig = figure;  % save the figure handle in a variable
    fname = 'Stage 2 Prop Mass v Delta V';

    hold on
    
    plot(hFinal, mProp2, "LineWidth", 1)
    title("Final Orbit Height v. Propellant Mass (Stage 2)")
    xlabel("Final Orbit Height [km]")
    ylabel("Propellant Mass [kg]")

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





