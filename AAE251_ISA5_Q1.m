%AAE 251 Fall 2024
%ISA 5
%AAE251_ISA5_Q1
%Author: Hudson Reynolds

%% Initializations:

plotVal = 0;

% conversions and constants:
ft2m = 0.3048;
g = 9.81;

%variables for aircraft:
mI = 8500;                        % maximum takeoff mass [kg]
mF = 5000;                        %dry mass [kg]
mFuel = 2750;                     % max fuel mass [kg]
crewNum = 2;                      % number of crew
crewMass = 75;                    % mass of crew member [kg]
passNum = 18;                     % max passenger number 
passMass = 80;                    % passenger weight [kg]
bagMass = 20;                     % mass of passenger bag [kg]
cT = 0.58;                        % specific fuel consumtion [lb of fuel / lb thrust-hr]
cT = cT / 3600;                   % convert to 1/s
S = 33.29;                        % wing area of the aircraft [m^2]
V = 201;                          % velocity of flight [m/s]
alt = 35000;                      %altitude [ft]
AR = 7.8;                         % aspect ratio of wing
e = 0.81;                         % oswold efficiency
cD0 = 0.03;                       % zero-lift drag coefficient
k = 1 / (pi * e * AR);            % induced drag coeff.
[~,~,~,rho] = atmosisa(alt*ft2m); % air density [kg/m^3]

% fuel fractions:
W1W0 = 0.995; % fuel fraction for aircraft start
W2W1 = 0.980; % fuel fraction for climb
W4W3 = 0.990; % fuel fraction for descent
W5W4 = 0.987; % fuel fraction for loiter
W6W5 = 0.992; % fuel fraction for landing


%% part A calcs:

mFuelA = mI - crewMass * crewNum - (passMass+bagMass)*passNum - mF;

cruiseStartMass = mI * W1W0 * W2W1;

W3W2 = (1 - (mFuelA / mI)) / (W1W0*W2W1*W4W3*W5W4*W6W5);

cruiseEndMass = mI * W1W0 * W2W1 * W3W2;

cruiseMassAvg = mean([cruiseStartMass, cruiseEndMass]);

cL = 2 * cruiseMassAvg * g / (rho * S * V^2);

cD = cD0 + k * cL^2;

L_D = cL/cD;

rangeA = (-L_D * V * log(W3W2) / cT);

mPayA = mI - crewMass * crewNum - mFuelA - mF;

fprintf("Range at point A %.2f km and payload is %.2f kg\n", rangeA/1e3, mPayA)


%% part B calcs

mFuelB = mFuel;

cruiseStartMass = mI * W1W0 * W2W1;

W3W2 = (1 - (mFuelB / mI)) / (W1W0*W2W1*W4W3*W5W4*W6W5);

cruiseEndMass = mI * W1W0 * W2W1 * W3W2;

cruiseMassAvg = mean([cruiseStartMass, cruiseEndMass]);

cL = 2 * cruiseMassAvg * g / (rho * S * V^2);

cD = cD0 + k * cL^2;

L_D = cL/cD;

rangeB = (-L_D * V * log(W3W2) / cT);

mPayB = mI - crewMass * crewNum - mFuel - mF;

maxPass = floor(mPayB / passMass);

luggageMass = (mPayB - (maxPass * passMass)) / maxPass;

fprintf("Range at point B is %.2f km and payload is %.2f kg\n", rangeB/1e3, mPayB)
fprintf("At point B, %d passengers with %.2f kg of luggage are allowed.\n", maxPass, luggageMass)

%% part C calcs:

mFuelC = mFuel;

mIC = crewMass * crewNum + mFuel + mF;

cruiseStartMass = mIC * W1W0 * W2W1;

W3W2 = (1 - (mFuelC / mIC)) / (W1W0*W2W1*W4W3*W5W4*W6W5);

cruiseEndMass = mIC * W1W0 * W2W1 * W3W2;

cruiseMassAvg = mean([cruiseStartMass, cruiseEndMass]);

cL = 2 * cruiseMassAvg * g / (rho * S * V^2);

cD = cD0 + k * cL^2;

L_D = cL/cD;

rangeC = (-L_D * V * log(W3W2) / cT);

mPayC = mIC - crewMass * crewNum - mFuel - mF;


fprintf("Range at point C is %.2f km and payload is %.2f kg\n", rangeC/1e3, mPayC)

%% Plotting:

ranges = [0, rangeA, rangeB, rangeC] / 1e3;

payloads = [mPayA, mPayA, mPayB, mPayC];

if plotVal == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'Payload Range Diagram';
    colorlist = ["---","#e41a1c", "#377eb8", "#4daf4a"];

    hold on   
    
    for i = 2:4
    plot(ranges(i), payloads(i), 'o', 'Color', colorlist(i)) 
    end
    plot(ranges, payloads, 'k-')
    title("Payload-Range Diagram")
    xlabel("Range [km]")
    ylabel("Payload Weight [kg]")   
    legend('A', 'B', 'C')

    ylim([0, mPayA * 1.1])

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










