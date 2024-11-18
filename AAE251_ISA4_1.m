%AAE 251 Fall 2024
%ISA 4
%AAE251_ISA4_1
%Author: Hudson Reynolds

% close all currently open
clc,clear, close all

%% Initializations:
OmegaEarth = 7.292e-5; % sidereal rate [rad s^-1]
radEarth = 6380;       % radius of earth [km]
muEarth = 3.986e5;     % gravitational parameter [km^3 s^-2]
dVLoss = 1.7;          % delta-V loss [km s^-1]
radOrbit = 220;        % orbital altitude [km]
latKSC = [28,31,27];   % latitude [degree, min, sec]
Az= 90;                % azimuth of launch [deg]
isp = 420;             % isp of rocket [s]
mPay = 4800;           % payload mass [kg]
g =9.807;              % gravitational acceleration [m/s^2]
plotVal = 0;           % turn plotting on or off


%% Calculations:

%start by finding the total delta V requirements of the vehicle:

latKSC = latKSC(1) + (latKSC(2)/60) + (latKSC(3)/3600);

dV_EarthHelp = OmegaEarth * radEarth * cosd(latKSC) * sind(Az);

a = radOrbit + radEarth;

dVLEO = sqrt(muEarth / a);

orbitalInc = acosd(sind(Az) * cosd(latKSC));

nu = abs(orbitalInc - 28);
dVIncChange = 2 * dVLEO * sind(nu/2);

dVtot = dVLEO + dVLoss + dVIncChange - dV_EarthHelp

% find the inert mass fraction

fInert = linspace(0.05, 0.25, 100);

C = isp * g;

massRatio = exp(dVtot*1000/C);

mProp = mPay * (massRatio - 1) * (1 - fInert) ./ (1 - fInert .* massRatio);

mInert = -fInert .* mProp ./ (fInert - 1);

mI = mPay + mProp + mInert;


if plotVal == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'fInert v Mass Plot';

    hold on   
    
    plot(fInert, mI, "LineWidth", 1)
    xline(1/massRatio, 'r--', 'LineWidth', 1)
    title("GLOW v. $f_{inert}$")
    xlabel("$f_{inert}$")
    ylabel("GLOW [kg]")
    legend('GLOW', 'Singularity')
    
    
    
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




