%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: AAE251 Porkchop
% Author: Hudson Reynolds - Created: 11/27/2024
% Last modified: 11/27/2024
%
% Description: This is the script that performs the calculations for the
% tranfer between planets. This is done by solving Lambert's problem
% between Earth and Venus using a patched conics model. This script outputs
% the delta-V's for each burn and the 

% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization:
% clear console and figures
clear;
clc;
close all;

earthPeriod = 365.256;
venusPeriod = 224.7;

syndoticPeriod = earthPeriod*venusPeriod / (earthPeriod - venusPeriod);


%define the number of days to look over for the departure
departureWindow = round(syndoticPeriod/6);
%departureWindow = 20;

kBound = 20;

lb = 150;

ub = 260;

tEndArray = linspace(lb, ub, kBound);

daysOffset = 360;

parfor i = 1:departureWindow
    %for j = [-1,1]

        for k = 1:kBound

        tStart = datetime('23-03-2025', 'InputFormat','dd-MM-yyyy'); % date at which the simulation starts
        
        tStart = tStart + i + daysOffset;
        
        tStartArray(i) = i + daysOffset;
        
        transferTime = tEndArray(k)
            
        tm = -1;
        
        orbitPlot = 0;
        
        dVtot = TransferCalculator(tStart, transferTime, tm, orbitPlot)

        dVArray(i, k) = dVtot;
        end
    %end
end


%% Plotting

    dVArray(dVArray>20)=20;

    hfig = figure;  % save the figure handle in a variable
    fname = 'Porkchop Plot';

    hold on   
    
    contourf(tStartArray, tEndArray, dVArray', 6, 'ShowText','on')
    
    colormap summer
    
    xlabel('Days after March 25, 2025')
    
    ylabel('Transfer Time [days]')

    title("$\Delta V$ v. Departure and Transit Time")

    ylim([lb ub])
    
    
    picturewidth = 20; % set the width of image in cm
    hw_ratio = 1.2; % aspect ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust font size
        
    set(findall(hfig,'-property','Box'),'Box','off') % turn off box
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-vector','-fillpage')
    print(hfig,fname,'-dpng','-r300')

    %%

    mindV = min(min(dVArray))
    [xIndex,yIndex] = find(dVArray==mindV)

    tStart = datetime('23-03-2025', 'InputFormat','dd-MM-yyyy');

    tStart = tStart + xIndex + daysOffset

    transferTime = round(tEndArray(yIndex))

    [dVtot, dV1, dV2] = TransferCalculator(tStart, transferTime, -1, 1)








