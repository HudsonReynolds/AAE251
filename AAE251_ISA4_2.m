%AAE 251 Fall 2024
%ISA 4
%AAE251_ISA4_1
%Author: Hudson Reynolds

%% Initializations:

radEarth = 6380;       % radius of earth [km]
muEarth = 3.986e5;     % gravitational parameter [km^3 s^-2]
hInit = 220;           % initial orbital altitude [km]
hFinal = 220:100:2020; % final height of orbit [km]
plotVal = 1;           % turn plotting on / off

aInit = radEarth + hInit;

aTransfer = radEarth + (hInit + hFinal)/2;

aFinal = radEarth + hFinal;

v0 = sqrt(muEarth / aInit);

v1 = sqrt(muEarth * (2 / aInit - 1 ./ aTransfer));

v1_2 = sqrt(muEarth * (2 ./ aFinal - 1 ./ aTransfer));

v2 = sqrt(muEarth ./ aFinal);

dV1 = v1-v0;

dV2 = v2-v1_2;

if plotVal == 1
    close all
    
    hfig = figure;  % save the figure handle in a variable
    fname = 'hf v Delta V';

    hold on   
    
    plot(hFinal, dV1, "LineWidth", 1)
    plot(hFinal, dV2, "LineWidth", 1)
    title("Final Orbit Height v. $\Delta V_1$ and $\Delta V_2$")
    xlabel("Final Orbit Height [km]")
    ylabel("$\Delta V$ [km/s]")
    legend('$\Delta V_1$', '$\Delta V_2$', 'Location','northwest')

    
    
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
