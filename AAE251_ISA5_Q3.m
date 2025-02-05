%AAE 251 Fall 2024
%ISA 4
%AAE251_ISA4_1
%Author: Hudson Reynolds

% aircraft data in following form:

% [fully loaded range [km], cruising mach [-]]
% source link

plotVal = 1;

%% Afterburners:
F22 = [2963, 1.5];
% assuming no drop tanks
%https://www.navy.mil/Resources/Fact-Files/Display-FactFiles/Article/2166130/f-16ab-fighting-falcon-fighter/

Concorde = [6667, 2.02];
% https://www.britishairways.com/content/en/it/information/about-ba/history-and-heritage/celebrating-concorde

%% Turbojets:
Boeing707 = [5200, 0.8];
% https://www.airliners.net/aircraft-data/boeing-707/87

GlosterMeteor = [970, 0.75];
% http://www.rafmuseum.org.uk/research/collections/gloster-meteor-t7/

%% Propeller:
P51D = [1600, 0.47];
% https://www.britannica.com/technology/P-51

Cessna172 = [1290, 0.18];
% https://archive.org/details/janesallworldsai7879fran

%% High Bypass Turbofan
LockheedTristar = [7410, 0.78];
% https://web.archive.org/web/20161115030506/http://rgl.faa.gov/Regulatory_and_Guidance_Library/rgMakeModel.nsf/0/5342e492984e8616862576ac0055feb3/$FILE/a23we.pdf

AirbusA320Neo = [6500, 0.78];
% https://www.airbus.com/aircraftfamilies/passengeraircraft/a320family/technology-and-innovation/

%% Low bypass
A7ECorsair = [1981, 0.9];
% https://www.nationalmuseum.af.mil/factsheets/factsheet?id=298

Convair880 = [4578, 0.76];
% https://books.google.com/books?id=vtsDAAAAMBAJ&pg=PA87#v=onepage&q&f=false

%% Turboprop:
Cessna208Caravan = [1980, 0.279];
% https://cessna.txtav.com/en/turboprop/caravan#_model-specs

PiaggioP180 = [2759, 0.53];
% http://www.avantievo.piaggioaerospace.it


list = vertcat(F22, Concorde, Boeing707, GlosterMeteor, P51D, Cessna172, ...
    LockheedTristar,AirbusA320Neo, A7ECorsair, Convair880, Cessna208Caravan, PiaggioP180);

markerList = ["o", "square", "diamond", "*", "x", "^"];

planeList = ["F-22", "Concorde", "Boeing 707", "Gloster Meteor", "P51D", "Cessna 172", ...
    "Lockheed Tristar", "Airbus A320 Neo", "A7E Corsair", "Convair 880", "Cessna 208 Caravan", "Piaggio P180"];


if plotVal == 1

    hfig = figure;  % save the figure handle in a variable

    fname = 'Aircraft Range-Mach';
    
    picturewidth = 20; % set this parameter and keep it forever
    hw_ratio = 1.2; % feel free to play with this ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust fontsize to your document

    title('Aircraft Range v. Mach for Various Aircraft Types')
        
    hold on
    
    for i = 1:12
        p(i) = plot(list(i,1), list(i,2), 'Marker', markerList(ceil(i/2)), 'Color', 'k');
        if i == 8
            text(list(i,1) + 100, list(i,2) + .03, planeList(i))
        elseif i == 7
             text(list(i,1) + 100, list(i,2) - .02, planeList(i))
        else
            text(list(i,1) + 100, list(i,2), planeList(i))
        end
       
    end
        
    xlabel('Range [km]')
    ylabel('Mach at Cruise [-]')

    xlim([0 9300]);    
    ylim([0 2.1])
    grid minor

    legend([p(1); p(3); p(5); p(7); p(9); p(11)], {'Afterburner', 'Turbojets', 'Propeller', 'High Bypass', 'Low Bypass', 'Turboprop'}, 'Location', 'northwest') 
    
    set(findall(hfig,'-property','Box'),'Box','off') % optional
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-painters','-fillpage')
    print(hfig,fname,'-dpng','-r300')

end






