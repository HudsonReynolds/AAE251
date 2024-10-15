%% 1a

rEarth = 6380; % km
muEarth = 3.986004418e5; %
rP = rEarth + 185;
rA = rEarth + 35786;
a = (rP + rA) / 2;

e = (rA / a) - 1;

P = 2 * pi * sqrt(a^3 / muEarth);

vP = sqrt(muEarth * (2 / rP - 1/a));

vA = sqrt(muEarth * (2 / rA - 1/a));

%% 1b

% inc change then circularize

dVi1 = 2 * vA * sind(27/2);

vCirc = sqrt(muEarth / rA);

dVcirc = vCirc - vA;

dVtot1 = dVi1 + dVcirc;

% inc change then circularize

dVi2 = 2 * vCirc * sind(27/2);

dVtot2 = dVi2 + dVcirc;

%% 1c

dVBudget = 3; %km/s

dVleftover1 = (dVBudget - dVtot1) * 1000

dVleftover2 = (dVBudget - dVtot2) * 1000

t1 = dVleftover1 / 45;
t2 = dVleftover2 / 45;

profit1 = (t1 - 15)* 50000

%% 1d
rP2 = rEarth + 250;
rA2 = rEarth + 60800
a2 = (rP2 + rA2) / 2;

vP2 = sqrt(muEarth * (2 / rP2 - 1/a2));
vA2 = sqrt(muEarth * (2 / rA2 - 1/a2));

rGeo = rEarth + 35786;

dVi3 = 2 * vA2 * sind(27/2);

aManuv1 = (rA2 + rGeo) / 2;

dVPerigeeRaise = sqrt(muEarth * ((2 / rA2) - (1 / aManuv1))) - sqrt(muEarth * ((2 / rA2) - (1 / a2)));

dVCirc2 = abs(sqrt(muEarth / rGeo) - sqrt(muEarth * ((2 / rGeo) - (1 / aManuv1))));

dVtot = dVi3 + dVPerigeeRaise + dVCirc2;

dVleftover3 = (3 - dVtot) * 1000;

t3 = dVleftover3 / 45;

profit2 = (t3 - 15) * 50000;

netProfit = profit2 - profit1;

fprintf("net profit is $%.6g \n", netProfit)
