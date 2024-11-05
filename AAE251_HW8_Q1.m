A = 88.2;
[~, ~, ~, rho] = atmosisa(0);
W = 33100;
cL0 = 0.02;
cLa = 0.12;
cD0 = 0.015;
cDa = 0.05;
V = linspace(100,275, 176);

[~, lift, drag] = LiftDragFunc(A, rho, cL0, cLa, cD0, cDa, V, W);

