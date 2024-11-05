%% LiftDragFunc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hudson Reynolds
% Description: function that finds results from the lift and drag equation
% assuming SLUF flight conditions
% and outputs the forces
%
% Inputs:
% A - wing / reference area [m^2]
% rho - density of air [kg/m^3]
% cL0 - coefficient of lift at zero AoA
% cLa - coefficient of lift multiple for 1 degree AoA
% cD0 - coefficient of drag at zero AoA
% cDa - coefficient of drag multiple for induced drag
% V - velocity [m/s]
% W - weight of aircraft [kg]
%
% Outputs:
% pos - initial position vector
% vel - initial velocity vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cL, liftForce, dragForce] = LiftDragFunc(A, rho, cL0, cLa, cD0, cDa, V, W)

% constants:
g = 9.81; % gravity constant [m/s^2]

% calculations

% calculate cL
cL = (2 * W * g) ./ (rho * V.^2 * A);
% calculate cD
cD = cD0 + cDa * (cL).^2;

% find the lift force using cL
liftForce = 1/2 * rho * V.^2 * A .* cL;

% find the drag force using cD
dragForce = 1/2 * rho * V.^2 * A .* cD;

end