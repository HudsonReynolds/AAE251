function meanLong = meanLongSolver(t1, t2, n, L0)
%% Mean Longitude Function
% function to solve the boundary value problem associated with finding the
% orbit between two points over a given time interval
%
% Inputs:
% t1- starting time at the epoch [date]
% t2 - end time [date]
% n - mean motion [rad/s]
% L0 - mean longitude at the epoch [deg]
%
% Outputs:
% meanLong - mean longitude at t2

timediff = abs(t1 - t2);

timediff = seconds(timediff);

L0rad = deg2rad(L0);

meanLong = L0rad + n * timediff;

meanLong = mod(rad2deg(meanLong), 360);


