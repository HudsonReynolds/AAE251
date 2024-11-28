%% Acceleration Function
% function describing the translational accelerations on the satellite
% inputs:
% t- input time [s]
% v - current velocity [-/s]
% r - distance from bodies
% burnStart - start time of burn [s]
% burnDuration - duration of burn [s]
% dV - total delta V of the burn [-/s]
% mu - gravitational constant of bodies (must be in same order as r)
% burn - 0 => no burn compute, 1 => burn compute
% outputs:
% dv - acceleration vector for integration
function dv = accel(t, v, r, mu)
    % gravity forces:
    gravAccel = 0;
    for i = 1:length(mu)
        r0 = (i - 1) * 3 + 1;
        rf =  i * 3;
        gravAccel = gravAccel + (mu(i) / (norm(r(r0:rf))^2) * (r(r0:rf) / norm(r(r0:rf))));
    end
    
    dv = gravAccel;
end