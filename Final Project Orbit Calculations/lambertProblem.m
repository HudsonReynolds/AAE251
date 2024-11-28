function [v0, vF] = lambertProblem(r0, rF, dT, tm, mu, tolerance, maxSteps, psi, psi_ub, psi_lb)
%% Lambert Problem Function
% function to solve the boundary value problem associated with finding the
% orbit between two points over a given time interval
%
% Inputs:
% r0 - starting point for the transfer, vector [m from sun]
% rF - ending point for the transfer, vector [m from sun]
% dT - time interval over which to transfer [s]
% tm - direction of transfer, 1 = short way, -1 = long way
% mu - gravitational parameter of parent body
% tolerance - tolerance in simulation before stopping
% maxSteps - maximum number of steps before solver times out
% psi - initial guess for the deltaT
% psi_ub - upper bound for the deltaT
% psi_lb - lower bound for the deltaT
% outputs:
% v0 - initial velocity vector
% vF - final velocity vector

% this script is based on the paper by Sharaf et. al.

    %define constants
    sqrtMu = sqrt(mu);
    r0Norm = norm(r0);
    rFNorm = norm(rF);

    %define gamma, the dot product of the initial values normalized (angle)
    gamma = dot(r0, rF) / (r0Norm * rFNorm);

    % define beta
    beta = tm * sqrt(1 - gamma^2);

    %define A
    A = tm * sqrt(r0Norm* rFNorm * (1 + gamma));

    %condition if A equals zero to avoid edge cases
    if A == 0
        A = [0,0,0];
    end

    % make initial guesses for the c2 and c3 values
    c2 = 0.5;
    c3 = 1/6;

    % iterate through values under convergance is achieved
    for i = 1:maxSteps
        B = r0Norm + rFNorm + (A * (psi * c3 - 1)) / sqrt(c2);

        if A > 0 & B < 0
            psi_lb = psi_lb + pi;

            B = B * -1;
        end

        %perform a variable substitution to solve
        chi3 = sqrt(B / c2)^3;

        deltaT = (chi3*c3+A*sqrt(B)) / sqrtMu;

        % break if solution is found
        if abs(dT - deltaT) < tolerance

            fprintf('Converged!\n')

            break
        end

        % modify upper and lower bounds to find solution
        if deltaT <= dT

            psi_lb = psi;
        else
            psi_ub = psi;
        end

        % perform computations for new values of c2 and c3
        psi = (psi_ub + psi_lb) / 2;

        c2 = (1 - cos(sqrt(psi))) / psi;

        c3 = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi));

        f = 1 - B / r0Norm;

        g = A * sqrt(B/mu);

        gdot = 1 - B/rFNorm;

        v0 = (rF - f*r0)/g;
        vF = (gdot*rF-r0)/g;

    end
end
