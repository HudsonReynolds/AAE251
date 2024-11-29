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

tStart = datetime('23-03-2025', 'InputFormat','dd-MM-yyyy'); % date at which the simulation starts

transferTime = 147;

tm = 1;

orbitPlot = 1;

TransferCalculator(tStart, transferTime, tm, orbitPlot)

