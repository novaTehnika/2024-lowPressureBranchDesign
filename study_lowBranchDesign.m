% analysis_dominateWECspeedAndAccel.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 07/15/2024
%
% PURPOSE/DESCRIPTION:
% This script performs a study on total head losses for a realistic intake
% circuit for a WEC-driven pump. Effects of check valves, intake filters,
% fluid inertia, and minor losses are taken into account. The study varies
% the nominal pipe diameter and includes realistic combinations of WEC
% speed and acceleration. Data for the WEC speed and acceleration are taken
% from the file "data_dominateVelAccel.mat", which is derived from the
% analysis performed by the script "analysis_dominateWECspeedAndAccel.m".
%
% FILE DEPENDENCY:
% data_dominateVelAccel.mat
%
% UPDATES:
% 07/15/2024 - Created.
%
% Copyright (C) 2024  Jeremy W. Simmons II
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load date for peak simultanrous WEC speed and velocity
load("data_dominateVelAccel.mat")

%% define system parameters
% WEC-driven pump
D = 0.23; % pump displacement

 % intake check valve of WEC-driven pump
kv_cv = 0.00048;
pc_cv = 1e5;
dp_cv = 1e5;

% Minor losses in the circuit
 % Long radius, 90-degree bends with flanged fittings
k_90degLongRadFlanged = 0.2; % loss coeff.
n_90degLongRadFlanged = 2; % number of elements of this type

 % Tee fitting with flow making a 90-degree turn (branched flow)
k_teeBranchedFlow = 1.0;
n_teeBranchedFlow = 1;

 % Intake modeled as flow from a resivoiur into a pipe with a well radiused
 % transition to the nominal pipe diameter
k_intake = 0.04;
n_intake = 1;

 % total loss coeff. for minor losses
k_total = k_90degLongRadFlanged*n_90degLongRadFlanged + ...
            k_teeBranchedFlow*n_teeBranchedFlow + ...
            k_intake*n_intake;

% resistance of intake filter, modeled as linear resitance
R_filter = 1e5/0.1; % [Pa/(m^3/s)]

%% study the effect of combined WEC speed and acceleration and sizing of
% the pipe
d_pipe = linspace(0.05,0.1,5);
for i = 1:numel(theta_dot_dom)
    for j = 1:numel(d_pipe)
        deltaP_major(i,j) = 0;%[darcyWiesbachEQ];
        deltaP_minor(i,j) = 0;%[minorLossEQ];
        deltaP_cv(i,j) = max(dp_cv/(kv_cv*sqrt(pc_cv+dp_cv)*q),(q/kv_cv)^2);
        deltaP_accel(i,j) = I*D*theta_ddot_dom(i);
    end
end