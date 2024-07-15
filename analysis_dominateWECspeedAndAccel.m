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
% This script finds the dominate speed and acceleration of a WEC from
% simulation results. Specifically, cases where the velocity and
% acceleration are in the same direction are saught. The intent is for this
% data to feed analysis of the minimum pressure in a WEC-driven pump.
% Velocity causes head losses and reduces static pressure while
% acceleration requires additional pressure differential to accelerate
% fluid having significant inertia.
%
% FILE DEPENDENCY:
% data_parPTO_accum_woPL_wPassiveRV_20240131_02_9090L.mat
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

% load data from representative simulation of WEC
load("data_parPTO_accum_woPL_wPassiveRV_20240131_02_9090L.mat")
theta_dot = out.dydt(:,par.iy.theta);
theta_ddot = out.dydt(:,par.iy.theta_dot);

% retain only velocity and acceleration in the same direction
iVec = find(theta_ddot.*theta_dot > 0);

% thin data to accepatible size
nLim = 10000; % target limit for number of data points
r = ceil(numel(iVec)/nLim);% factor for thining timeseries data
iVec = iVec(1:r:numel(iVec));
theta_dot = theta_dot(iVec);
theta_ddot = theta_ddot(iVec);

% determine diminate data points
iVec_dom = paretoFront2D(abs(theta_dot),'max', ...
                         abs(theta_ddot),'max');
theta_dot_dom = abs(theta_dot(iVec_dom));
theta_ddot_dom = abs(theta_ddot(iVec_dom));

% plot dominate data points
figure
scatter(theta_dot_dom,theta_ddot_dom)
xlabel('angular speed (rad/s)')
ylabel('angular acceleration (rad/s^2)')
title('Dominate Simultaneous WEC Speed and Acceleration')