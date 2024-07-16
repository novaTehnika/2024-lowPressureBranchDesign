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
% Fluid properties
par.rho = 1023; % [kg/m3] density of air
par.mu = 9.4e-4; % [Pa-s]  Dynamic (absolute) viscosity
par.beta = 2.2e9; % [Pa]  Bulk Modulus of air free fluid
par.p_vap = 0.037e5; % [Pa] vapour pressure of seawater
par.R = 0.0001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
par.p_o = 101.3e3; % [Pa]  Atmospheric pressure (reference)
par.gamma = 1.4; % [-]ratio of specific heats for air

% WEC-driven pump
D = 0.23; % pump displacement

 % intake check valve of WEC-driven pump
kv_cv = 0.00048;
pc_cv = 1e4;
dp_cv = 1e4;

% Piping
 % intermediate line between WEC-driven pump port and check valve
L_intermediate = 5; % [m]
 % intake pipeline
L_intake = 0; % [m]

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
d_pipe = linspace(0.05,0.2,16);

ni = numel(theta_dot_dom);
nj = numel(d_pipe);
deltaP_major = zeros(ni,nj);
deltaP_minor = zeros(ni,nj);
deltaP_cv = zeros(ni,nj);
deltaP_accel = zeros(ni,nj);
deltaP_total = zeros(ni,nj);

for i = 1:ni
    q = D*theta_dot_dom(i);
    dq_dt = D*theta_ddot_dom(i);
    for j = 1:nj
        A_pipe = pi/4*d_pipe(j)^2;
        
        deltaP_major(i,j) = (flowR(q,d_pipe(j),L_intermediate,par) + ...
                             flowR(q,d_pipe(j),L_intake,par))*q;
        
        deltaP_minor(i,j) = k_total*par.rho/2*(q/A_pipe)^2;
        
        deltaP_cv(i,j) = max(dp_cv/(kv_cv*sqrt(pc_cv+dp_cv)*q), ...
                             (q/kv_cv)^2);

        I_intermediate = par.rho*L_intermediate/A_pipe;
        I_intake = par.rho*L_intake/A_pipe;
        deltaP_accel(i,j) = (I_intermediate+I_intake)*dq_dt;

        deltaP_total(i,j) = deltaP_major(i,j) + deltaP_minor(i,j) + ...
                            deltaP_cv(i,j) + deltaP_accel(i,j);
    end
end

%% acceleration head at peak acceleration
[~, i] = max(theta_ddot_dom);
figure
plot(d_pipe*100,1e-5*deltaP_accel(i,:))
xlabel('pipe diameter (cm)')
ylabel('pressure drop (bar)')
title('Acceleration Pressure Drop at Peak Acceleration')

%% Major and minor losses at peak velocity
[~, i] = max(theta_dot_dom);
figure
plot(d_pipe*100,1e-5*deltaP_major(i,:))
hold on
plot(d_pipe*100,1e-5*deltaP_minor(i,:))
plot(d_pipe*100,1e-5*deltaP_cv(i,:))
plot(d_pipe*100,1e-5*(deltaP_major(i,:) + ...
                      deltaP_minor(i,:) + deltaP_cv(i,:)))
legend('major','minor','check valve','combined')
xlabel('pipe diameter (cm)')
ylabel('pressure drop (bar)')
title('Flow Rate Dependent Pressure Loss at Peak Velocity')


%% total pressure drop for all dominate speed and accel
figure
for i = 1:ni
plot(d_pipe*100,1e-5*deltaP_total(i,:))
hold on
end
xlabel('pipe diameter (cm)')
ylabel('pressure drop (bar)')
title('Total Pressure Drop')

%% total pressure drop, worst case as func of pipe ID
max_deltaP = max(deltaP_total,[],1);
figure
plot(d_pipe*100,1e-5*max_deltaP)
xlabel('pipe diameter (cm)')
ylabel('pressure drop (bar)')
title('Maximum Total Pressure Drop')