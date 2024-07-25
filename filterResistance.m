%% Shelco string wound filters

filtSize = [3
5
10
20
30
50
75];

R = [0.998840622978554
0.507517074676804
0.330591973978289
0.228994770518389
0.182111279525073
0.149136500173677
0.122774616726664];

filtSizeInterp = linspace(0.5,100,1e3);
Rinterp_linear = interp1(filtSize,R,filtSizeInterp,'linear');
Rinterp_makima = interp1(filtSize,R,filtSizeInterp,'makima');
Rinterp_linearExp = exp(interp1(log(filtSize),log(R),log(filtSizeInterp),'linear'));
Rinterp_makimaExp = exp(interp1(log(filtSize),log(R),log(filtSizeInterp),'makima'));

figure
% loglog(filtSizeInterp,Rinterp_linear)
% hold on
% loglog(filtSizeInterp,Rinterp_makima)
loglog(filtSizeInterp,Rinterp_linearExp)
hold on
loglog(filtSizeInterp,Rinterp_makimaExp)
scatter(filtSize,R);

legend('linear','makima','linear, exp','makima, exp','data')


%% Shelco High flow filters

flow_ShelcoHFC5u = [0
102.857142857143
202.285714285714
300
399.428571428571
497.142857142857
589.714285714285]/15850.3; % [(gpm) -> m^3/s]

pdiff_ShelcoHFC5u = [0
0.264350453172205
0.876132930513595
1.63141993957704
2.59818731117825
3.74622356495468
5]*6894.757; % [(psi) -> Pa]

flow_interp = linspace(0, 600, 1e3);
pdiff_interp = interp1(flow_ShelcoHFC5u,pdiff_ShelcoHFC5u,flow_interp,'spline')

figure

plot(flow_interp,pdiff_interp)
hold on
scatter(flow_gpm,pdiff_ShelcoHFC5u);

%%
rho = 1000;
D = 0.23;
theta_dot = 1; % [rad/s]
q = D*theta_dot; % [m^3/s]
R_10 = (R(3))*6894.757*15850.3; % [(psi/gpm) -> Pa/(m^3/s)]
d_1 = 0.2; % [m]
A_1 = pi/4*d_1^2;

d_N = 1*2.54/100; % [m]
A_N = pi/4*d_N^2;

% minor losses
k_enteranceFlush = 0.5;
k_exitFlush = 1;

N = 1:20;
dp_filter = zeros(numel(N),1);
dp_minor = zeros(numel(N),1);
dp_total = zeros(numel(N),1);
for iN = 1:numel(N)
    dp_filter(iN) = R_10*q/N(iN);
    dp_minor(iN) = rho/2*(k_exitFlush*(q/A_1)^2 + k_enteranceFlush*(q/N(iN)/A_N)^2);
    dp_total(iN) = dp_minor(iN) + dp_filter(iN);
end

figure
semilogy([],[])
hold on
scatter(N,1e-6*dp_filter,100,'kx','LineWidth',2)
scatter(N,1e-6*dp_minor,100,'ks','LineWidth',2)
scatter(N,1e-6*dp_total,100,'r^','LineWidth',2)
