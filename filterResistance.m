%% Shelco string wound filters
if 0
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
end

%% Shelco High flow filters

% 5um filters, 60"
flow_ShelcoHFC5u60in = [0
102.857142857143
202.285714285714
300
399.428571428571
497.142857142857
589.714285714285]/15850.3; % [(gpm) -> m^3/s]

pdiff_ShelcoHFC5u60in = [0
0.264350453172205
0.876132930513595
1.63141993957704
2.59818731117825
3.74622356495468
5]*6894.757; % [(psi) -> Pa]

flow_interp = linspace(0,max(flow_ShelcoHFC5u60in), 1e3);
pdiff_interp = interp1(flow_ShelcoHFC5u60in,pdiff_ShelcoHFC5u60in,flow_interp,'spline');

figure

plot(flow_interp,1e-6*pdiff_interp)
hold on
scatter(flow_ShelcoHFC5u60in,1e-6*pdiff_ShelcoHFC5u60in);

%% Shelco High Flow Multi-Cartridge Housings (3, 4, 7, or 12 cartridges)
R_6F = 0.00450839456109*6894.757*15850.3; % [(psi/gpm) -> Pa/(m^3/s)]
R_8F = 0.001718715182745*6894.757*15850.3; % [(psi/gpm) -> Pa/(m^3/s)]
R_10F = 0.001170833178906*6894.757*15850.3; % [(psi/gpm) -> Pa/(m^3/s)]
R_12F = 0.000673426200783*6894.757*15850.3; % [(psi/gpm) -> Pa/(m^3/s)]


%%
rho = 1000;
D = 0.23;
theta_dot = 1; % [rad/s]
q = D*theta_dot; % [m^3/s]
%R_10 = (R(3))*6894.757*15850.3; % [(psi/gpm) -> Pa/(m^3/s)]
d_1 = 0.2; % [m]
A_1 = pi/4*d_1^2;

d_N = 1*2.54/100; % [m]
A_N = pi/4*d_N^2;

% minor losses
k_enteranceFlush = 0.5;
k_exitFlush = 1;

N_cart = [3,4,7,12];
N_housing = 1:5;
dp_filter = zeros(numel(N_cart),numel(N_housing));
dp_housing = zeros(numel(N_cart),numel(N_housing));
dp_minor = zeros(numel(N_cart),numel(N_housing));
dp_total = zeros(numel(N_cart),numel(N_housing));
for i = 1:numel(N_cart)
    for j = 1:numel(N_housing)
        dp_filter(i,j) = interp1(flow_ShelcoHFC5u60in, ...
                                 pdiff_ShelcoHFC5u60in, ...
                                 q/N_housing(j)/N_cart(i),'spline');
        if q/N_housing(j)/N_cart(i) > 500/15850.3
            dp_filter(i,j) = NaN;
        end
        dp_housing(i,j) = R_12F*q/N_housing(j);
        % dp_minor(i,j) = rho/2*(k_exitFlush*(q/A_1)^2 + k_enteranceFlush*(q/N(iN)/A_N)^2);
        dp_total(i,j) = dp_minor(i,j) + dp_filter(i,j) + dp_housing(i,j);
    end
end

%%

figure
subplot(2,1,1)
for j = 1:numel(N_housing) 
    scatter(N_cart,1e-5*dp_total(:,j),100,'LineWidth',2)
    hold on
    leg(j) = {[num2str(N_housing(j)),' housings']};
end
ylabel("pressure loss (bar)")
xlabel('cartridges per housing')
legend(leg)

subplot(2,1,2)
hold on
for j = 1:numel(N_housing) 
    scatter(N_cart.*N_housing(j),1e-5*dp_total(:,j),100,'LineWidth',2)
    hold on
    leg(j) = {[num2str(N_housing(j)),' housings']};
end
ylabel("pressure loss (bar)")
xlabel('total cartridges')
legend(leg)

figure
subplot(2,1,1)
scatter(N_housing,1e-5*dp_housing(1,:),100,'LineWidth',2)
ylabel("pressure loss (bar)")
xlabel('number of housings')
title('Housing Losses')

subplot(2,1,2)
hold on
for j = 1:numel(N_housing) 
    scatter(N_cart.*N_housing(j),1e-5*dp_filter(:,j),100,'LineWidth',2)
    hold on
    leg(j) = {[num2str(N_housing(j)),' housings']};
end
ylabel("pressure loss (bar)")
xlabel('total cartridges')
legend(leg)
title('Filter Losses')