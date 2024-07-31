v_throughScreen = (0.5)*0.3048; % [(ft/s) -> m/s]
q = 0.23; % [m^3/s]
A_screen = q/v_throughScreen;

a = logspace(log10(0.1),log10(1),10); % aspect ratio

R = sqrt(A_screen./(2*pi*a));
h = a.*R;

figure
for i = 1:numel(a)
    rectangle('Position',[-R(i), 0, 2*R(i), h(i)])
    hold on
end