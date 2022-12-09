function cl_vs_aoa(v)
% CL_VS_AOA  Comparison of cl between all approaches, for v16.

if v == 16
	cfg = 1:3;
elseif v == 25
	cfg = 4:6;
else
	error("v must be 16 or 25");
end

% Import the wind tunnel experiment setup.
lab_res = load('group_5.mat');

aoas_lab = lab_res.AoA(cfg);
aoas_num = -5:5:15;

% init cl.
cl_hs    = zeros(size(aoas_num));
cl_lab   = zeros(size(aoas_lab));


% cl Xfoil (inviscid). See screenshots.
xfoil_inv = [-0.6316, 0.0001, 0.6317,  1.2586,  1.8760];

% cl Xfoil (viscous). See screenshots.
xfoil_visc = [ ...
	-0.5213, 0, 0.5214, 1.0992, 1.2486; ...  % v = 16 m/s
	-0.5320, 0, 0.5321, 1.1265, 1.3256];     % v = 25 m/s

% cl hs (our panel code).
for aoa = 1:length(aoas_num)
	[~, ~, cl_hs(aoa), ~] = hess_smith('0018', 200, [aoas_num(aoa), v]);
end

% cl lab.
for aoa = 1:length(aoas_lab)
	[cl_lab(aoa), ~] = wind_tunnel(cfg(aoa));
end

% cl from thin airfoil theory.
cl_tat = @(x) 2*pi*deg2rad(x);

% cl from conformal mapping.
tc = 0.18;
cl_cm = @(x) 2*pi * (1 + 4*tc/(3*sqrt(3))) * sin(deg2rad(x));

% Plot these cl.
figure('WindowStyle','docked');
hold on;
plot(aoas_num, xfoil_inv(1, :), 'Marker','o');
if v == 16
	plot(aoas_num, xfoil_visc(1, :), 'LineStyle', '-.', 'Marker','+');
else
	plot(aoas_num, xfoil_visc(2, :), 'LineStyle', '-.', 'Marker','+');
end
plot(aoas_num, cl_hs, 'Marker','x');
plot(aoas_lab, cl_lab, 'Marker','^');
fplot(cl_tat, [-5, 16]);
fplot(cl_cm, [-5, 16]);
title(['Cl vs aoa for v = ', num2str(v), ' m/s']);
xlabel("aoa (deg)");
ylabel("cl");
grid;
legend("Xfoil(inv)", "Xfoil (visc)", "Our code", "Wind tunnel", "Thin airfoil theory", "Conformal mapping", 'Location','northwest');
