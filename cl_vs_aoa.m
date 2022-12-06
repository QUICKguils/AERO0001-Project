% CL_VS_AOA  Comparison of cl between all approaches, for v16.

clear; close all;

% Import the wind tunnel experiment setup.
lab_set = load('setup.mat');
lab_res = load('group_5.mat');
% Choose a configuration (aoa and speed, see group_6 data).
cfg = 1:3;

aoas = lab_res.AoA(cfg);

% init cl.
cl_hs    = zeros(1, 3);
cl_xfoil = zeros(1, 3);
cl_lab   = zeros(1, 3);

% cl xfoil (viscous). See screenshots.
cl_xfoil = [0.5214, 1.0992, 1.2486];

% cl hs (our panel code).
for aoa = 1:length(aoas)
	[~, ~, cl_hs(aoa), ~] = hess_smith('0018', 200, cfg(aoa));
end

% cl lab.
for aoa = 1:length(aoas)
	% TODO implement the cl calculations.
end

% cl from thin airfoil theory.
cl_tat = @(x) 2*pi*deg2rad(x);


% Plot these cl.
figure('WindowStyle','docked');
hold on;

plot(aoas, [cl_xfoil; cl_hs]);
fplot(cl_tat, [-5, 18]);
title("Cl vs aoa.");
xlabel("aoa (deg)");
ylabel("Cl");
grid;
legend("Xfoil(viscous)", "Hess and Smith", "Thin airfoil theory", 'Location','northwest');