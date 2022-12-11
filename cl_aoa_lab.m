function cl_aoa_lab(v, opts)
% CL_VS_AOA  Comparison of cl vs aoa between all methods.

%% Parameter settings.

% Set default optional arguments and assert inputs.
if nargin == 0
	v = 16;
end
if nargin <= 1
	opts = 'p';
end
if nargin > 1
	assert( ...
		lab_v == 16 || lab_v == 25, ...
		"Wind tunnel tests have only been performed for 16 m/s and 25 m/s.")
end

% Determine the lab configuration from the desired speed v.
if v == 16
	cfg = 1:3;
elseif v == 25
	cfg = 4:6;
end

% Import the wind tunnel experiment setup.
lab_res = load('group_5.mat');

aoas_lab = lab_res.AoA(cfg);
aoas_num = -5:5:15;

% init cl.
cl_hs    = zeros(size(aoas_num));
cl_lab   = zeros(size(aoas_lab));

%% Lift coefficient computations.

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
% [cl_lab, ~] = arrayfun(@wind_tunnel, cfg(aoas_lab));

% cl from conformal mapping.
tc = 0.18;
cl_cm = @(x) 2*pi * (1 + 4*tc/(3*sqrt(3))) * sin(deg2rad(x));

%% Results plotting.

% Plot these cl.
if contains(opts, 'p')
	% Instantiate the figure.
	figure('WindowStyle','docked');
	hold on;

	% Plot all the methods.
	plot(aoas_num, xfoil_inv(1, :), 'Marker','o');
	if v == 16
		plot(aoas_num, xfoil_visc(1, :), 'LineStyle', '-.', 'Marker','+');
	else
		plot(aoas_num, xfoil_visc(2, :), 'LineStyle', '-.', 'Marker','+');
	end
	plot(aoas_num, cl_hs, 'Marker','x');
	plot(aoas_lab, cl_lab, 'Marker','^');
	fplot(cl_cm, [-5, 16]);

	% Dress the plot.
	title(['Cl vs aoa for v = ', num2str(v), ' m/s']);
	xlabel("aoa (deg)");
	ylabel("cl");
	grid;
	legend( ...
		"Xfoil (inv.)", "Xfoil (visc.)", "Matlab code", "Wind tunnel", "Conformal mapping", ...
		'Location', 'northwest');
end

%% Results registration.

if contains(opts, 'w')
	
end
end
